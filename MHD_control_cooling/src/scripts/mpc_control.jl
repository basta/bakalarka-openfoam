#=
Main script to run Dense MPC simulation using DenseMPCLib.jl
Infers system dimensions N and m from loaded A, B matrices.
Requires m_phys (physical input dimension) to be set manually.
Includes initialization from historical data.
=#

using Pkg
Pkg.activate(".") # Activate the project environment relative to this script

# Add necessary packages if not already present
# Pkg.add(["Revise", "LinearAlgebra", "GLMakie", "Logging", "JuMP", "OSQP", "BlockArrays", "SparseArrays", "Statistics", "JLD2"])

using Revise
using LinearAlgebra
using GLMakie # Ensure GLMakie is available for plotting
using Logging
using Statistics
using JLD2
using SparseArrays
using BlockArrays
using JuMP

# --- Include the library module ---
# Assumes mpc_lib.jl is in the parent directory relative to this script's location
# Adjust the path if necessary
try
    includet("../mpc_lib.jl") # Changed path assumption
    using .DenseMPCLib # Bring exported names into scope
catch e
    error("Failed to include mpc_lib.jl. Ensure the file exists at the specified path. Error: $e")
end

includet("../dmd_lib.jl")

# --- 1. User-Defined Parameters ---
const T = Float64 # Set numeric type
const history_path = "./data/dataset-live2.jld2"

# --- Physical System Parameters (MUST BE SET MANUALLY) ---
const m_phys = 16

# --- MPC Problem Definition ---
const Np = 5              # Prediction horizon
const Q_state_val = 1.0   # Cost weight on original state x_k (assumed to be z_1)Q
const R_phys_val = 0.05   # Cost weight on physical input u_phys_k
const Q_term_val = 10.0   # Terminal cost weight on original state x_k (z_1)
const u_phys_max = 2.0    # Physical input constraint: |u_phys_j| <= u_phys_max
const target_state_val = 0.0 # Target value for the original state x_k (z_1)

# --- Simulation Setup ---
const sim_steps = 1000
# const start_state_val = 2000 # Initial value if NOT using history - REMOVED

# --- History Initialization Parameters (NEW) ---
const ny_init = 200 # ASSUMPTION: Number of past outputs assumed to be in z state [y_k, y_{k-1}, ..., y_{k-ny+1}]
# Set this based on your system identification / state definition.
# If your state z doesn't explicitly contain past y's like this, adjust the initialization function.

# --- Load or define REAL historical data (NEW) ---
# These should contain the actual measurements/inputs right before time k=0
# Order: Oldest first, newest last.
# y_history_real[end] should be y_0
# u_phys_history_real[end] should be u_{phys,-1}

X_data, U, N_X, N_U, TOTAL_SAMPLES = load_and_preprocess_data(history_path)

const z_history_real = X_data[:, 1:10]
# u_phys_history_real needs shape (m_phys, num_steps). We need the last N_input_delays steps.
const u_phys_history_real = U[:, 1:10]
# --- End of Placeholders ---


# --- 2. Load System Matrices and Infer Dimensions ---
println("Loading system matrices from ./ABC.jld2...")
global A::Matrix{T}, B::Matrix{T} # Ensure type stability
try
    # Use JLD2.load (corrected function name)
    AB_data = JLD2.load("./ABC.jld2")
    if !haskey(AB_data, "A") || !haskey(AB_data, "B")
        error("File ./ABC.jld2 must contain variables named 'A' and 'B'.")
    end
    global A = T.(AB_data["A"]) # Ensure type T
    global B = T.(AB_data["B"]) # Ensure type T
catch e
    error("Failed to load or access data from ./ABC.jld2. Error: $e")
end

# --- Infer Dimensions from loaded matrices ---
const N = size(A, 1) # Inferred state dimension (lifted state z)
const m = size(B, 2) # Inferred model input dimension (potentially large input u)

@info "--- System Dimensions ---"
@info "Inferred State Dimension N = $N (from A)"
@info "Inferred Model Input Dimension m = $m (from B)"
@info "User-defined Physical Input Dimension m_phys = $m_phys"

# --- Sanity Checks ---
if m < m_phys
    error("Inferred model input dimension m ($m) cannot be smaller than physical input dimension m_phys ($m_phys). Check m_phys setting or the loaded B matrix.")
end
if N <= 0 || m <= 0
    error("Inferred dimensions N ($N) and m ($m) must be positive.")
end
if size(A, 2) != N
    error("Loaded matrix A must be square (N x N). Found dimensions $(size(A)).")
end
if size(B, 1) != N
    error("Loaded matrix B must have N rows (N x m). Found dimensions $(size(B)).")
end
# (NEW) Check if assumed ny_init is compatible with state dimension N
if ny_init > N
    @warn "Assumed number of output delays ny_init ($ny_init) is larger than state dimension N ($N). Adjust ny_init or check state definition."
    # Depending on the strictness needed, you might make this an error.
end


# --- Infer Number of Input Delays (nu for initialization) ---
local N_input_delays::Int # Corresponds to 'nu' in initialization context
if m == m_phys
    N_input_delays = 0
    @info "Model input m matches physical input m_phys. Assuming no input delays (nu=0)."
elseif m % m_phys == 0
    N_input_delays = div(m, m_phys) - 1
    @info "Inferred N_input_delays (nu) = $N_input_delays based on m=$m and m_phys=$m_phys."
    if N_input_delays < 0
        error("Calculated N_input_delays is negative. This should not happen if m >= m_phys.") # Should be caught earlier
    end
else
    # This case makes initialization difficult, as the structure of B is unclear.
    @error """Model input dimension m ($m) is not an integer multiple of m_phys ($m_phys).
              Cannot reliably infer input delay structure (nu) for initialization or simulation.
              Check the loaded B matrix and m_phys definition."""
    # Stop execution as initialization and simulation loop logic depend on this structure.
    error("Ambiguous input delay structure (m vs m_phys). Stopping.")
    # N_input_delays = -1 # Flag indicating ambiguity / assuming no delay handling needed (OLD, now error)
end
const nu_init = N_input_delays # Assign to const for clarity in initialization


@info "--- MPC Parameters ---"
@info "Prediction Horizon Np = $Np"
@info "Cost Weights: Q_state=$Q_state_val, R_phys=$R_phys_val, Q_term=$Q_term_val"
@info "Input Constraint |u_phys| <= $u_phys_max"
@info "--- Initialization Parameters ---" # (NEW)
@info "Assumed Output Delays in State ny = $ny_init"
@info "Inferred Input Delays nu = $nu_init"

# --- 3. Prepare Sparse Inputs for the Library's Setup Function ---
# (Code identical to before - preparing Q, R, q, r, E, F, b vectors)
# --- Cost Function Matrices ---
Q_stage = spzeros(T, N, N);
if N > 0
    Q_stage[1, 1] = Q_state_val
else
    @warn "State dimension N=0, cannot set Q_stage cost."
end
Q_term = spzeros(T, N, N);
if N > 0
    Q_term[1, 1] = Q_term_val
else
    @warn "State dimension N=0, cannot set Q_term cost."
end
R_phys_block = sparse(Diagonal(fill(R_phys_val, m_phys)))
if m > m_phys
    R_stage = blockdiag(R_phys_block, spzeros(T, m - m_phys, m - m_phys))
elseif m == m_phys
    R_stage = R_phys_block
else
    error("Logic error: m < m_phys detected again.") # Should be caught earlier
end
q_stage = zeros(T, N);
q_term = zeros(T, N);
r_stage = zeros(T, m);
if N > 0 && Q_stage[1, 1] != 0 # Avoid multiplying by zero if Q_state_val is 0
    q_stage[1] = -target_state_val * Q_stage[1, 1] # Simplified linear term q = -y_target * Q
end
if N > 0 && Q_term[1, 1] != 0
    q_term[1] = -target_state_val * Q_term[1, 1]
end
Q_vec = AbstractMatrix{T}[copy(Q_stage) for _ = 1:Np];
push!(Q_vec, Q_term);
R_vec = AbstractMatrix{T}[copy(R_stage) for _ = 1:Np]
q_vec = AbstractVector{T}[copy(q_stage) for _ = 1:Np];
push!(q_vec, q_term);
r_vec = AbstractVector{T}[copy(r_stage) for _ = 1:Np]

# --- Constraint Matrices: |u_phys_k| <= u_phys_max ---
E_stage = spzeros(T, 2 * m_phys, N)
F_phys_block = sparse([I(m_phys); -I(m_phys)]) # 2*m_phys x m_phys
if m > m_phys
    F_stage = [F_phys_block spzeros(T, 2 * m_phys, m - m_phys)] # 2*m_phys x m
elseif m == m_phys
    F_stage = F_phys_block # 2*m_phys x m
else
    error("Logic error: m < m_phys detected again.") # Should be caught earlier
end
b_stage = ones(T, 2 * m_phys) * u_phys_max # RHS of constraints
E_term = spzeros(T, 0, N);
b_term = zeros(T, 0);
E_vec = AbstractMatrix{T}[copy(E_stage) for _ = 1:Np];
push!(E_vec, E_term);
F_vec = AbstractMatrix{T}[copy(F_stage) for _ = 1:Np]
b_vec = AbstractVector{T}[copy(b_stage) for _ = 1:Np];
push!(b_vec, b_term);


# --- 4. Offline Step: Setup Controller using the Library ---
# This performs the dense matrix computations and pre-allocates the QP solver model
println("Setting up MPC controller...")
controller = setup_dense_mpc(A, B, Np, Q_vec, R_vec, q_vec, r_vec, E_vec, F_vec, b_vec, m_phys, solver_time_limit=1.0);
println("MPC controller setup complete.")


# --- 5. Simulation Setup ---


# --- Initialize Simulation Variables ---
states = Matrix{T}(undef, N, sim_steps + 1)
phys_inputs = Matrix{T}(undef, m_phys, sim_steps)

# --- Initialize State and Input History using the new function ---
println("Initializing state and input history from data...")
try
    global initial_state_z0, initial_input_history # Make available outside try block
    initial_state_z0, initial_input_history = initialize_system_state_and_history(
        z_history_real,
        u_phys_history_real,
        nu_init,
        N,
        m_phys
    )
    states[:, 1] = initial_state_z0
    # Make a mutable copy for the simulation loop to update
    global input_history = copy(initial_input_history)

    log_state_val = (N > 0) ? round(states[1, 1], digits=3) : "N/A (N=0)"
    @info "Initialization complete. Initial State z₀[1] = $log_state_val"
    if nu_init > 0
        log_input_hist_val = (m_phys > 0 && nu_init > 0) ? round(input_history[1, 1], digits=3) : "N/A"
        @info "Initial input history buffer[1,1] (u_phys{-1}): $log_input_hist_val"
    end

catch e
    @error "Failed during state/history initialization: $e"
    @error "Check historical data arrays (y_history_real, u_phys_history_real) and dimensions (ny_init, nu_init)."
    error("Stopping due to initialization failure.") # Stop script
end


# --- Timing and Logging Variables ---
online_times_ns = zeros(Int64, sim_steps)
successful_solves = 0
total_online_time_ns = 0

# --- 6. Online Simulation Loop ---
println("\nStarting MPC simulation loop (sim_steps=$sim_steps)...")
for k = 1:sim_steps
    t_start_ns = time_ns()

    current_state_z = states[:, k]

    # --- Compute Control Action using the Library ---
    u_phys_k, status = compute_control_action(controller, current_state_z)
    # --- (End of Library Call) ---

    t_end_ns = time_ns() # End timing online step
    online_times_ns[k] = t_end_ns - t_start_ns
    global total_online_time_ns += online_times_ns[k]

    if status == JuMP.OPTIMAL || status == JuMP.ALMOST_OPTIMAL
        global successful_solves += 1
        # (Optional) Add handling for non-optimal solves, e.g., reuse previous input
        # else
        #    @warn "Step $k: QP solver did not find optimal solution (Status: $status). Reusing previous input."
        #    if k > 1
        #        u_phys_k = phys_inputs[:, k-1]
        #    else
        #        # Handle first step failure - maybe apply zero input?
        #        u_phys_k = zeros(T, m_phys)
        #        @warn "Step 1: QP failed, applying zero input."
        #    end
    end

    # Clamp physical input (redundant if constraints are well-enforced, but good practice)
    u_phys_k .= clamp.(u_phys_k, -u_phys_max, u_phys_max)

    # Store physical input
    phys_inputs[:, k] = u_phys_k

    # --- Construct full delayed input vector u_k for system simulation ---
    local uk::Vector{T}
    # Use nu_init (N_input_delays) which was validated earlier
    if nu_init > 0
        # Construct full input uk = [u_phys_k; u_phys_{k-1}; ... ; u_{phys,k-nu}]
        uk = zeros(T, m) # m = m_phys * (nu_init + 1)
        uk[1:m_phys] = u_phys_k
        # input_history holds [u_{k-1}, u_{k-2}, ..., u_{k-nu}]
        uk[m_phys+1:end] = vec(input_history) # Flatten history matrix column-major

        # Update history buffer for next step: shift columns right, add new input at the start
        # input_history becomes [u_k, u_{k-1}, ..., u_{k-nu+1}]
        if size(input_history, 2) == nu_init # Double check size
            if nu_init >= 1 # Only shift if there's history
                # Shift existing history
                # Use view to avoid temporary allocation if possible, though copy might be clearer
                input_history[:, 2:end] .= @view input_history[:, 1:end-1]
                # Insert new input
                input_history[:, 1] = u_phys_k
            end
        else
            # This should not happen if initialization was correct
            @error "Mismatch between nu_init ($nu_init) and input_history size ($(size(input_history))) during simulation. Simulation step might be incorrect."
        end

    elseif nu_init == 0 # m == m_phys
        # No delays modeled in B, use only current physical input
        uk = u_phys_k
        # else # nu_init < 0 case already caused an error earlier
    end

    # --- Apply input u_k to the system dynamics ---
    if size(uk, 1) != size(B, 2)
        @error "CRITICAL DIMENSION MISMATCH during simulation step $k:"
        @error "  State vector z_k size: $(size(states[:, k]))"
        @error "  Input vector uk size: $(size(uk)) (Should be $m)"
        @error "  Matrix A size: $(size(A))"
        @error "  Matrix B size: $(size(B))"
        @error "  Required: size(uk, 1) == size(B, 2) ($(size(uk, 1)) != $(size(B, 2)))"
        error("Stopping simulation due to critical dimension mismatch.") # Stop execution
    end
    states[:, k+1] = A * states[:, k] + B * uk

    # --- Logging ---
    if k % 50 == 0 || k == sim_steps # Log less frequently
        log_state_val = (N > 0) ? round(states[1, k+1], digits=3) : "N/A (N=0)"
        log_input_val = (m_phys > 0) ? round(u_phys_k[1], digits=3) : "N/A (m_phys=0)"
        @info "Step $k/$sim_steps: State[1] = $log_state_val, Input[1] = $log_input_val, Solve Status: $status, Time: $(round(online_times_ns[k]/1e6, digits=2)) ms"
    end
end

@info "Simulation complete. Successful QP solves: $successful_solves / $sim_steps"

# --- 7. Report Timing Statistics ---
# (Identical to before)
if successful_solves > 0 && sim_steps > 0
    valid_times_ns = online_times_ns[1:sim_steps]
    avg_online_time_ms = mean(valid_times_ns) / 1e6
    max_online_time_ms = maximum(valid_times_ns) / 1e6
    min_online_time_ms = minimum(valid_times_ns) / 1e6
    @info "Online MPC Step Timing (avg / max / min): $(round(avg_online_time_ms, digits=4)) ms / $(round(max_online_time_ms, digits=4)) ms / $(round(min_online_time_ms, digits=4)) ms"
else
    @warn "No timing statistics available (no successful solves or zero simulation steps)."
end


# --- 8. Plotting ---
# (Modified slightly to handle potential plotting issues if GLMakie not fully working)
println("Preparing plots...")
try
    fig = Figure(size=(800, 700))

    # State Plot
    ax_state = Axis(fig[1, 1], title="State Trajectory (Element 1, N=$N)", xlabel="Time step k", ylabel="z₁(k)")
    if N > 0
        lines!(ax_state, 0:sim_steps, states[1, :], label="State z₁")
        hlines!(ax_state, [target_state_val], color=:red, linestyle=:dash, label="Target")
    else
        text!(ax_state, "State dimension N=0,\ncannot plot state trajectory.", position=(sim_steps / 2, 0), align=(:center, :center))
    end
    axislegend(ax_state, position=:rt)

    # Input Plot
    ax_input = Axis(fig[2, 1], title="Input Trajectory (Element 1, m_phys=$m_phys)", xlabel="Time step k", ylabel="u_phys,₁(k)")
    if m_phys > 0 && sim_steps > 0 # Check sim_steps > 0 for stairs
        # Plot only the first physical input for clarity if m_phys is large
        plot_m_idx = 1
        stairs!(ax_input, 1:sim_steps, phys_inputs[plot_m_idx, :], step=:post, label="Input u_phys, $plot_m_idx")
        hlines!(ax_input, [u_phys_max, -u_phys_max], color=:red, linestyle=:dash, label="|u| <= $u_phys_max")
        # Auto limits usually work well, but can set manually if needed
        # ylims!(ax_input, -u_phys_max * 1.1, u_phys_max * 1.1)
    else
        text!(ax_input, "Physical input dimension m_phys=0\nor sim_steps=0, cannot plot input.", position=(sim_steps / 2, 0), align=(:center, :center))
    end
    axislegend(ax_input, position=:rt)


    # Display or Save Plot
    display(fig) # Try to display interactive plot

    filename = "mpc_inferred_N$(N)_m$(m)_mp$(m_phys)_Np$(Np)_init.png"
    try
        save(filename, fig)
        println("Figure saved to $filename")
    catch e_save
        println("Failed to save figure: $e_save")
    end

catch e_plot
    @error "Failed during plotting: $e_plot"
    @warn "Ensure GLMakie and a suitable backend are working correctly."
end


println("\n--- Final Values ---")
final_state_val = (N > 0) ? states[1, end] : "N/A (N=0)"
println("Final state (element 1): ", final_state_val)
# Avoid norm calculation if N=0
final_norm_val = (N > 0) ? norm(states[:, end]) : "N/A (N=0)"
println("Final state norm (lifted): ", final_norm_val)

println("\nScript finished.")
