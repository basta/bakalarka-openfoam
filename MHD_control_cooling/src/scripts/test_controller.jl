# =============== SIMULATION SCRIPT (Updated) ===============
using LinearAlgebra
using Plots
using CircularArrayBuffers # Needed by DMDC_MPCController
using JLD2               # Needed by DMDC_MPCController & dmd_lib
using SparseArrays       # Needed by DMDC_MPCController
using JuMP               # Needed by DMDC_MPCController & mpc_lib
using Statistics         # Needed by dmd_lib
using BlockArrays        # Needed by mpc_lib
using OSQP               # Needed by mpc_lib (as the default solver)
using JLD2
# using Pkg; Pkg.add(["LinearAlgebra", "Plots", "CircularArrayBuffers", "JLD2", "SparseArrays", "JuMP", "Statistics", "BlockArrays", "OSQP"])

println("Loading controller and library code...")

# --- Load Actual Libraries ---
# NOTE: dmd_lib.jl defines global constants and a main() function which might
#       cause conflicts or unexpected behavior when included as a library.
#       Consider refactoring dmd_lib.jl into a module containing only the
#       required functions (e.g., create_abc_system_from_data).
try
    include("../dmd_lib.jl")
    println("Loaded dmd_lib.jl")
catch e
    @error "Failed to load dmd_lib.jl" exception = (e, catch_backtrace())
    rethrow(e)
end

try
    # Note: mpc_lib.jl is defined as a module (DenseMPCLib)
    include("../mpc/mpc_lib.jl")
    println("Loaded mpc_lib.jl")
    # Make the module available
    using .DenseMPCLib
catch e
    @error "Failed to load mpc_lib.jl" exception = (e, catch_backtrace())
    rethrow(e)
end


# --- Controller Interface ---
try
    include("../controllers/interface.jl")
    println("Loaded interface.jl")
catch e
    @error "Failed to load interface.jl" exception = (e, catch_backtrace())
    rethrow(e)
end

# --- Identification Controller ---
try
    include("../controllers/identification.jl")
    println("Loaded identification.jl")
catch e
    @error "Failed to load identification.jl" exception = (e, catch_backtrace())
    rethrow(e)
end
# --- Use the actual MPC controller file ---
println("Loading DMDC_MPCController code...")
# Note: mpc_dmdc.jl includes dmd_lib.jl and mpc_lib.jl itself.
# Including them globally here might be redundant but usually okay.
# It also uses the OuterProductOptimPlaceholder defined above.
try
    include("../controllers/mpc_dmdc.jl")
    println("Loaded mpc_dmdc.jl")
catch e
    @error "Failed to load mpc_dmdc.jl" exception = (e, catch_backtrace())
    rethrow(e)
end

# ==============================================
# --- System Definition ---
println("Defining system...")
# Example: A simple discrete-time linear system x_k+1 = A*x_k + B*u_k
# Let state dim n=5, input dim m=8 (matching default m_phys in DMDC_MPC)
# Adjust N_X_System based on the number of states your DMDc model uses (from RELEVANT_STATES in dmd_lib.jl)
# If RELEVANT_STATES = [6], then N_X_System = 1
# If RELEVANT_STATES = [1, 2, 3], then N_X_System = 3
# This MUST match the N_X inferred/used by DMDC_MPCController
# N_X_System = 1 # Set based on RELEVANT_STATES in dmd_lib.jl # Commented out - n determined by loaded A
# n = N_X_System   # State dimension for the simulation # Commented out - n determined by loaded A
# m = 16           # Input dimension (matches m_phys) # Commented out - m determined by loaded B

# Example stable A matrix (replace with your actual system)
# Dimensions must match n x n

local A::Matrix{Float64}, B::Matrix{Float64}
try
    jldopen("ABC.jld2", "r") do file
        global A = file["A"]
        global B = file["B"]
        # Optionally load delay and mean_vec if needed elsewhere, though simulation loop doesn't use them directly
        # global delay = file["delay"]
        global mean_vec = file["mean_vec"]
    end
catch e
    @error "Failed to load A and B matrices from ABC.jld2. Ensure the file exists and contains 'A' and 'B'." exception = (e, catch_backtrace())
    # Provide dummy matrices to allow script to potentially continue, or rethrow
    @warn "Using dummy A and B matrices."
    global A = Matrix(Diagonal([0.9])) # Dummy 1x1
    global B = zeros(1, 16)            # Dummy 1x16
    # rethrow(e) # Uncomment to stop execution if file load fails
end


# Determine dimensions from loaded matrices
n = size(A, 1) # Actual state dimension from loaded A
m = size(B, 2) # Actual input dimension from loaded B

println("System Matrices Loaded/Defined:")
println("A ($n x $n)")
println("B ($n x $m)")

# ==============================================
# --- Simulation Parameters ---
println("Setting up simulation parameters...")
T_sim = 1000.0 # Total simulation time
dt = 5.0     # Time step (assuming controllers work in discrete time)
num_steps = Int(round(T_sim / dt))
# Initial state dimension must match n (from loaded A)
x0 = zeros(n) # Initial state
# x0[1:floor(Int,n/2)+1] .= mean_vec
@info mean_vec
# x0[1] = mean_vec[1]
println("T_sim = $T_sim, dt = $dt, Steps = $num_steps")

# ==============================================
# --- Controller Instantiation ---
println("Instantiating controllers...")

# 1. Identification Controller
# Ensure num_inputs matches the actual input dimension 'm' from loaded B
ident_controller = IdentificationController(100.0, m) # T_stabilize=100, num_inputs=m

# 2. DMDC MPC Controller
#    Uses the actual dmd_lib.jl and mpc_lib.jl now.
#    Still requires outer_product_optim.jl (using placeholder).
#    !! Adjust this config based on your actual setup !!
mpc_config = Dict(
    "controller" => Dict(
        "type" => "DMDC_MPC",
        "params" => Dict(
            "DMDC_MPC" => Dict(
                "Np" => 10,             # Prediction horizon
                "m_phys" => 8,          # Physical input dimension MUST MATCH system 'm' (from loaded B)
                # Use the DATA_PATH defined in dmd_lib.jl or override here
                "data_path" => "./data.jld2", # Path for DMDc ID (ensure this file exists and is readable)
                "abc_path" => "ABC.jld2", # Path to load/save A/B (will be created by dmd_lib if load_abc=false)
                "load_abc" => true,    # Set to true since we loaded A/B above for the simulation
                # Use the DELAY defined in dmd_lib.jl or override here (MUST match delay used to generate ABC.jld2 if load_abc=true)
                "delay" => 0,           # System delay steps
                # Add placeholder paths for costs/constraints if needed by your actual MPC setup
                # Example cost values (used by the placeholder mpc_dmdc.jl)
                "Q_state_val" => 1.0,
                "R_phys_val" => 0.1,
                "Q_term_val" => 10.0,
                # Example constraints (used by placeholder mpc_dmdc.jl)
                "u_phys_max" => 1.5,
                # N_X is inferred during identification in dmd_lib.jl based on RELEVANT_STATES
                # or loaded from ABC.jld2 if load_abc=true
            )
        )
    )
)

# Ensure data file exists if not loading ABC (though load_abc is now true)
data_file_path = mpc_config["controller"]["params"]["DMDC_MPC"]["data_path"]
abc_file_path = mpc_config["controller"]["params"]["DMDC_MPC"]["abc_path"]
if !mpc_config["controller"]["params"]["DMDC_MPC"]["load_abc"] && !isfile(data_file_path)
    @error "Data file '$data_file_path' not found, which is required for DMDc identification (load_abc=false)."
    # Consider stopping execution or providing instructions
end
# Clean up any old ABC file if recomputing (now disabled as load_abc=true)
# if !mpc_config["controller"]["params"]["DMDC_MPC"]["load_abc"] && isfile(abc_file_path)
#     @info "Removing old ABC file specified in config: $abc_file_path (will be regenerated by dmd_lib)"
#     try
#         rm(abc_file_path)
#     catch e
#         @warn "Could not remove old ABC file: $e"
#     end
# end


# Instantiate the MPC controller
# Wrap in try-catch as it involves file I/O and potential dimension mismatches
local mpc_controller # Use local for safety within try-catch
try
    # Pass the config dictionary to the constructor
    global mpc_controller = DMDC_MPCController(mpc_config)
catch e
    @error "Failed to instantiate DMDC_MPCController. Check configuration, library paths, data files, and RELEVANT_STATES in dmd_lib.jl." exception = (e, catch_backtrace())
    global mpc_controller = nothing # Indicate failure
end

# ==============================================
# --- Simulation Function ---
function simulate_system(controller::AbstractController, A::Matrix{Float64}, B::Matrix{Float64}, x0::Vector{Float64}, dt::Float64, num_steps::Int)
    n_sim = size(A, 1) # Augmented state dimension
    m_sim = size(B, 2) # Model input dimension (N_U = 16)
    x_hist = zeros(n_sim, num_steps + 1)
    # Determine physical input dimension expected by the controller (should be 8 for MPC)
    m_phys_ctrl = hasfield(typeof(controller), :m_phys) ? controller.m_phys : m_sim
    u_ctrl_hist = zeros(m_phys_ctrl, num_steps) # Store the 8-dim control from MPC
    u_model_hist = zeros(m_sim, num_steps) # Store the 16-dim input used in model

    t_hist = zeros(num_steps + 1)

    x_hist[:, 1] = x0
    t_hist[1] = 0.0
    x_current = copy(x0) # Augmented state z_k

    N_X_ctrl = hasfield(typeof(controller), :N_X) ? controller.N_X : n_sim # Physical state dim

    println("\n--- Starting Simulation: $(typeof(controller)) ---")
    # ... (keep print statements) ...

    sim_start_time = time()

    for k = 1:num_steps
        t_current = k * dt
        t_hist[k+1] = t_current

        # --- State Formatting for Controller ---
        state_dict = Dict{Symbol,Any}()
        state_dict[:time] = t_current
        state_dict[:state_vector] = x_current # Pass augmented state z_k
        # Map relevant part of physical state (first N_X elements of z_k, which are normalized)
        # estimate_physical_state needs raw data, but here we only have z_k.
        # We need to pass info so the controller can reconstruct its required history.
        # The controller *should* use its internal history buffers.
        # Let's assume estimate_physical_state works correctly if given the right keys from a real system.
        # For this simulation, the controller MUST rely on its internal history updated in the previous step.
        # We can pass the current *normalized* physical state part under a suitable key if needed.
        state_dict[:normalized_physical_state] = x_current[1:N_X_ctrl] # Pass z_k[1:N_X]

        # --- Compute Control Action ---
        # Controller should return the 8-dimensional physical vector [u_e; u_m]
        local u_k_ctrl::Vector{Float64} # This is [u_e; u_m]
        try
            u_k_ctrl = compute_control_action(controller, t_current, state_dict)
        catch err
            @error "Error during compute_control_action for $(typeof(controller)) at step $k" exception = (err, catch_backtrace())
            @warn "Applying zero control input as fallback."
            # Determine fallback based on controller type if needed
            u_k_ctrl = zeros(m_phys_ctrl) # Fallback to zeros (size 8 for MPC)
        end
        u_ctrl_hist[:, k] = u_k_ctrl

        # --- Reconstruct Model Input u'_k from Physical Control u_k_ctrl --- ## <--- NEW SECTION
        local u_k_sim::Vector{Float64} # This is the 16-dim u'_k = vec(u_e * u_m')
        if m_phys_ctrl == 8 && m_sim == 16 # Specific check for the outer product case
             if length(u_k_ctrl) == 8
                 u_e_sim = u_k_ctrl[1:4]
                 u_m_sim = u_k_ctrl[5:8]
                 u_k_sim = vec(u_e_sim * u_m_sim') # Reconstruct the 16-dim u'
             else
                 @warn "[Simulation Step $k] Controller returned unexpected dimension $(length(u_k_ctrl)) for physical input. Expected 8. Using zeros for model input."
                 u_k_sim = zeros(m_sim) # m_sim is 16
             end
        else
             # Fallback for non-MPC controllers or different dimensions
             @warn "[Simulation Step $k] Handling non-standard input dimensions (m_phys=$(m_phys_ctrl), m_sim=$(m_sim)). Padding/truncating control output."
             u_k_sim = zeros(m_sim)
             len_copy = min(length(u_k_ctrl), m_sim)
             u_k_sim[1:len_copy] = u_k_ctrl[1:len_copy]
        end
        u_model_hist[:,k] = u_k_sim # Store the 16-dim input used by model


        # --- Update System State ---
        # Use the 16-dimensional u_k_sim (reconstructed u')
        try
            x_next = A * x_current + B * u_k_sim # z_{k+1} = A_full*z_k + B_full*u'_k
            x_hist[:, k+1] = x_next
            x_current = x_next # Update augmented state for next loop
        catch err
            @error "Error during state update z_k+1 = A*z_k + B*u'_k at step $k" exception = (err, catch_backtrace())
            @warn "State update failed. Simulation might become unstable."
            break # Stop the simulation loop on error
        end

        # --- Update Controller Internal State (if applicable) ---
        # This typically requires the *new* physical state measurement.
        # In simulation, we provide the newly calculated state `x_current` (z_{k+1}).
        new_state_dict = Dict{Symbol,Any}()
        new_state_dict[:time] = t_current # Time corresponding to measurement leading to x_next
        new_state_dict[:state_vector] = x_current # Pass z_{k+1}
        new_state_dict[:normalized_physical_state] = x_current[1:N_X_ctrl] # Pass z_{k+1}[1:N_X]

        try
            update_controller_state!(controller, t_current, new_state_dict)
        catch err
            @error "Error during update_controller_state! for $(typeof(controller)) at step $k" exception = (err, catch_backtrace())
        end

        # --- Progress ---
        # ... (keep progress print) ...
    end

    sim_end_time = time()
    println("--- Simulation Finished: $(typeof(controller)) ---")
    println("Total time: $(round(sim_end_time - sim_start_time, digits=2)) seconds")

    # Return history of augmented state, physical control, and model input used
    return t_hist, x_hist, u_ctrl_hist, u_model_hist
end

# ==============================================
# --- Run Simulations ---

# --- Simulation 1: Identification Controller ---
global t_id, x_id, u_ctrl_id, u_model_id
try
    # Run simulation for Identification Controller
    # This controller directly outputs m=16 values in this setup
    global t_id, x_id, u_ctrl_id, u_model_id = simulate_system(ident_controller, A, B, x0, dt, num_steps)
catch e
    @error "Simulation failed for IdentificationController" exception=(e, catch_backtrace())
    global t_id, x_id, u_ctrl_id, u_model_id = ntuple(_ -> nothing, 4)
end

global t_mpc, x_mpc, u_ctrl_mpc, u_model_mpc
if !isnothing(mpc_controller)
    try
        # Run simulation for MPC Controller
        # Should return 8-dim physical control in u_ctrl_mpc
        global t_mpc, x_mpc, u_ctrl_mpc, u_model_mpc = simulate_system(mpc_controller, A, B, x0, dt, num_steps)
    catch e
        @error "Simulation failed for DMDC_MPCController" exception=(e, catch_backtrace())
        global t_mpc, x_mpc, u_ctrl_mpc, u_model_mpc = ntuple(_ -> nothing, 4)
    end
else
     @warn "Skipping DMDC_MPC simulation because instantiation failed."
     global t_mpc, x_mpc, u_ctrl_mpc, u_model_mpc = ntuple(_ -> nothing, 4)
end



# --- Simulation 2: DMDC MPC Controller ---
global t_mpc, x_mpc, u_mpc # Ensure scope
if !isnothing(mpc_controller)
    try
        # Reset initial state if desired, or continue from x0
        global t_mpc, x_mpc, u_mpc = simulate_system(mpc_controller, A, B, x0, dt, num_steps)
    catch e
        @error "Simulation failed for DMDC_MPCController" exception = (e, catch_backtrace())
        global t_mpc, x_mpc, u_mpc = nothing, nothing, nothing
    end
else
    @warn "Skipping DMDC_MPC simulation because instantiation failed."
    global t_mpc, x_mpc, u_mpc = nothing, nothing, nothing # Mark as not run
end


# ==============================================
# --- Plotting Results ---
println("Plotting results...")

# Check if simulation results exist before plotting
if !isnothing(x_id) && (!isnothing(u_ctrl_id) || !isnothing(u_ctrl_mpc))

    # --- Plot States ---
    # Option 1: Plot NORMALIZED physical states (first N_X components of x_id/x_mpc)
    local N_X::Int
    # Determine N_X (try from controller, then from ABC file)
    if !isnothing(mpc_controller); N_X = mpc_controller.N_X
    else try; jldopen("ABC.jld2", "r") do f; N_X = length(f["mean_vec"]) end catch; N_X=0 end end

    if N_X > 0
        num_states_to_plot = min(N_X, 3)
        p_state_norm = plot(xlabel="Time (steps*dt)", ylabel="Normalized State Value", title="Normalized Physical State Trajectories (z_k[1:$N_X])", legend=:outertopright)
        for i in 1:num_states_to_plot
             plot!(p_state_norm, t_id, x_id[i, :], label="Ident Ctrl - Norm x$i", lw=1.5)
             if !isnothing(t_mpc) && !isnothing(x_mpc) && size(x_mpc, 1) >= i
                 plot!(p_state_norm, t_mpc, x_mpc[i, :], label="MPC Ctrl - Norm x$i", lw=1.5, linestyle=:dash)
             end
        end
    else
         p_state_norm = plot(title="State Plot (N_X unknown)") # Placeholder if N_X unknown
    end

    # --- Plot Controls ---
    # Plot the 8-dimensional PHYSICAL control computed by MPC
    # And the 16-dim control from Identification controller (first few elements)
    p_input_phys = plot(xlabel="Time (steps*dt)", ylabel="Input Value", title="Physical Control Signals ([u_e;u_m] for MPC)", legend=:outertopright)
    time_input = t_id[1:num_steps]

    # Plot first 3 elements of Ident controller output (likely size 16)
    num_ident_inputs_to_plot = min(size(u_ctrl_id, 1), 3)
    for i in 1:num_ident_inputs_to_plot
         plot!(p_input_phys, time_input, u_ctrl_id[i, :], label="Ident Ctrl - u$i", lw=1.5, drawstyle=:steps)
    end

    # Plot first 3 elements of MPC controller output (physical [u_e;u_m], size 8)
    if !isnothing(t_mpc) && !isnothing(u_ctrl_mpc)
        num_mpc_inputs_to_plot = min(size(u_ctrl_mpc, 1), 3) # Max 8, plot up to 3
        for i in 1:num_mpc_inputs_to_plot
            plot!(p_input_phys, time_input, u_ctrl_mpc[i, :], label="MPC Ctrl - Phys u$i", lw=1.5, linestyle=:dash, drawstyle=:steps)
        end
    end

    # Combine plots
    p_combined = plot(p_state_norm, p_input_phys, layout=(2, 1), size=(800, 700))
    display(p_combined)
    try savefig(p_combined, "simulation_results_norm_physCtrl.png"); println("Plot saved.") catch e @warn "Plot save failed." end

else
    @warn "One or both simulations failed or N_X unknown. Skipping plotting."
end

# Keep plots open in interactive mode (Optional)
# println("Press Enter to exit...")
# readline()

println("\nScript finished.")
