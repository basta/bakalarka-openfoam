# Import necessary packages
using LinearAlgebra # For matrix operations
using OSQP        # For solving the Quadratic Program (QP)
using Plots       # For plotting results
using SparseArrays # For sparse matrix operations
using JLD2        # For loading data

include("../optim/outer_product_optim.jl")

# Define a structure to hold MPC data
struct MPCData
    A::Matrix{Float64}
    # --- FIX 1: Changed B to Matrix{Float64} ---
    B::Matrix{Float64} # Store B as a Matrix
    nx::Int
    nu::Int
    Np::Int
    Nc::Int
    Q::Diagonal{Float64, Vector{Float64}}
    R::Diagonal{Float64, Vector{Float64}}
    xref::Vector{Float64}
    umin::Vector{Float64}
    umax::Vector{Float64}
    xmin::Vector{Float64}
    xmax::Vector{Float64}
    Sx::Matrix{Float64}
    Su::Matrix{Float64}
    Qbar::SparseMatrixCSC{Float64, Int64}
    Rbar::SparseMatrixCSC{Float64, Int64}
    P_sparse::SparseMatrixCSC{Float64, Int64}
    A_in::SparseMatrixCSC{Float64, Int64} # Note: Was SparseMatrixCSC{Bool, Int64} in error trace, ensure Float64
    A_st::SparseMatrixCSC{Float64, Int64} # Store Su as sparse A_st
    A_con::SparseMatrixCSC{Float64, Int64}
    l_in::Vector{Float64}
    u_in::Vector{Float64}
    Xmin::Vector{Float64}
    Xmax::Vector{Float64}
    Xref::Vector{Float64}
    model::OSQP.Model # Store the OSQP model instance
end

"""
    setup_mpc(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax; osqp_settings...)

Performs offline calculations for the dense MPC controller and initializes the OSQP solver.

Args:
    A, B: State-space matrices (B should be nx x nu matrix).
    Np, Nc: Prediction and control horizons.
    Q, R: Weighting matrices.
    xref: Reference state vector.
    umin, umax: Input constraints (vectors).
    xmin, xmax: State constraints (vectors).
    osqp_settings: Optional keyword arguments for OSQP.setup! (e.g., verbose=false).

Returns:
    An MPCData struct containing pre-calculated matrices and the initialized OSQP model.
"""
function setup_mpc(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax; osqp_settings...)
    nx = size(A, 1)
    # --- FIX 2: Use B directly, assuming it's already a Matrix ---
    # B_mat = reshape(B, nx, size(B, 2)) # No longer needed if B is passed as Matrix
    nu = size(B, 2) # Get nu directly from input Matrix B

    # --- Dimension Checks ---
    @assert size(A, 1) == size(A, 2) "Matrix A must be square ($nx x $nx)."
    @assert size(B, 1) == nx "Matrix B must have $nx rows (dimension mismatch with A)."
    @assert Np > 0 "Prediction horizon Np must be positive."
    @assert Nc > 0 "Control horizon Nc must be positive."
    # @assert Nc <= Np "Control horizon Nc should typically be <= Np." # Optional check
    @assert size(Q) == (nx, nx) "Matrix Q must be square with dimension $nx x $nx."
    @assert size(R) == (nu, nu) "Matrix R must be square with dimension $nu x $nu."
    @assert length(xref) == nx "Reference vector xref must have length $nx."
    @assert length(umin) == nu "Control lower bound vector umin must have length $nu."
    @assert length(umax) == nu "Control upper bound vector umax must have length $nu."
    @assert length(xmin) == nx "State lower bound vector xmin must have length $nx."
    @assert length(xmax) == nx "State upper bound vector xmax must have length $nx."


    # --- Dense MPC Formulation Calculations ---
    Sx = zeros(nx * Np, nx)
    Su = zeros(nx * Np, nu * Nc)

    # Calculate Sx
    temp_A = A
    for i = 1:Np
        rows = (i-1)*nx+1 : i*nx
        Sx[rows, :] = temp_A
        temp_A = temp_A * A
    end

    # Calculate Su using the input Matrix B directly
    for i = 1:Np
        rows_x = (i-1)*nx+1 : i*nx
        for j = 1:min(i, Nc)
            cols_u = (j-1)*nu+1 : j*nu
            if i - j >= 0
                Su[rows_x, cols_u] = (A^(i - j)) * B # Use B directly
            end
        end
    end

    # Build QP matrices
    Qbar = kron(sparse(I(Np)), sparse(Q))
    Rbar = kron(sparse(I(Nc)), sparse(R))

    Hessian_calc = Su' * Qbar * Su + Rbar
    P_sparse = sparse(Symmetric(Hessian_calc))

    # --- Constraint Formulation ---
    Umin = repeat(umin, Nc)
    Umax = repeat(umax, Nc)
    # Ensure A_in is Float64 if needed, though Identity is usually okay as Bool for structure
    A_in = sparse(Matrix{Float64}(I, nu * Nc, nu * Nc)) # Explicitly Float64 Identity
    l_in = Umin
    u_in = Umax

    Xmin = repeat(xmin, Np)
    Xmax = repeat(xmax, Np)
    A_st = sparse(Su) # Store Su as sparse matrix A_st

    A_con = sparse([A_in; A_st])

    # Stacked reference
    Xref = repeat(xref, Np)

    # --- OSQP Solver Setup ---
    model = OSQP.Model()
    # Default OSQP settings
    default_settings = Dict(:verbose => false, :eps_abs => 1e-4, :eps_rel => 1e-4, :max_iter => 5000)
    # Merge default with user-provided settings
    settings = merge(default_settings, Dict(osqp_settings))

    OSQP.setup!(model; P=P_sparse, q=zeros(nu*Nc), A=A_con, l=zeros(size(A_con,1)), u=zeros(size(A_con,1)), settings...)

    # Create and return the MPCData struct (passing the original Matrix B)
    return MPCData(A, B, nx, nu, Np, Nc, Q, R, xref, umin, umax, xmin, xmax,
                   Sx, Su, Qbar, Rbar, P_sparse, A_in, A_st, A_con,
                   l_in, u_in, Xmin, Xmax, Xref, model)
end

"""
    compute_control_input(xk, mpc_data::MPCData)

Computes the optimal control input for the current state using the pre-configured MPC setup.

Args:
    xk: Current state vector.
    mpc_data: The MPCData struct returned by setup_mpc.

Returns:
    The first optimal control input vector uk (or nothing if solver fails).
"""
function compute_control_input(xk, mpc_data::MPCData)
    # Unpack necessary data from the struct
    (; Sx, Su, Qbar, Xref, Xmin, Xmax, l_in, u_in, A_con, model, nu, nx) = mpc_data

    # Update linear term q based on current state xk and reference Xref
    q_k = Su' * Qbar * (Sx * xk - Xref)

    # Update state constraint bounds based on current state xk
    l_st_k = Xmin - Sx * xk
    u_st_k = Xmax - Sx * xk

    # Combine constraint bounds
    l_k = [l_in; l_st_k]
    u_k = [u_in; u_st_k]

    # Update the QP problem in the OSQP model
    OSQP.update!(model; q=q_k, l=l_k, u=u_k)

    # Solve the QP
    results = OSQP.solve!(model)

    # Check solver status
    if results.info.status_val ∉ (1, 2) # 1: solved, 2: solved inaccurate
        println("Warning: OSQP solver failed with status: $(results.info.status)")
        if isnothing(results.x)
             println("Error: No solution found.")
             return nothing # Indicate failure
        else
            println("Warning: Using potentially suboptimal solution.")
            # Fallthrough to use the suboptimal solution
        end
    end

    # Extract the optimal control sequence U
    U_opt = results.x

    # Extract the first control input
    uk_vec = U_opt[1:nu] # This is a Vector

    return uk_vec # Return the control input vector
end

# --- Main Simulation Script ---

# --- 1. System Definition ---
println("Loading system data from ABC.jld2...")
data = jldopen("./ABC.jld2")
dt = 0.1 # Assuming dt is fixed or not needed if A, B are already discrete
A = data["A"]
B = data["B"] # Load B (expected to be Matrix)
close(data) # Close the JLD2 file
println("System data loaded.")
println("Size of A: ", size(A))
println("Size of B: ", size(B))

# --- 2. MPC Configuration ---
println("Configuring MPC parameters...")
Np = 50
Nc = 10
nx_loaded = size(A, 1)
nu_loaded = size(B, 2)
println("nx = $nx_loaded, nu = $nu_loaded")

# Ensure Q is nx x nx Diagonal
Q_diag = zeros(nx_loaded)
Q_diag[1] = 1.0 # Penalize only the first state deviation
Q = Diagonal(Q_diag)

# Ensure R is nu x nu Diagonal
R = Diagonal(0.1 * ones(nu_loaded)) # Example: Penalize all inputs equally

# Ensure xref has length nx
xref = zeros(nx_loaded)
xref[1] = +1e7 # Set reference for the first state

# Ensure constraints match nu
umin = -5.0 * ones(nu_loaded)
umax = 5.0 * ones(nu_loaded)

# Ensure state constraints match nx
xmin = -Inf * ones(nx_loaded)
xmax = Inf * ones(nx_loaded)
println("MPC parameters configured.")


# --- 3. Perform Offline MPC Setup ---
println("Setting up MPC controller...")
mpc_data = setup_mpc(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax; verbose=false)
println("MPC setup complete.")

# --- 4. Simulation Loop ---
Tsim = 500
x0 = zeros(mpc_data.nx) # Initial state based on nx from mpc_data
xk = copy(x0)

# Store results
X_hist = zeros(mpc_data.nx, Tsim + 1)
U_hist = zeros(mpc_data.nu, Tsim)
X_hist[:, 1] = xk

println("Starting simulation...")
for k = 1:Tsim
    # Compute control input for the current state
    uk_vec = compute_control_input(xk, mpc_data)

    # TODO outer product test
    u, v = find_closest_outer_product(uk_vec, 4, 4)
    @info "Closest outer product found: $(u) * $(v)'"
    uk_vec = vec(u * v')

    if isnothing(uk_vec)
        println("Stopping simulation due to solver failure at step $k.")
        # Decide on fallback action, e.g., apply zero input
        uk_vec = zeros(mpc_data.nu) # Apply zero input as fallback
        # Or break the loop:
        # break
    end

    # --- FIX 3: Use Matrix * Vector multiplication always ---
    # Apply the control input to the system using matrix B
    global xk = mpc_data.A * xk + mpc_data.B * uk_vec # B is Matrix, uk_vec is Vector

    # Store results
    X_hist[:, k+1] = xk
    U_hist[:, k] = uk_vec # Store the vector uk_vec

    if k % 10 == 0
        println("Sim step $k/$Tsim")
    end
end
println("Simulation finished.")

# --- 5. Plotting Results ---
println("Plotting results...")
time = 0:dt:(Tsim*dt) # Use the dt assumed or defined earlier

# Plot only the first few states/inputs if the system is large
plot_states = min(mpc_data.nx, 1) # Plot up to 3 states
plot_inputs = min(mpc_data.nu, 2) # Plot up to 2 inputs

# State plots
plot_list_states = []
for i = 1:plot_states
    p_state = plot(time, X_hist[i, :], label="x$i", xlabel="Time (s)", ylabel="State Value", title="State Trajectories")
    # Plot reference only if it's for the current state index i
    if i <= length(xref) && isfinite(xref[i]) # Check if reference exists and is finite
         plot!(p_state, time, fill(xref[i], length(time)), label="x$i Ref", linestyle=:dash)
    end
    push!(plot_list_states, p_state)
end

# Input plots
plot_list_inputs = []
for i = 1:plot_inputs
     p_input = plot(time[1:Tsim], U_hist[i, :], label="u$i", xlabel="Time (s)", ylabel="Input Value", title="Control Inputs")
     # Plot constraints only if they exist for the current input index i
     if i <= length(umin) && isfinite(umin[i])
         plot!(p_input, time[1:Tsim], fill(umin[i], Tsim), label="u$i min", linestyle=:dash, color=:red)
     end
     if i <= length(umax) && isfinite(umax[i])
         plot!(p_input, time[1:Tsim], fill(umax[i], Tsim), label="u$i max", linestyle=:dash, color=:red)
     end
     push!(plot_list_inputs, p_input)
end

# Combine plots
# Determine layout based on number of plots
total_plots = length(plot_list_states) + length(plot_list_inputs)
ncols = max(1, min(total_plots, 2)) # Arrange in 1 or 2 columns
# Combine all plots into one layout
final_plot = plot(plot_list_states..., plot_list_inputs..., layout=(ncols, ceil(Int, total_plots / ncols)), legend=:best)


# Display the plot
# display(final_plot) # Uncomment to display interactively
savefig("mpc_results_refactored.png")
println("Results saved to mpc_results_refactored.png")
