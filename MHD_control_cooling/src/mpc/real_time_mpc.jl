# real_time_mpc.jl
# This script defines a function to compute MPC control inputs in real-time,
# assuming the computationally expensive setup has been done offline.

using LinearAlgebra
using JuMP
# Assuming OSQP is the solver used in mpc_lib.jl, otherwise adjust
# using OSQP

# Include the library containing the controller struct and setup/compute functions
# Adjust the path if mpc_lib.jl is located elsewhere relative to this script
try
    include("mpc_lib.jl") # Or "../mpc_lib.jl" etc.
    using .DenseMPCLib
catch e
    error("Failed to include mpc_lib.jl. Ensure the file exists at the specified path. Error: $e")
end

"""
    get_control_input(controller::DenseMPCController,
                      current_state_unnormalized::AbstractVector{T},
                      mean_vec::AbstractVector{T},
                      N_X::Int,
                      N_U::Int,
                      DELAY::Int) where {T<:Real}

Computes the optimal physical control input for the current time step using a pre-configured DenseMPCController.

This function takes the *unnormalized* augmented state vector as input, normalizes the
state components (x_k, x_{k-1}, ...), and then calls the `compute_control_action`
function from DenseMPCLib.

# Arguments
- `controller::DenseMPCController`: The pre-configured controller object obtained from `setup_dense_mpc`.
- `current_state_unnormalized::AbstractVector{T}`: The current augmented state vector
  received from the simulation *before* normalization. It must follow the structure:
  `[x_k; x_{k-1}; ...; x_{k-DELAY}; u_{k-1}; ...; u_{k-DELAY}]`, where `x` has dimension `N_X`
  and `u` has dimension `N_U`. The total length must match `controller.N`.
- `mean_vec::AbstractVector{T}`: The mean vector used for normalizing the original state `x`
  (dimension `N_X`). This is the `mean_vec` saved in `ABC.jld2`.
- `N_X::Int`: The dimension of the original (physical) state vector `x`.
- `N_U::Int`: The dimension of the original (physical) input vector `u` (before outer product/embedding).
             Note: This `N_U` corresponds to the dimension of the *preprocessed* `U` in `create_ABC.jl`
             (e.g., 16 after outer product if original u_e/u_m were 4-dim). It must match `controller.m / (DELAY + 1)` if there are input delays, or `controller.m` if `DELAY=0`.
- `DELAY::Int`: The number of delay steps used when creating the augmented system matrices A_full, B_full.

# Returns
- `u_phys_k::Vector{T}`: The computed optimal physical control input vector (dimension `controller.m_phys`).
                        Returns zeros if the QP solver fails.
- `status`: The termination status of the QP solver (e.g., `JuMP.OPTIMAL`).
"""
function get_control_input(
    controller::DenseMPCController,
    current_state_unnormalized::AbstractVector{T},
    mean_vec::AbstractVector{T},
    N_X::Int,
    N_U::Int,
    DELAY::Int
    ) where {T<:Real}

    # --- Input Validation ---
    expected_state_len = N_X * (DELAY + 1) + N_U * DELAY
    if length(current_state_unnormalized) != expected_state_len
        error("Dimension mismatch: Length of current_state_unnormalized ($(length(current_state_unnormalized))) does not match expected length N_X*(DELAY+1) + N_U*DELAY ($expected_state_len) for N_X=$N_X, N_U=$N_U, DELAY=$DELAY.")
    end
    if expected_state_len != controller.N
         error("Dimension mismatch: Calculated expected state length ($expected_state_len) does not match controller's state dimension controller.N ($(controller.N)). Check N_X, N_U, DELAY or controller setup.")
    end
     if length(mean_vec) != N_X
        error("Dimension mismatch: Length of mean_vec ($(length(mean_vec))) does not match N_X ($N_X).")
    end

    # --- State Normalization ---
    # Construct the full mean vector to subtract from the augmented state.
    # It applies to the x_k, x_{k-1}, ..., x_{k-DELAY} parts.
    # The u_{k-1}, ..., u_{k-DELAY} parts are not normalized.

    # Create a vector of zeros with the correct type and length for the input part
    zero_input_part = zeros(T, N_U * DELAY)

    # Repeat mean_vec DELAY+1 times for the state parts
    mean_state_part = repeat(mean_vec, DELAY + 1)

    # Combine the parts
    full_mean_vector = [mean_state_part; zero_input_part]

    # Perform normalization
    normalized_state_z0 = current_state_unnormalized - full_mean_vector

    # --- Compute Control Action ---
    # Call the function from the library using the pre-configured controller
    # and the now-normalized state.
    u_phys_k, status = compute_control_action(controller, normalized_state_z0)

    return u_phys_k, status
end

# --- Optional: Example demonstrating how to load and call (requires JLD2) ---
# using JLD2
# function example_run()
#     println("Running example...")
#     # --- Load Pre-computed Data ---
#     try
#         # Load controller components (assuming they were saved)
#         # This is just a placeholder - you'd load the actual controller object if saved,
#         # or reconstruct it using setup_dense_mpc as shown in the docstring.
#         # For this example, we'll load A, B etc. and call setup.
#         AB_data = JLD2.load("ABC.jld2") # Make sure this file exists
#         A_full = AB_data["A"]
#         B_full = AB_data["B"]
#         mean_vec = AB_data["mean_vec"]
#         DELAY = AB_data["delay"]
#         N_X = size(mean_vec, 1)

#         # Infer N_U (Important: Matches the preprocessed U dimension)
#         # Example inference assuming B structure B = [B0; I; 0; ...]
#         # This might need adjustment based on your exact B_full structure
#         if DELAY == 0
#             N_U = size(B_full, 2)
#         else
#             # Assuming the identity block for u(k) starts at row N_X*(DELAY+1)+1
#             # and has N_U rows. The number of columns in B_full should be N_U.
#             N_U = size(B_full, 2)
#             # Sanity check: Does the total state dim match?
#             expected_N = N_X * (DELAY + 1) + N_U * DELAY
#             if size(A_full, 1) != expected_N
#                 @warn "Inferred N_U ($N_U) might be incorrect based on A_full size vs expected size."
#             end
#         end
#         println("Inferred N_X=$N_X, N_U=$N_U, DELAY=$DELAY")

#         # Define MPC parameters (replace with your actual values)
#         m_phys = 16 # MUST match your physical system
#         Np = 5
#         T = Float64
#         Q_state_val = 1.0; R_phys_val = 0.05; Q_term_val = 10.0; target_state_val = 0.0
#         u_phys_max = 2.0

#         # Create dummy cost/constraint vectors (replace with actual ones)
#         # (Copying simplified setup logic from mpc_control.jl for demonstration)
#         N = size(A_full, 1)
#         m = size(B_full, 2)
#         Q_stage = spzeros(T, N, N); Q_stage[1, 1] = Q_state_val
#         Q_term_mat = spzeros(T, N, N); Q_term_mat[1, 1] = Q_term_val
#         R_phys_block = sparse(Diagonal(fill(R_phys_val, m_phys)))
#         R_stage = (m > m_phys) ? blockdiag(R_phys_block, spzeros(T, m - m_phys, m - m_phys)) : R_phys_block
#         q_stage = zeros(T, N); q_stage[1] = -target_state_val * Q_stage[1, 1]
#         q_term_vec = zeros(T, N); q_term_vec[1] = -target_state_val * Q_term_mat[1, 1]
#         r_stage = zeros(T, m);
#         Q_vec = AbstractMatrix{T}[copy(Q_stage) for _ = 1:Np]; push!(Q_vec, Q_term_mat);
#         R_vec = AbstractMatrix{T}[copy(R_stage) for _ = 1:Np]
#         q_vec = AbstractVector{T}[copy(q_stage) for _ = 1:Np]; push!(q_vec, q_term_vec);
#         r_vec = AbstractVector{T}[copy(r_stage) for _ = 1:Np]
#         E_stage = spzeros(T, 2 * m_phys, N)
#         F_phys_block = sparse([I(m_phys); -I(m_phys)])
#         F_stage = (m > m_phys) ? [F_phys_block spzeros(T, 2 * m_phys, m - m_phys)] : F_phys_block
#         b_stage = ones(T, 2 * m_phys) * u_phys_max
#         E_term = spzeros(T, 0, N); b_term = zeros(T, 0);
#         E_vec = AbstractMatrix{T}[copy(E_stage) for _ = 1:Np]; push!(E_vec, E_term);
#         F_vec = AbstractMatrix{T}[copy(F_stage) for _ = 1:Np]
#         b_vec = AbstractVector{T}[copy(b_stage) for _ = 1:Np]; push!(b_vec, b_term);

#         # Setup controller
#         println("Setting up controller for example...")
#         controller = setup_dense_mpc(A_full, B_full, Np, Q_vec, R_vec, q_vec, r_vec, E_vec, F_vec, b_vec, m_phys)
#         println("Controller setup complete.")

#         # --- Create Dummy State ---
#         # Create a dummy state vector matching the expected structure and size
#         dummy_state_unnormalized = randn(T, controller.N)
#         println("Created dummy state vector of size $(controller.N)")

#         # --- Call the Function ---
#         println("Calling get_control_input...")
#         u_phys_next, status = get_control_input(controller, dummy_state_unnormalized, mean_vec, N_X, N_U, DELAY)

#         println("\n--- Example Results ---")
#         println("QP Status: ", status)
#         println("Computed Physical Input (u_phys_k): ", u_phys_next)
#         println("Input Dimension: ", length(u_phys_next))

#     catch e
#         println("Error during example run: $e")
#         println("Ensure 'ABC.jld2' exists and contains A, B, mean_vec, delay.")
#         println("Ensure JLD2 package is installed (`using Pkg; Pkg.add(\"JLD2\")`).")
#     end
# end

# Uncomment to run the example:
# example_run()
