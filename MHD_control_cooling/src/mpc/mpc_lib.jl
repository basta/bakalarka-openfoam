module DenseMPCLib

using LinearAlgebra
using BlockArrays
using SparseArrays
using JuMP
using OSQP
using Logging

export DenseMPCController, setup_dense_mpc, compute_control_action, create_delay_state

#-----------------------------------------------------------------------------
# Helper Function (Internal) - Computes Dense Matrices
#-----------------------------------------------------------------------------
"""
    compute_dense_mpc_matrices(A, B, Np, Q, R, q, r, E, F, b)

Computes the matrices for the dense Model Predictive Control (MPC)
Quadratic Programming (QP) formulation based on the Korda & Mezić paper
(Automatica 2018, Appendix).

Internal function used by setup_dense_mpc.
(Docstring adapted from previous version for context)
"""
function compute_dense_mpc_matrices(
    A::AbstractMatrix{T}, B::AbstractMatrix{T}, Np::Int,
    Q::Vector{<:AbstractMatrix{T}}, R::Vector{<:AbstractMatrix{T}},
    q::Vector{<:AbstractVector{T}}, r::Vector{<:AbstractVector{T}},
    E::Vector{<:AbstractMatrix{T}}, F::Vector{<:AbstractMatrix{T}},
    b::Vector{<:AbstractVector{T}}
) where {T<:Real}

    N = size(A, 1) # Dimension of lifted state
    m = size(B, 2) # Dimension of model input vector (potentially large)

    # --- Input Validation (Basic) ---
    @assert length(Q) == Np + 1 "Length of Q must be Np+1"
    @assert length(R) == Np "Length of R must be Np"
    @assert length(q) == Np + 1 "Length of q must be Np+1"
    @assert length(r) == Np "Length of r must be Np"
    @assert length(E) == Np + 1 "Length of E must be Np+1"
    @assert length(F) == Np "Length of F must be Np"
    @assert length(b) == Np + 1 "Length of b must be Np+1"

    # --- 1. Construct Block System Matrices A_bl (script A) and B_bl (script B) ---
    A_bl = zeros(T, N * (Np + 1), N)
    B_bl = zeros(T, N * (Np + 1), m * Np)

    A_pow_i = Matrix{T}(I, N, N) # A^0 = Identity
    A_bl[1:N, :] = A_pow_i

    for i = 1:Np
        A_pow_i = A * A_pow_i # A^i
        A_bl[i*N+1:(i+1)*N, :] = A_pow_i

        # Construct row i of B_bl (corresponding to z_i)
        term = B
        B_bl[i*N+1:(i+1)*N, (i-1)*m+1:i*m] = term # A^(0)B u_{i-1} term
        A_pow_j = A
        for j = 1:(i-1)
            term = A * term # A^j * B
            B_bl[i*N+1:(i+1)*N, (i-1-j)*m+1:(i-j)*m] = term # A^(j)B u_{i-1-j} term
        end
    end

    # --- 2. Construct Block Cost Matrices Q_bl, R_bl, q_bl, r_bl ---
    Q_bl = blockdiag([sparse(Q[i+1]) for i in 0:Np]...)
    R_bl = blockdiag([sparse(R[i+1]) for i in 0:Np-1]...)
    q_bl = vcat(q...)
    r_bl = vcat(r...)

    # --- 3. Construct Block Constraint Matrices E_bl, F_bl, b_bl ---
    # Determine total number of constraints and check consistency
    num_constraints_per_step = [size(E[i+1], 1) for i in 0:Np]
    total_constraints = sum(num_constraints_per_step)

    # Pad F with zeros corresponding to terminal constraint E[Np+1]
    F_list_padded = [sparse(F[i+1]) for i in 0:Np-1]
    push!(F_list_padded, spzeros(T, num_constraints_per_step[end], m))

    E_bl = blockdiag([sparse(E[i+1]) for i in 0:Np]...)
    F_bl_sparse = blockdiag(F_list_padded...) # Aligns with Z = [z0', ..., zNp']' and U_padded = [u0', ..., uNp']'

    # Extract the part corresponding to U = [u0', ..., u{Np-1}']'
    F_bl = F_bl_sparse[:, 1:(m*Np)]

    b_bl = vcat(b...)

    # --- 4. Compute Dense QP Matrices H, G, h, L, M, c ---
    H = R_bl + B_bl' * Q_bl * B_bl
    G = 2 * A_bl' * Q_bl * B_bl
    h = r_bl + B_bl' * q_bl

    L = E_bl * B_bl + F_bl
    M = E_bl * A_bl
    c = b_bl

    # Ensure H is numerically symmetric
    H = (Matrix(H) + Matrix(H)') / 2 # Convert to dense for factorizations

    # Convert output matrices to dense for typical QP solver interfaces
    # Note: This can consume significant memory if m*Np is large
    return Matrix(H), Matrix(G), Vector(h), Matrix(L), Matrix(M), Vector(c)
end


#-----------------------------------------------------------------------------
# Controller Struct
#-----------------------------------------------------------------------------
"""
    DenseMPCController{T<:Real}

Holds the pre-computed matrices and parameters for a dense MPC controller.

Fields:
- H, G, h, L, M, c: Dense QP matrices.
- Np: Prediction horizon.
- N: Lifted state dimension.
- m: Model input dimension (including delays).
- m_phys: Physical input dimension.
- qp_model: Pre-allocated JuMP model for the QP.
- U_var: JuMP optimization variable reference.
- cons_ref: JuMP constraint reference.
"""
struct DenseMPCController{T<:Real}
    H::Matrix{T}
    G::Matrix{T}
    h::Vector{T}
    L::Matrix{T}
    M::Matrix{T}
    c::Vector{T}
    Np::Int
    N::Int
    m::Int
    m_phys::Int
    qp_model::Model
    U_var::Vector{VariableRef}
    cons_ref::Vector{ConstraintRef}
end

#-----------------------------------------------------------------------------
# Setup Function (Offline Computation)
#-----------------------------------------------------------------------------
"""
    setup_dense_mpc(A, B, Np, Q, R, q, r, E, F, b, m_phys; solver_time_limit=0.1)

Performs the offline computation for dense MPC and returns a controller object.

Arguments:
- A, B: Lifted system matrices (N x N, N x m).
- Np: Prediction horizon.
- Q, R, q, r: Cost function matrices/vectors (Vectors of length Np+1 or Np).
- E, F, b: Constraint matrices/vectors (Vectors of length Np+1 or Np).
- m_phys: Dimension of the physical control input.
- solver_time_limit: Time limit in seconds for the online QP solve.

Returns:
- An instance of `DenseMPCController`.
"""
function setup_dense_mpc(
    A::AbstractMatrix{T}, B::AbstractMatrix{T}, Np::Int,
    Q::Vector{<:AbstractMatrix{T}}, R::Vector{<:AbstractMatrix{T}},
    q::Vector{<:AbstractVector{T}}, r::Vector{<:AbstractVector{T}},
    E::Vector{<:AbstractMatrix{T}}, F::Vector{<:AbstractMatrix{T}},
    b::Vector{<:AbstractVector{T}},
    m_phys::Int;
    solver_time_limit::Float64=0.1 # Default time limit 100ms
) where {T<:Real}

    N = size(A, 1)
    m = size(B, 2)
    mNp = m * Np

    @info "Setting up Dense MPC Controller (Np=$Np, N=$N, m=$m, m_phys=$m_phys)..."
    @info "Computing dense matrices..."
    # --- Call internal function for offline computation ---
    H_dense, G_dense, h_vec, L_dense, M_dense, c_vec = @time compute_dense_mpc_matrices(A, B, Np, Q, R, q, r, E, F, b)
    @info "Dense matrices computed."

    if mNp > 5000 # Heuristic threshold
        @warn """Dense formulation Hessian H has dimension $(size(H_dense)).
                 Solving a QP of this size online might be very slow or infeasible,
                 even with pre-computation. Consider sparse MPC if performance is inadequate."""
    end

    # --- Pre-allocate JuMP model ---
    @info "Pre-allocating JuMP QP model..."
    qp_model = Model(() -> OSQP.Optimizer()) # Use factory for fresh solver state
    set_silent(qp_model)

    @variable(qp_model, U[1:mNp])
    # Define objective structure (linear part is placeholder)
    @objective(qp_model, Min, 0.5 * U' * H_dense * U + h_vec' * U) # Use h_vec as placeholder linear part
    # Define constraint structure (RHS is placeholder)
    @constraint(qp_model, cons[j=1:size(L_dense, 1)], (L_dense*U)[j] <= c_vec[j]) # Use c_vec as placeholder RHS

    @info "Dense MPC Controller setup complete."

    return DenseMPCController{T}(
        H_dense, G_dense, h_vec, L_dense, M_dense, c_vec,
        Np, N, m, m_phys,
        qp_model, U, cons # Store references to model, variable, constraint
    )
end

#-----------------------------------------------------------------------------
# Control Action Function (Online Computation)
#-----------------------------------------------------------------------------
"""
    compute_control_action(controller::DenseMPCController, z0::AbstractVector)

Computes the optimal physical control action using the pre-computed dense MPC controller.

Arguments:
- controller: An instance of `DenseMPCController`.
- z0: The current lifted state vector (dimension N).

Returns:
- u_phys_k: The optimal physical control input vector (dimension m_phys). Returns zeros if QP fails.
- status: The termination status of the QP solver.
"""
function compute_control_action(
    controller::DenseMPCController{T},
    z0::AbstractVector{T}
) where {T<:Real}

    @assert length(z0) == controller.N "Dimension mismatch: length(z0) != controller.N"

    # --- Update QP parameters based on current state z0 ---
    linear_cost_vec = controller.h + controller.G' * z0
    constraint_rhs = controller.c - controller.M * z0

    # --- Update JuMP model efficiently ---
    set_objective_coefficient.(controller.qp_model, controller.U_var, linear_cost_vec)
    set_normalized_rhs.(controller.cons_ref, constraint_rhs)

    # --- Solve the QP ---
    optimize!(controller.qp_model)
    status = termination_status(controller.qp_model)

    # --- Get control input ---
    u_prime_optimal_k0 = zeros(T, controller.m) # Initialize model input u'_0 (size N_U = 16)
    if status == JuMP.OPTIMAL || status == JuMP.ALMOST_OPTIMAL
        U_optimal = value.(controller.U_var)
        # Extract the *model* input for the current step (first m elements of U_optimal)
        u_prime_optimal_k0 .= U_optimal[1:controller.m] # <--- MODIFIED: Get full u'_0 (size m=N_U=16)
    else
        @warn "MPC QP solve failed. Status: $status. Returning zero model input."
        # u_prime_optimal_k0 remains zeros
    end

    # Note: No clamping applied here, as this is the model input u'.
    # Clamping should happen on the physical u_e, u_m after projection if needed.

    @info "mpc_lib returning optimal model input u'_0 size $(size(u_prime_optimal_k0))" # Modified log
    return u_prime_optimal_k0, status # <--- MODIFIED: Return u'_0 and status
end


"""
    initialize_system_state_and_history(z_hist_real, u_phys_hist_real, nu, N, m_phys)

Initializes the state vector z0 and the input history b# CRITICAL ASSUMPTIONS:
1. State Initialization: Assumes the last column of `z_hist_real` (`z_hist_real[:, end]`)
   represents the desired initial state vector `z0`. This could come from a previous
   simulation, state estimation, or direct measurement if the full state is available.
2. Input History: Assumes the `input_history` buffer should store past physical inputs
   `[u_{phys,k-1}, u_{phys,k-2}, ..., u_{phys,k-nu}]^T` used by the simulation loop,
   extracted from `u_phys_hist_real`.

# Arguments
- `z_hist_real::Matrix{T}`: Real state history, shape (N, steps), ending with column z_0.
                            Needs at least 1 column.
- `u_phys_hist_real::Matrix{T}`: Real physical input history, shape (m_phys, steps),
                                 ending with column for u_{phys,-1}. Needs at least `nu` columns.
- `nu::Int`: Number of past physical inputs (N_input_delays).
- `N::Int`: Total dimension of the state vector `z`.
- `m_phys::Int`: Dimension of the physical input vector.

# Returns
- `z0::Vector{T}`: The initialized state vector (size N).
- `init_input_hist::Matrix{T}`: The initialized input history buffer (size m_phys x nu).
"""

function create_delay_state(X_data::Matrix{T}, U_data::Matrix{T}, delay::Int) where {T<:Number}
    # --- Input Validation ---
    if delay < 0
        throw(ArgumentError("delay must be non-negative, got $delay"))
    end

    num_snapshots_x = size(X_data, 2)
    num_snapshots_u = size(U_data, 2)
    required_x_cols = delay + 1
    required_u_cols = delay # Need columns up to index `end-1` for u_{k-1}, down to `end-delay` for u_{k-delay}

    if num_snapshots_x < required_x_cols
        throw(DimensionMismatch("X_data requires at least $required_x_cols columns for delay=$delay, but has only $num_snapshots_x"))
    end
    # Only check U_data columns if delay > 0, as no U columns are needed for delay=0
    if delay > 0 && num_snapshots_u < required_u_cols
        throw(DimensionMismatch("U_data requires at least $required_u_cols columns for delay=$delay, but has only $num_snapshots_u"))
    end

    # --- Extract Data Slices ---

    # Extract x_k, x_{k-1}, ..., x_{k-delay}
    # Columns: end, end-1, ..., end-delay
    x_slice = X_data[:, end:-1:end-delay] # Shape (state_dim, delay+1)

    # Extract u_{k-1}, u_{k-2}, ..., u_{k-delay}
    # Columns: end-1, end-2, ..., end-delay
    local u_slice::Matrix{T}
    if delay > 0
        u_slice = U_data[:, end-1:-1:end-delay] # Shape (input_dim, delay)
    else
        # Create an empty matrix with correct element type if delay is 0
        input_dim = size(U_data, 1) # Get input dim even if not used, for type stability? Or assume non-empty U_data?
        # Let's create 0 columns, input_dim rows if possible
        u_slice = Matrix{T}(undef, size(U_data, 1), 0) # Shape (input_dim, 0)
    end

    # --- Concatenate ---
    # Flatten the slices column by column and concatenate
    # vec() flattens column-major, which is exactly what we want:
    # [col1; col2; ...] = [x_k; x_{k-1}; ...]
    z_k = vcat(vec(x_slice), vec(u_slice))

    return z_k
end


end # end module DenseMPCLib
