# src/controllers/alternating_mpc_ipopt.jl
using LinearAlgebra
# using OSQP # REMOVED
using JuMP
using Ipopt # NEW
using SparseArrays
using JLD2 # If loading A, B from file
using TOML # To read config inside constructor
using Logging # Added for @debug, @warn, @error

# Assuming AbstractController is defined in interface.jl
include("interface.jl")

# --- Helper Functions for Alternating Minimization (Unchanged) ---

"""
Builds the mapping matrix for a single time step theta_i = M * x_i,
where x_i is either u_i or v_i.
"""
function build_M_matrix_single(fixed_vec::Vector{Float64}, m_out::Int, n_out::Int, is_V_fixed::Bool)
    mn = m_out * n_out
    if is_V_fixed
        # Theta = vec(u * v') = (v kron I_m) * u
        @assert length(fixed_vec) == n_out "Fixed vector v must have length n=$n_out"
        # Resulting matrix maps u (m x 1) to theta (mn x 1)
        return kron(fixed_vec, sparse(I, m_out, m_out)) # Returns mn x m sparse matrix
    else
        # Theta = vec(u * v') = (I_n kron u) * v
        @assert length(fixed_vec) == m_out "Fixed vector u must have length m=$m_out"
        # Resulting matrix maps v (n x 1) to theta (mn x 1)
        return kron(sparse(I, n_out, n_out), fixed_vec) # Returns mn x n sparse matrix
    end
end

"""
Builds the block-diagonal sparse matrix M_V such that Theta = M_V * U.
"""
function build_M_V(V_sequence::Vector{Float64}, Nc::Int, m::Int, n::Int)
    @assert length(V_sequence) == n * Nc "V_sequence length mismatch"
    block_matrices = Vector{SparseMatrixCSC{Float64, Int64}}(undef, Nc)
    nu = m * n # e.g., 16
    for i = 1:Nc
        v_i = V_sequence[(i-1)*n+1 : i*n]
        block_matrices[i] = build_M_matrix_single(v_i, m, n, true) # nu x m matrix
    end
    # Ensure the result is sparse
    return blockdiag(block_matrices...)::SparseMatrixCSC{Float64, Int64} # (nu*Nc) x (m*Nc) sparse matrix
end

"""
Builds the block-diagonal sparse matrix M_U such that Theta = M_U * V.
"""
function build_M_U(U_sequence::Vector{Float64}, Nc::Int, m::Int, n::Int)
     @assert length(U_sequence) == m * Nc "U_sequence length mismatch"
    block_matrices = Vector{SparseMatrixCSC{Float64, Int64}}(undef, Nc)
    nu = m * n # e.g., 16
    for i = 1:Nc
        u_i = U_sequence[(i-1)*m+1 : i*m]
        block_matrices[i] = build_M_matrix_single(u_i, m, n, false) # nu x n matrix
    end
     # Ensure the result is sparse
    return blockdiag(block_matrices...)::SparseMatrixCSC{Float64, Int64} # (nu*Nc) x (n*Nc) sparse matrix
end


# --- MPCDataAM Struct Definition ---
# Stores data for the Alternating Minimization MPC
# MODIFIED: Removed OSQP model fields
struct MPCDataAM
    # System Matrices and Dimensions
    A::Matrix{Float64}
    B::Matrix{Float64} # B is nx x nu (e.g., nx x 16)
    nx::Int         # Augmented state dimension
    nu::Int         # Control input dimension (e.g., 16 for theta)
    m::Int          # Dimension of physical input u (e.g., 4)
    n::Int          # Dimension of physical input v (e.g., 4)
    Np::Int         # Prediction horizon
    Nc::Int         # Control horizon

    # Weighting Matrices (Original Theta-based)
    Q::Diagonal{Float64, Vector{Float64}} # nx x nx
    R::Diagonal{Float64, Vector{Float64}} # nu x nu (e.g., 16x16)

    # Reference and Constraints (Original Theta-based)
    xref::Vector{Float64} # nx
    umin::Vector{Float64} # nu (e.g., 16) - Bounds on theta
    umax::Vector{Float64} # nu (e.g., 16) - Bounds on theta
    xmin::Vector{Float64} # nx
    xmax::Vector{Float64} # nx

    # Physical Constraints for U and V sequences
    U_phys_min::Vector{Float64} # m*Nc
    U_phys_max::Vector{Float64} # m*Nc
    V_phys_min::Vector{Float64} # n*Nc
    V_phys_max::Vector{Float64} # n*Nc

    # Precomputed MPC Matrices (Original Theta-based)
    Sx::Matrix{Float64}   # (nx*Np) x nx
    Su::Matrix{Float64}   # (nx*Np) x (nu*Nc)
    Qbar::SparseMatrixCSC{Float64, Int64} # (nx*Np) x (nx*Np)
    Rbar::SparseMatrixCSC{Float64, Int64} # (nu*Nc) x (nu*Nc)
    P_sparse::SparseMatrixCSC{Float64, Int64} # Hessian for Theta QP: (nu*Nc) x (nu*Nc)
    A_con::SparseMatrixCSC{Float64, Int64}    # Constraint matrix for Theta QP (input + state)
    A_in_theta::SparseMatrixCSC{Float64, Int64} # Input part of A_con
    A_st_theta::SparseMatrixCSC{Float64, Int64} # State part of A_con

    # Stacked Constraints/References (Original Theta-based)
    l_in::Vector{Float64} # Lower bounds for theta inputs (nu*Nc)
    u_in::Vector{Float64} # Upper bounds for theta inputs (nu*Nc)
    Xmin::Vector{Float64} # Stacked state lower bounds (nx*Np)
    Xmax::Vector{Float64} # Stacked state upper bounds (nx*Np)
    Xref::Vector{Float64} # Stacked state reference (nx*Np)

    # Store constraint dimensions
    num_theta_cons::Int # Number of original theta constraints (input + state)
    num_u_phys_cons::Int # Number of physical U constraints
    num_v_phys_cons::Int # Number of physical V constraints
    num_total_cons::Int # Total constraints for OSQP setup (kept name for consistency, but not used by JuMP setup)

    # REMOVED OSQP Models
    # model_U::OSQP.Model
    # model_V::OSQP.Model
end


# --- setup_mpc_am function ---
# MODIFIED: Removed OSQP setup
function setup_mpc_am(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax, m, n,
                      U_phys_min, U_phys_max, V_phys_min, V_phys_max; osqp_settings...) # osqp_settings ignored now
    nx = size(A, 1)
    nu = size(B, 2) # Should be m*n (e.g., 16)
    @info "Setting up Alternating MPC with nx=$nx, nu=$nu (m=$m, n=$n)"
    @assert nu == m * n "Dimension mismatch: nu ($nu) != m*n ($(m*n))"

    # --- Dimension Checks (similar to setup_mpc) ---
    @assert size(A, 1) == size(A, 2) "Matrix A must be square ($nx x $nx)."
    @assert size(B, 1) == nx "Matrix B must have $nx rows (dimension mismatch with A)."
    @assert Np > 0 "Prediction horizon Np must be positive."
    @assert Nc > 0 "Control horizon Nc must be positive."
    @assert size(Q) == (nx, nx) "Matrix Q must be square with dimension $nx x $nx."
    @assert size(R) == (nu, nu) "Matrix R must be square with dimension $nu x $nu."
    @assert length(xref) == nx "Reference vector xref must have length $nx."
    @assert length(umin) == nu "Control lower bound vector umin (for theta) must have length $nu."
    @assert length(umax) == nu "Control upper bound vector umax (for theta) must have length $nu."
    @assert length(xmin) == nx "State lower bound vector xmin must have length $nx."
    @assert length(xmax) == nx "State upper bound vector xmax must have length $nx."
    @assert length(U_phys_min) == m * Nc "U_phys_min length mismatch"
    @assert length(U_phys_max) == m * Nc "U_phys_max length mismatch"
    @assert length(V_phys_min) == n * Nc "V_phys_min length mismatch"
    @assert length(V_phys_max) == n * Nc "V_phys_max length mismatch"


    # --- Dense MPC Formulation Calculations (Theta-based) ---
    Sx = zeros(nx * Np, nx)
    Su = zeros(nx * Np, nu * Nc) # nu=m*n here

    # Calculate Sx
    temp_A = A
    for i = 1:Np
        rows = (i-1)*nx+1 : i*nx
        Sx[rows, :] = temp_A
        temp_A = temp_A * A
    end

    # Calculate Su using the input Matrix B directly (B is nx x nu)
    for i = 1:Np
        rows_x = (i-1)*nx+1 : i*nx
        for j = 1:min(i, Nc)
            cols_u = (j-1)*nu+1 : j*nu # nu=m*n here
            if i - j >= 0
                Su[rows_x, cols_u] = (A^(i - j)) * B
            end
        end
    end

    # Build Theta-based QP matrices
    Qbar = kron(sparse(I(Np)), sparse(Q))
    Rbar = kron(sparse(I(Nc)), sparse(R)) # R is nu x nu

    # Regularization for P_sparse
    small_reg = 1e-6 # Default value, might need tuning
    @info "Using P_sparse regularization: $small_reg"
    Hessian_calc = Su' * Qbar * Su + Rbar + small_reg * I
    P_sparse = sparse(Symmetric(Hessian_calc)) # nu*Nc x nu*Nc
    # Check if P_sparse is positive definite
    try
        cholesky(P_sparse)
        @info "P_sparse (Theta Hessian) is positive definite."
    catch e
        @warn "P_sparse (Theta Hessian) is NOT positive definite: $e"
        # Consider increasing small_reg further if this occurs
    end


    # --- Theta-based Constraint Formulation ---
    Umin_theta = repeat(umin, Nc) # umin has length nu
    Umax_theta = repeat(umax, Nc) # umax has length nu
    A_in_theta = sparse(I, nu * Nc, nu * Nc) # Identity for direct theta bounds
    l_in = Umin_theta
    u_in = Umax_theta

    Xmin = repeat(xmin, Np)
    Xmax = repeat(xmax, Np)
    A_st_theta = sparse(Su) # State constraint matrix related to Theta

    A_con = sparse([A_in_theta; A_st_theta]) # Combined constraint matrix nu*Nc + nx*Np x nu*Nc

    # Stacked reference
    Xref = repeat(xref, Np)

    # --- Calculate Constraint Dimensions (still useful for indexing) ---
    num_theta_cons = size(A_con, 1) # Number of original constraints (input + state)
    num_u_phys_cons = m * Nc        # Number of physical U constraints
    num_v_phys_cons = n * Nc        # Number of physical V constraints
    num_total_cons = num_theta_cons + num_u_phys_cons + num_v_phys_cons
    @info "Constraint Dimensions: num_theta_cons=$num_theta_cons, num_u_phys_cons=$num_u_phys_cons, num_v_phys_cons=$num_v_phys_cons"

    # --- REMOVED OSQP Solver Setup ---

    # Create and return the MPCDataAM struct
    return MPCDataAM(A, B, nx, nu, m, n, Np, Nc, Q, R, xref, umin, umax, xmin, xmax,
                     U_phys_min, U_phys_max, V_phys_min, V_phys_max,
                     Sx, Su, Qbar, Rbar, P_sparse, A_con, A_in_theta, A_st_theta,
                     l_in, u_in, Xmin, Xmax, Xref,
                     num_theta_cons, num_u_phys_cons, num_v_phys_cons, num_total_cons)
end


# --- The AlternatingMPCController Struct (Unchanged) ---
mutable struct AlternatingMPCController <: AbstractController
    mpc_data::MPCDataAM # Use the new data structure
    sensor_indices::Union{Vector{Int}, Nothing} # Optional: Indices to map system_state to parts of xk
    mean_vec::Vector{Float64}
    last_xk::Union{Vector{Float64}, Nothing}      # Stores augmented state x_{k-1}
    last_uk::Union{Vector{Float64}, Nothing}      # Stores *effective* theta_{k-1} = vec(u_{k-1}v_{k-1}')
    input_delays::Int                             # Number of input delays in state augmentation

    # --- Alternating Minimization Specific ---
    N_iter_am::Int # Number of Alternating Minimization iterations per time step
    last_U_opt::Union{Vector{Float64}, Nothing} # Stores full U sequence [u_k,...,u_k+Nc-1] from last MPC step for warm-start
    last_V_opt::Union{Vector{Float64}, Nothing} # Stores full V sequence [v_k,...,v_k+Nc-1] from last MPC step for warm-start
end

# --- Constructor for the AlternatingMPCController (Unchanged logic, just calls modified setup) ---
function AlternatingMPCController(config::Dict)
    mpc_params = config["controller"]["params"]["AlternatingMPC"] # Use specific section in config

    matrix_file = mpc_params["matrix_file"]
    local mean_vec, A, B
    try
        jldopen(matrix_file, "r") do file
            A = file["A"]
            B = file["B"] # B should be nx x nu (e.g., nx x 16)
            mean_vec = file["mean_vec"]
            @info "Loaded A [$(size(A))] and B [$(size(B))] from $matrix_file"
        end
    catch e
        @error "Failed to load matrices from $matrix_file" exception=(e, catch_backtrace())
        rethrow(e)
    end

    # Extract other MPC parameters from config
    Np = mpc_params["Np"]
    Nc = mpc_params["Nc"]
    nx = size(A, 1) # Augmented state dimension
    nu = size(B, 2) # Should be m*n (e.g., 16)

    # --- Get AM parameters ---
    N_iter_am = get(mpc_params, "N_iter_am", 5) # Default to 5 iterations
    m_decomp = get(mpc_params, "m_decomp", 4)   # Dimension of u
    n_decomp = get(mpc_params, "n_decomp", 4)   # Dimension of v
    @info "Using Alternating Minimization with N_iter=$N_iter_am, u_dim=$m_decomp, v_dim=$n_decomp"

    if nu != m_decomp * n_decomp
        @error "Loaded B matrix nu=$(nu) does not match m*n = $(m_decomp*n_decomp). Check model/config."
        throw(ArgumentError("Dimension mismatch nu != m*n"))
    end

    # Get number of input delays
    input_delays = get(mpc_params, "input_delays", -1)
    if input_delays < 0
        @error "Missing or invalid 'input_delays' parameter in controller config."
        throw(ArgumentError("Missing 'input_delays' in config"))
    elseif input_delays * nu >= nx && input_delays > 0 # Allow input_delays = 0
         @error "Invalid delay structure: input_delays*nu ($(input_delays*nu)) >= nx ($nx). Check config/model."
         throw(ArgumentError("Invalid delay structure"))
    end
    @info "Using input_delays (du) = $input_delays"

    # Construct Q
    Q_base = Float64(get(mpc_params, "Q", 1.0))
    Q_diag_override = get(mpc_params, "Q_diag", nothing)
    local Q_diag
    if isnothing(Q_diag_override)
        Q_diag = ones(nx) * Q_base
        # TODO Override
        Q_diag[1] = Q_base
    else
        Q_diag = Float64.(Q_diag_override)
    end
    @assert length(Q_diag) == nx "Length of Q_diag must match state dimension nx=$nx"

    # Construct R (for theta)
    R_diag_vec = Float64.(mpc_params["R_diag"]) # Length should be nu (e.g., 16)
    @assert length(R_diag_vec) == nu "Length of R_diag must match input dimension nu=$nu"

    # Construct xref
    x0ref_val = Float64(mpc_params["x0ref"])
    xref = zeros(nx)
    dy = nx - nu * input_delays # Calculate number of output delays
    if dy < 0
        @error "Calculated number of output delays dy=$dy is negative. Check nx, nu, input_delays."
        throw(ArgumentError("Invalid delay structure leading to dy < 0"))
    end
    xref[1:dy] .= x0ref_val
    @info "Constructed xref with first element: $(xref[1])"

    # Theta constraints
    umin_theta = Float64.(mpc_params["umin"]) # Length nu
    umax_theta = Float64.(mpc_params["umax"]) # Length nu
    @assert length(umin_theta) == nu && length(umax_theta) == nu "Theta constraints must match nu=$nu"
    if all(umin_theta .== 0.0) && all(umax_theta .== 0.0)
        @warn "Theta constraints umin and umax are both all zeros."
    end

    # State constraints
    xmin = fill(-Inf, nx)
    xmax = fill(Inf, nx)
    # Optional: Load finite state bounds if needed from config

    # Load physical constraints for u, v
    u_phys_min_vec = Float64.(get(mpc_params, "u_phys_min", fill(-Inf, m_decomp)))
    u_phys_max_vec = Float64.(get(mpc_params, "u_phys_max", fill(+Inf, m_decomp)))
    v_phys_min_vec = Float64.(get(mpc_params, "v_phys_min", fill(-Inf, n_decomp)))
    v_phys_max_vec = Float64.(get(mpc_params, "v_phys_max", fill(+Inf, n_decomp)))

    @assert length(u_phys_min_vec) == m_decomp "u_phys_min length must match m_decomp=$m_decomp"
    @assert length(u_phys_max_vec) == m_decomp "u_phys_max length must match m_decomp=$m_decomp"
    @assert length(v_phys_min_vec) == n_decomp "v_phys_min length must match n_decomp=$n_decomp"
    @assert length(v_phys_max_vec) == n_decomp "v_phys_max length must match n_decomp=$n_decomp"
    @info "Loaded Physical U bounds: [$(minimum(u_phys_min_vec)), $(maximum(u_phys_max_vec))]"
    @info "Loaded Physical V bounds: [$(minimum(v_phys_min_vec)), $(maximum(v_phys_max_vec))]"

    # Repeat for Nc horizon
    U_phys_min = repeat(u_phys_min_vec, Nc)
    U_phys_max = repeat(u_phys_max_vec, Nc)
    V_phys_min = repeat(v_phys_min_vec, Nc)
    V_phys_max = repeat(v_phys_max_vec, Nc)

    # Ipopt settings can be passed via JuMP `set_optimizer_attribute` if needed
    # osqp_settings = get(mpc_params, "osqp_settings", Dict()) # Ignored

    @info "Setting up Alternating MPC Controller..."
    mpc_data = setup_mpc_am(A, B, Np, Nc, Diagonal(Q_diag), Diagonal(R_diag_vec), xref,
                              umin_theta, umax_theta, xmin, xmax, m_decomp, n_decomp,
                              U_phys_min, U_phys_max, V_phys_min, V_phys_max) # Removed osqp_settings
    @info "MPC setup complete."

    sensor_indices = get(mpc_params, "sensor_indices", nothing)

    # Initialize controller state
    last_xk = nothing
    last_uk = nothing
    last_U_opt = nothing
    last_V_opt = nothing

    return AlternatingMPCController(mpc_data, sensor_indices, mean_vec, last_xk, last_uk, input_delays,
                                    N_iter_am, last_U_opt, last_V_opt)
end


# --- Implement the required interface function ---
# MODIFIED: Uses JuMP with Ipopt
function compute_control_action(controller::AlternatingMPCController, time::Float64, system_state::Dict{Symbol, Any})
    # --- Extract Data ---
    mpc_data = controller.mpc_data
    (; Sx, Su, Qbar, Rbar, P_sparse, A_con, A_in_theta, A_st_theta, Xref, Xmin, Xmax, l_in, u_in,
       U_phys_min, U_phys_max, V_phys_min, V_phys_max,
       num_theta_cons, num_u_phys_cons, num_v_phys_cons, num_total_cons,
       # model_U, model_V, # REMOVED
       nx, nu, Nc, m, n, Np) = mpc_data
    N_iter_am = controller.N_iter_am
    mean_vec = controller.mean_vec
    du = controller.input_delays
    dy = nx - nu * du

    # --- State Reconstruction ---
    prev_xk = controller.last_xk === nothing ? zeros(nx) : controller.last_xk
    prev_theta = controller.last_uk === nothing ? zeros(nu) : controller.last_uk
    xk = zeros(nx)
    # ... (rest of state reconstruction logic as before) ...
    # 1. Get current measurement y_k (centered)
    local yk_centered::Float64
    measurement_key = :wallheatflux # Or load from config if variable
    if haskey(system_state, measurement_key)
        measurement_val_raw = system_state[measurement_key] ./ 1e6
        # Check if measurement_val_raw is empty or not a vector/scalar
        if !isa(measurement_val_raw, AbstractVector) && !isa(measurement_val_raw, Number)
             @error "Measurement '$measurement_key' is not a vector or number, type=$(typeof(measurement_val_raw)). Returning safe input."
             controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
        end
        if isempty(measurement_val_raw)
             @error "Measurement '$measurement_key' is an empty vector. Returning safe input."
             controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
        end
        measurement_val = isa(measurement_val_raw, AbstractVector) ? measurement_val_raw[1] : measurement_val_raw
        @info "Measurement '$measurement_key' value: $measurement_val"
        if length(mean_vec) >= 1
            yk_centered = Float64(measurement_val) - mean_vec[1]
        else
            @error "mean_vec missing or too short, cannot center measurement. Returning safe input."
            controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
        end
    else
        @error "Missing measurement key '$measurement_key' in system_state. Returning safe input."
        controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
    end
    xk[1] = yk_centered # y_k
    # 2. Past outputs: y_{k-1} down to y_{k-dy+1}
    if dy > 1
        if length(prev_xk) >= (dy - 1)
            xk[2:dy] = prev_xk[1:dy-1]
        else
            @error "prev_xk too short for output history reconstruction. length=$(length(prev_xk)), needed=$(dy-1). Returning safe input."
            controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
        end
    elseif dy < 1
         @warn "Number of output delays dy=$dy is not positive. State reconstruction might be incorrect."
    end
    # 3. Past inputs: theta_{k-1} down to theta_{k-du}
    if du > 0
        start_idx_theta_km1 = dy + 1
        end_idx_theta_km1 = dy + nu
        if length(prev_theta) == nu && end_idx_theta_km1 <= nx && start_idx_theta_km1 > 0
            xk[start_idx_theta_km1 : end_idx_theta_km1] = prev_theta
        else
            @error "prev_theta length mismatch or index out of bounds for theta_{k-1}. length=$(length(prev_theta)), nu=$nu, start_idx=$start_idx_theta_km1, end_idx=$end_idx_theta_km1, nx=$nx. Returning safe input."
            controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
        end
        if du > 1
            start_idx_prev_thetas = dy + 1
            num_older_thetas = (du - 1) * nu
            end_idx_prev_thetas = start_idx_prev_thetas + num_older_thetas - 1
            start_idx_new_thetas = end_idx_theta_km1 + 1
            end_idx_new_thetas = start_idx_new_thetas + num_older_thetas - 1
            if start_idx_prev_thetas > 0 && end_idx_prev_thetas <= length(prev_xk) && start_idx_new_thetas > 0 && end_idx_new_thetas <= length(xk)
                 xk[start_idx_new_thetas : end_idx_new_thetas] = prev_xk[start_idx_prev_thetas : end_idx_prev_thetas]
            else
                 @error "Index out of bounds during state assembly for older thetas. Check dimensions/delays. Returning safe input."
                 controller.last_xk = xk; controller.last_U_opt = zeros(m * Nc); controller.last_V_opt = zeros(n * Nc); controller.last_uk = zeros(nu); return zeros(m + n)
            end
        end
    end
    @debug "Constructed state xk (norm, first 5): $(round(norm(xk), digits=4)), $(round.(xk[1:min(5, length(xk))], digits=4))"

    # --- Calculate State-Dependent QP Terms (Unchanged) ---
    state_error = Sx * xk - Xref
    q_theta_k = Su' * Qbar * state_error
    l_st_k = Xmin - Sx * xk
    u_st_k = Xmax - Sx * xk
    l_theta_k = [l_in; l_st_k]
    u_theta_k = [u_in; u_st_k]
    @debug "State bounds check: any(l_st_k .> u_st_k) = $(any(l_st_k .> u_st_k))"
    # ... (debug prints for bounds if needed) ...

    # --- Initialize U and V for AM (Unchanged) ---
    init_val = 0.01
    warm_start_threshold = 1e-6
    if controller.last_U_opt !== nothing && length(controller.last_U_opt) == m * Nc && norm(controller.last_U_opt) > warm_start_threshold
        U_current = [controller.last_U_opt[m+1 : m*Nc]; zeros(m)]
    else
        U_current = rand(m * Nc) * init_val
    end
    if controller.last_V_opt !== nothing && length(controller.last_V_opt) == n * Nc && norm(controller.last_V_opt) > warm_start_threshold
        V_current = [controller.last_V_opt[n+1 : n*Nc]; zeros(n)]
    else
        V_current = rand(m * Nc) * init_val
    end

    # --- Alternating Minimization Loop ---
    solve_success = true
    qp_regularization = 1e-4

    # Pre-allocate sparse identity matrices
    I_u_bounds_template = sparse(I, num_u_phys_cons, m * Nc)
    I_v_bounds_template = sparse(I, num_v_phys_cons, n * Nc)

    local U_current_final, V_current_final
    U_current_final = copy(U_current) # Keep track of the final U
    V_current_final = copy(V_current) # Keep track of the final V


    for iter = 1:N_iter_am
        #@debug "AM Iteration $iter/$N_iter_am"

        # --- 1. Optimize for U (given V_current) ---
        try
            M_V = build_M_V(V_current, Nc, m, n)
            # REMOVED: Check for near-zero V_current / M_V
            # if norm(V_current) < 1e-9 || norm(M_V) < 1e-9
            #      @warn "AM Iter $iter: V_current or M_V matrix is near zero. Skipping U optimization."
            # else
            # Calculate U-QP Hessian and Gradient
            P_U_dense = M_V' * P_sparse * M_V
            P_U = sparse(0.5 * (P_U_dense + P_U_dense')) # Ensure symmetry
            P_U += (qp_regularization * I)

            q_U = M_V' * q_theta_k

            # Map original theta constraints to U space
            A_con_U = A_con * M_V

            # Build Augmented Constraints (Matrix and Bounds)
            Z_v_bounds_for_u_step = spzeros(num_v_phys_cons, m * Nc)
            A_U_aug = sparse([A_con_U; I_u_bounds_template; Z_v_bounds_for_u_step])
            l_U_aug = [l_theta_k; U_phys_min; fill(-Inf, num_v_phys_cons)]
            u_U_aug = [u_theta_k; U_phys_max; fill(+Inf, num_v_phys_cons)]

            # --- Solve U-Step using JuMP + Ipopt ---
            model_u = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model_u, "print_level", 0) # Suppress Ipopt output
            # Add other Ipopt attributes if needed, e.g., tolerances:
            # set_optimizer_attribute(model_u, "tol", 1e-6)
            # set_optimizer_attribute(model_u, "constr_viol_tol", 1e-6)

            @variable(model_u, U_vars[1:(m*Nc)])
            # Set warm start if available
            set_start_value.(U_vars, U_current)

            @objective(model_u, Min, 0.5 * U_vars' * P_U * U_vars + q_U' * U_vars)

            # Add constraints row by row, handling infinities
            for i in 1:num_total_cons
                # Skip if row is all zeros (can happen with sparse structure)
                if nnz(A_U_aug[i,:]) == 0
                    continue
                end
                row_expr = @expression(model_u, dot(A_U_aug[i,:], U_vars))
                # row_expr = @expression(model_u, sum(A_U_aug[i, j] * U_vars[j] for j in 1:(m*Nc))) # Alternative slower way
                if isfinite(l_U_aug[i]) && isfinite(u_U_aug[i])
                     @constraint(model_u, l_U_aug[i] <= row_expr <= u_U_aug[i])
                elseif isfinite(l_U_aug[i])
                    @constraint(model_u, row_expr >= l_U_aug[i])
                elseif isfinite(u_U_aug[i])
                    @constraint(model_u, row_expr <= u_U_aug[i])
                end
                 # If both are infinite, no constraint needed for this row
            end

            optimize!(model_u)

            term_status = termination_status(model_u)
            prim_status = primal_status(model_u)

            if term_status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED) && prim_status == MOI.FEASIBLE_POINT
                U_current = value.(U_vars)
                if !(term_status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED))
                     @warn "AM Iter $iter: Ipopt for U solved with approximate status $term_status."
                end
                #@debug "  U-Step solved. New U_current norm: $(round(norm(U_current), digits=4))"
            else
                @warn "AM Iter $iter: Ipopt for U failed. Termination: $term_status, Primal: $prim_status. Keeping previous U."
                # Log objective value if available
                if has_values(model_u)
                     @warn "  Ipopt U-step objective value: $(objective_value(model_u))"
                end
                solve_success = false
                break # Exit AM loop on failure
            end
            # --- End U-Step Solve ---
            # end # REMOVED: end for the if norm(...) check
        catch e
             @error "AM Iter $iter: Error during U optimization (JuMP/Ipopt): $e" stacktrace=stacktrace(catch_backtrace())
             solve_success = false
             break
        end

        # --- 2. Optimize for V (given U_current) ---
        if !solve_success break end # Check again before V-step

        try
            M_U = build_M_U(U_current, Nc, m, n)
            #  if norm(U_current) < 1e-9 || norm(M_U) < 1e-9
            #      @warn "AM Iter $iter: U_current or M_U matrix is near zero. Skipping V optimization."
            #  else
            # Calculate V-QP Hessian and Gradient
            P_V_dense = M_U' * P_sparse * M_U
            P_V = sparse(0.5 * (P_V_dense + P_V_dense')) # Ensure symmetry
            P_V += (qp_regularization * I)

            q_V = M_U' * q_theta_k

            # Map original theta constraints to V space
            A_con_V = A_con * M_U

            # Build Augmented Constraints (Matrix and Bounds)
            Z_u_bounds_for_v_step = spzeros(num_u_phys_cons, n * Nc)
            A_V_aug = sparse([A_con_V; Z_u_bounds_for_v_step; I_v_bounds_template])
            l_V_aug = [l_theta_k; fill(-Inf, num_u_phys_cons); V_phys_min]
            u_V_aug = [u_theta_k; fill(+Inf, num_u_phys_cons); V_phys_max]

            # --- Solve V-Step using JuMP + Ipopt ---
            model_v = Model(Ipopt.Optimizer)
            set_optimizer_attribute(model_v, "print_level", 0)

            @variable(model_v, V_vars[1:(n*Nc)])
            set_start_value.(V_vars, V_current) # Warm start

            @objective(model_v, Min, 0.5 * V_vars' * P_V * V_vars + q_V' * V_vars)

            # Add constraints row by row
            for i in 1:num_total_cons
                 if nnz(A_V_aug[i,:]) == 0
                    continue
                end
                row_expr = @expression(model_v, dot(A_V_aug[i,:], V_vars))
                # row_expr = @expression(model_v, sum(A_V_aug[i, j] * V_vars[j] for j in 1:(n*Nc)))
                if isfinite(l_V_aug[i]) && isfinite(u_V_aug[i])
                     @constraint(model_v, l_V_aug[i] <= row_expr <= u_V_aug[i])
                elseif isfinite(l_V_aug[i])
                    @constraint(model_v, row_expr >= l_V_aug[i])
                elseif isfinite(u_V_aug[i])
                    @constraint(model_v, row_expr <= u_V_aug[i])
                end
            end

            optimize!(model_v)

            term_status = termination_status(model_v)
            prim_status = primal_status(model_v)

            if term_status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL, MOI.ALMOST_LOCALLY_SOLVED) && prim_status == MOI.FEASIBLE_POINT
                V_current = value.(V_vars)
                 if !(term_status in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED))
                     @warn "AM Iter $iter: Ipopt for V solved with approximate status $term_status."
                end
                #@debug "  V-Step solved. New V_current norm: $(round(norm(V_current), digits=4))"
            else
                @warn "AM Iter $iter: Ipopt for V failed. Termination: $term_status, Primal: $prim_status. Keeping previous V."
                 if has_values(model_v)
                     @warn "  Ipopt V-step objective value: $(objective_value(model_v))"
                end
                solve_success = false
                break
            end
            # --- End V-Step Solve ---
            #  end # REMOVED: end for the if norm(...) check
        catch e
             @error "AM Iter $iter: Error during V optimization (JuMP/Ipopt): $e" stacktrace=stacktrace(catch_backtrace())
             solve_success = false
             break
        end

        U_current_final = copy(U_current)
        V_current_final = copy(V_current)
    end # End AM loop

    # --- >>> START: Calculate Expected Trajectory <<< ---
    local predicted_states::Matrix{Float64} # Define type for clarity
    if solve_success
        @info "AM converged. Calculating predicted trajectory."

        # 1. Reconstruct the optimal Theta sequence from final U and V
        Theta_opt = zeros(Float64, nu * Nc)
        for i = 1:Nc
            # Extract u_i and v_i for time step k+i-1
            u_i_idx = (i-1)*m+1 : i*m
            v_i_idx = (i-1)*n+1 : i*n

            # Ensure indices are within bounds (safety check)
            if maximum(u_i_idx) > length(U_current_final) || maximum(v_i_idx) > length(V_current_final)
                 @error "Index out of bounds during Theta_opt reconstruction. i=$i, Nc=$Nc"
                 # Handle error appropriately, maybe set solve_success = false?
                 solve_success = false
                 break # Exit the Theta_opt calculation
            end

            u_i = U_current_final[u_i_idx]
            v_i = V_current_final[v_i_idx]

            # Calculate theta_i = vec(u_i * v_i')
            theta_i = vec(u_i * v_i')

            # Place it in the stacked vector
            theta_opt_idx = (i-1)*nu+1 : i*nu
            if maximum(theta_opt_idx) > length(Theta_opt)
                 @error "Index out of bounds for Theta_opt placement. i=$i, Nc=$Nc"
                 solve_success = false
                 break
            end
            Theta_opt[theta_opt_idx] = theta_i
        end

        if solve_success # Recalculate only if Theta_opt was built successfully
             # 2. Calculate the stacked predicted state vector
             X_pred_stacked = Sx * xk + Su * Theta_opt # Use mpc_data.Sx and mpc_data.Su

             # 3. Reshape into a matrix (nx rows, Np columns)
             # Each column is a predicted state x_{k+i|k}
             predicted_states = reshape(X_pred_stacked, nx, Np)

             # 4. Log or store the result
             # Example: Log the norm and the first predicted state vector
             @info "Predicted Trajectory (norm): $(round(norm(predicted_states), digits=4))"
             @info "Predicted x_{k+1|k} (first 5): $(round.(predicted_states[1:min(5, nx), 1], digits=4))"
             @info "Predicted x_{k+Np|k} (first 5): $(round.(predicted_states[1:min(5, nx), Np], digits=4))"

             # --- TODO: Store predicted_states if needed for external analysis ---
             # Option 1: Add a field to the controller struct
             # controller.last_predicted_trajectory = predicted_states
             # Option 2: Return it along with the control action (requires changing the function signature and caller)
             # return physical_input, predicted_states
        end
    end
    # --- >>> END: Calculate Expected Trajectory <<< ---


    # --- Extract Final Control Action ---
    if !solve_success
        @error "Alternating minimization failed to converge or encountered an error. Returning zero input."
        controller.last_xk = xk
        controller.last_U_opt = zeros(m * Nc)
        controller.last_V_opt = zeros(n * Nc)
        controller.last_uk = zeros(nu)
        return zeros(m + n)
    end

    uk_final = U_current[1:m]
    vk_final = V_current[1:n]

    # --- Store State and Input for Next Step (Unchanged) ---
    controller.last_xk = xk
    controller.last_U_opt = U_current
    controller.last_V_opt = V_current
    theta_k_effective = vec(uk_final * vk_final')
    if length(theta_k_effective) != nu
        @error "Internal error: Effective theta dimension mismatch. Expected $nu, got $(length(theta_k_effective))."
         controller.last_uk = zeros(nu)
    else
         controller.last_uk = theta_k_effective
    end

    # --- Return Physical Input (Unchanged) ---
    physical_input = [uk_final; vk_final]
    expected_physical_dim = m + n
    if any(isnan, physical_input)
        @error "NaN detected in final control input. Returning zeros."
        controller.last_U_opt = zeros(m * Nc)
        controller.last_V_opt = zeros(n * Nc)
        controller.last_uk = zeros(nu)
        return zeros(expected_physical_dim)
    end
     if length(physical_input) != expected_physical_dim
         @warn "Final physical input length incorrect ($(length(physical_input)) vs $expected_physical_dim). Returning zeros."
         return zeros(expected_physical_dim)
     end

    @info "Sending input: $(round.(physical_input, digits=4))"
    return physical_input
end

# --- update_controller_state! function (Unchanged) ---
function update_controller_state!(controller::AlternatingMPCController, time::Float64, new_state_data::Dict{Symbol, Any})
    # ... (no changes needed) ...
    return nothing
end
