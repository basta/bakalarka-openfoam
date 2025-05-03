# src/controllers/alternating_mpc.jl
using LinearAlgebra
using OSQP
using SparseArrays
using JLD2 # If loading A, B from file
using TOML # To read config inside constructor
using Logging # Added for @debug, @warn, @error

# Assuming AbstractController is defined in interface.jl
include("interface.jl")

# --- Helper Functions for Alternating Minimization ---

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

    # NEW: Physical Constraints for U and V sequences
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

    # MODIFIED: Store constraint dimensions
    num_theta_cons::Int # Number of original theta constraints (input + state)
    num_u_phys_cons::Int # Number of physical U constraints
    num_v_phys_cons::Int # Number of physical V constraints
    num_total_cons::Int # Total constraints for OSQP setup

    # OSQP Models for Alternating Minimization
    model_U::OSQP.Model # OSQP model for optimizing U (m*Nc variables)
    model_V::OSQP.Model # OSQP model for optimizing V (n*Nc variables)
end


# --- setup_mpc_am function ---
# MODIFIED: Added physical bounds arguments
function setup_mpc_am(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax, m, n,
                      U_phys_min, U_phys_max, V_phys_min, V_phys_max; osqp_settings...)
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
    # NEW: Check physical bounds dimensions
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

    # MODIFIED: Increase regularization for P_sparse
    small_reg = 1e-6 # Increased from 1e-8
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

    # --- OSQP Solver Setup for U and V ---
    num_vars_u = m * Nc
    num_vars_v = n * Nc

    # MODIFIED: Calculate constraint dimensions including physical bounds
    num_theta_cons = size(A_con, 1) # Number of original constraints (input + state)
    num_u_phys_cons = m * Nc        # Number of physical U constraints
    num_v_phys_cons = n * Nc        # Number of physical V constraints
    num_total_cons = num_theta_cons + num_u_phys_cons + num_v_phys_cons # Total rows for OSQP setup
    @info "OSQP Setup: num_theta_cons=$num_theta_cons, num_u_phys_cons=$num_u_phys_cons, num_v_phys_cons=$num_v_phys_cons, num_total_cons=$num_total_cons"

    # Default OSQP settings
    # Enable OSQP scaling, increase max_iter
    # MODIFIED: Set verbose=true for debugging
    default_settings = Dict(:verbose => true, :eps_abs => 1e-4, :eps_rel => 1e-4, :max_iter => 20000, :check_termination => 25, :adaptive_rho => true, :scaling => 10) # Increased max_iter
    settings = merge(default_settings, Dict(osqp_settings))
    @info "Using OSQP settings: $settings"


    # MODIFIED: Placeholder matrices for setup (use num_total_cons)
    P_placeholder_u = spzeros(num_vars_u, num_vars_u)
    q_placeholder_u = zeros(num_vars_u)
    A_placeholder_u = spzeros(num_total_cons, num_vars_u) # Use total constraint size

    P_placeholder_v = spzeros(num_vars_v, num_vars_v)
    q_placeholder_v = zeros(num_vars_v)
    A_placeholder_v = spzeros(num_total_cons, num_vars_v) # Use total constraint size

    # Use +/- Inf for placeholder bounds where constraints aren't active initially
    l_placeholder = fill(-Inf, num_total_cons) # Use total constraint size
    u_placeholder = fill(+Inf, num_total_cons) # Use total constraint size

    # Setup model_U
    model_U = OSQP.Model()
    try
        OSQP.setup!(model_U; P=P_placeholder_u, q=q_placeholder_u, A=A_placeholder_u,
                      l=l_placeholder, u=u_placeholder, settings...)
    catch e
        @error "Failed to setup OSQP model_U: $e"
        rethrow(e)
    end

    # Setup model_V
    model_V = OSQP.Model()
     try
        OSQP.setup!(model_V; P=P_placeholder_v, q=q_placeholder_v, A=A_placeholder_v,
                      l=l_placeholder, u=u_placeholder, settings...)
    catch e
        @error "Failed to setup OSQP model_V: $e"
        rethrow(e)
    end

    # Create and return the MPCDataAM struct
    # MODIFIED: Store physical bounds and constraint dimensions
    return MPCDataAM(A, B, nx, nu, m, n, Np, Nc, Q, R, xref, umin, umax, xmin, xmax,
                     U_phys_min, U_phys_max, V_phys_min, V_phys_max, # NEW fields
                     Sx, Su, Qbar, Rbar, P_sparse, A_con, A_in_theta, A_st_theta,
                     l_in, u_in, Xmin, Xmax, Xref,
                     num_theta_cons, num_u_phys_cons, num_v_phys_cons, num_total_cons, # NEW fields
                     model_U, model_V)
end


# --- The AlternatingMPCController Struct ---
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

    # Physical constraints are now handled inside MPCDataAM and setup
end

# --- Constructor for the AlternatingMPCController ---
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
    else
        Q_diag = Float64.(Q_diag_override)
    end
    @assert length(Q_diag) == nx "Length of Q_diag must match state dimension nx=$nx"
    # @debug "Q diag min/max: $(minimum(Q_diag)) / $(maximum(Q_diag))"


    # Construct R (for theta)
    R_diag_vec = Float64.(mpc_params["R_diag"]) # Length should be nu (e.g., 16)
    @assert length(R_diag_vec) == nu "Length of R_diag must match input dimension nu=$nu"
    # @debug "R diag min/max: $(minimum(R_diag_vec)) / $(maximum(R_diag_vec))"


    # Construct xref
    x0ref_val = Float64(mpc_params["x0ref"])
    xref = zeros(nx)
    dy = nx - nu * input_delays # Calculate number of output delays (assuming nu = m*n used for delays)
    if dy < 0
        @error "Calculated number of output delays dy=$dy is negative. Check nx, nu, input_delays."
        throw(ArgumentError("Invalid delay structure leading to dy < 0"))
    end
    xref[1:dy] .= x0ref_val
    @info "Constructed xref with first element: $(xref[1])"

    # Theta constraints (will be transformed for U/V QPs)
    umin_theta = Float64.(mpc_params["umin"]) # Length should be nu (e.g., 16)
    umax_theta = Float64.(mpc_params["umax"]) # Length should be nu (e.g., 16)
    @assert length(umin_theta) == nu && length(umax_theta) == nu "Theta constraints must match nu=$nu"
    # @debug "Theta umin/umax first element: $(umin_theta[1]) / $(umax_theta[1])"
    if all(umin_theta .== 0.0) && all(umax_theta .== 0.0)
        @warn "Theta constraints umin and umax are both all zeros. Controller will likely output zeros."
    end


    # State constraints
    xmin = fill(-Inf, nx)
    xmax = fill(Inf, nx)
    # Optional: Load finite state bounds if needed from config
    # Example: xmin_config = get(mpc_params, "xmin", nothing) ... then fill xmin

    # --- NEW: Load physical constraints for u, v ---
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


    osqp_settings = get(mpc_params, "osqp_settings", Dict())

    @info "Setting up Alternating MPC Controller..."
    # MODIFIED: Pass physical bounds to setup_mpc_am
    mpc_data = setup_mpc_am(A, B, Np, Nc, Diagonal(Q_diag), Diagonal(R_diag_vec), xref,
                              umin_theta, umax_theta, xmin, xmax, m_decomp, n_decomp,
                              U_phys_min, U_phys_max, V_phys_min, V_phys_max; # NEW arguments
                              osqp_settings...)
    @info "MPC setup complete."

    sensor_indices = get(mpc_params, "sensor_indices", nothing)

    # Initialize controller state
    last_xk = nothing
    last_uk = nothing # Stores effective theta_k-1
    last_U_opt = nothing
    last_V_opt = nothing

    return AlternatingMPCController(mpc_data, sensor_indices, mean_vec, last_xk, last_uk, input_delays,
                                    N_iter_am, last_U_opt, last_V_opt) # Add new fields
end


# --- Implement the required interface function ---
function compute_control_action(controller::AlternatingMPCController, time::Float64, system_state::Dict{Symbol, Any})
    # --- Extract Data ---
    mpc_data = controller.mpc_data
    # MODIFIED: Extract new fields from mpc_data
    (; Sx, Su, Qbar, Rbar, P_sparse, A_con, A_in_theta, A_st_theta, Xref, Xmin, Xmax, l_in, u_in,
       U_phys_min, U_phys_max, V_phys_min, V_phys_max, # NEW
       num_theta_cons, num_u_phys_cons, num_v_phys_cons, num_total_cons, # NEW
       model_U, model_V, nx, nu, Nc, m, n, Np) = mpc_data
    N_iter_am = controller.N_iter_am
    mean_vec = controller.mean_vec
    du = controller.input_delays
    dy = nx - nu * du # Number of output delays

    # --- State Reconstruction ---
    # Get previous state x_{k-1} and previous *effective* input theta_{k-1}
    prev_xk = controller.last_xk === nothing ? zeros(nx) : controller.last_xk
    prev_theta = controller.last_uk === nothing ? zeros(nu) : controller.last_uk # nu = m*n
    #@debug "Time: $time, prev_xk[1]=$(round(prev_xk[1], digits=4)), prev_theta[1]=$(round(prev_theta[1], digits=4))"


    # Construct current state xk = [y_k, y_{k-1}, ..., theta_{k-1}, theta_{k-2}, ...]
    xk = zeros(nx)

    # 1. Get current measurement y_k (centered)
    local yk_centered::Float64
    measurement_key = :wallheatflux # Or load from config if variable
    if haskey(system_state, measurement_key)
        measurement_val_raw = system_state[measurement_key] ./ 1e6
        # Check if measurement_val_raw is empty or not a vector/scalar
        if !isa(measurement_val_raw, AbstractVector) && !isa(measurement_val_raw, Number)
             @error "Measurement '$measurement_key' is not a vector or number, type=$(typeof(measurement_val_raw)). Returning safe input."
             # Store current state (even if bad) and return zeros
             controller.last_xk = xk # Store partially constructed state
             controller.last_U_opt = zeros(m * Nc)
             controller.last_V_opt = zeros(n * Nc)
             controller.last_uk = zeros(nu)
             return zeros(m + n)
        end
        if isempty(measurement_val_raw)
             @error "Measurement '$measurement_key' is an empty vector. Returning safe input."
              # Store current state (even if bad) and return zeros
             controller.last_xk = xk # Store partially constructed state
             controller.last_U_opt = zeros(m * Nc)
             controller.last_V_opt = zeros(n * Nc)
             controller.last_uk = zeros(nu)
             return zeros(m + n)
        end
        # Take the first element if it's a vector
        measurement_val = isa(measurement_val_raw, AbstractVector) ? measurement_val_raw[1] : measurement_val_raw

        if length(mean_vec) >= 1
            yk_centered = Float64(measurement_val) - mean_vec[1]
            # @debug "Raw measurement: $(round(measurement_val, digits=4)), Centered y_k = $(round(yk_centered, digits=4))"
        else
            @error "mean_vec missing or too short, cannot center measurement. Returning safe input."
             # Store current state (even if bad) and return zeros
             controller.last_xk = xk # Store partially constructed state
             controller.last_U_opt = zeros(m * Nc)
             controller.last_V_opt = zeros(n * Nc)
             controller.last_uk = zeros(nu)
            return zeros(m + n) # Return safe physical input [u; v]
        end
    else
        @error "Missing measurement key '$measurement_key' in system_state. Returning safe input."
         # Store current state (even if bad) and return zeros
         controller.last_xk = xk # Store partially constructed state
         controller.last_U_opt = zeros(m * Nc)
         controller.last_V_opt = zeros(n * Nc)
         controller.last_uk = zeros(nu)
        return zeros(m + n)
    end
    xk[1] = yk_centered # y_k

    # 2. Past outputs: y_{k-1} down to y_{k-dy+1}
    if dy > 1
        if length(prev_xk) >= (dy - 1)
            xk[2:dy] = prev_xk[1:dy-1]
        else
            @error "prev_xk too short for output history reconstruction. length=$(length(prev_xk)), needed=$(dy-1). Returning safe input."
             # Store current state (even if bad) and return zeros
             controller.last_xk = xk # Store partially constructed state
             controller.last_U_opt = zeros(m * Nc)
             controller.last_V_opt = zeros(n * Nc)
             controller.last_uk = zeros(nu)
            return zeros(m + n)
        end
    elseif dy < 1 # Should not happen if nx, nu, du are correct, but check
         @warn "Number of output delays dy=$dy is not positive. State reconstruction might be incorrect."
    end

    # 3. Past inputs: theta_{k-1} down to theta_{k-du}
    if du > 0
        # theta_{k-1} is the prev_theta calculated in the last step
        start_idx_theta_km1 = dy + 1
        end_idx_theta_km1 = dy + nu
        if length(prev_theta) == nu && end_idx_theta_km1 <= nx && start_idx_theta_km1 > 0
            xk[start_idx_theta_km1 : end_idx_theta_km1] = prev_theta
        else
            @error "prev_theta length mismatch or index out of bounds for theta_{k-1}. length=$(length(prev_theta)), nu=$nu, start_idx=$start_idx_theta_km1, end_idx=$end_idx_theta_km1, nx=$nx. Returning safe input."
             # Store current state (even if bad) and return zeros
             controller.last_xk = xk # Store partially constructed state
             controller.last_U_opt = zeros(m * Nc)
             controller.last_V_opt = zeros(n * Nc)
             controller.last_uk = zeros(nu)
            return zeros(m + n)
        end

        # theta_{k-2} down to theta_{k-du}
        if du > 1
            # These come from prev_xk, starting after the outputs in prev_xk
            start_idx_prev_thetas = dy + 1 # Start after y delays in prev_xk
            num_older_thetas = (du - 1) * nu
            end_idx_prev_thetas = start_idx_prev_thetas + num_older_thetas - 1

            # Where they go in xk (after theta_{k-1})
            start_idx_new_thetas = end_idx_theta_km1 + 1
            end_idx_new_thetas = start_idx_new_thetas + num_older_thetas - 1

            # Check bounds before assignment
            if start_idx_prev_thetas > 0 && end_idx_prev_thetas <= length(prev_xk) &&
               start_idx_new_thetas > 0 && end_idx_new_thetas <= length(xk)
                 xk[start_idx_new_thetas : end_idx_new_thetas] = prev_xk[start_idx_prev_thetas : end_idx_prev_thetas]
            else
                 @error "Index out of bounds during state assembly for older thetas. Check dimensions/delays. Returning safe input."
                 @debug "Indices: prev=[$start_idx_prev_thetas:$end_idx_prev_thetas], new=[$start_idx_new_thetas:$end_idx_new_thetas], len(prev_xk)=$(length(prev_xk)), len(xk)=$(length(xk))"
                  # Store current state (even if bad) and return zeros
                  controller.last_xk = xk # Store partially constructed state
                  controller.last_U_opt = zeros(m * Nc)
                  controller.last_V_opt = zeros(n * Nc)
                  controller.last_uk = zeros(nu)
                 return zeros(m + n)
            end
        end
    end
    # --- End of State Construction ---

    # DEBUG: Check the reconstructed state
    @debug "Constructed state xk (norm, first 5): $(round(norm(xk), digits=4)), $(round.(xk[1:min(5, length(xk))], digits=4))"


    # --- Calculate State-Dependent QP Terms (for original Theta-QP constraints) ---
    state_error = Sx * xk - Xref
    q_theta_k = Su' * Qbar * state_error # Linear term for the Theta QP objective contribution

    # --- Calculate Original Theta constraints bounds (state-dependent part) ---
    l_st_k = Xmin - Sx * xk
    u_st_k = Xmax - Sx * xk

    # DEBUG: Check if state bounds are inverted
    @debug "State bounds check: any(l_st_k .> u_st_k) = $(any(l_st_k .> u_st_k))"
    if any(isfinite.(l_st_k)) && any(isfinite.(u_st_k)) # Avoid errors if all are Inf
        @debug "l_st_k finite min/max: $(minimum(filter(isfinite, l_st_k))) / $(maximum(filter(isfinite, l_st_k)))"
        @debug "u_st_k finite min/max: $(minimum(filter(isfinite, u_st_k))) / $(maximum(filter(isfinite, u_st_k)))"
    else
        @debug "l_st_k or u_st_k contains no finite values."
    end

    # Combine with fixed input bounds for theta
    l_theta_k = [l_in; l_st_k] # Combined lower bounds for A_con * Theta (size: num_theta_cons)
    u_theta_k = [u_in; u_st_k] # Combined upper bounds for A_con * Theta (size: num_theta_cons)


    #@debug "q_theta_k norm: $(round(norm(q_theta_k), digits=4))"
    #@debug "l_theta_k min/max: $(round(minimum(l_theta_k), digits=4)) / $(round(maximum(l_theta_k), digits=4))" # Less useful now, checked state part above
    #@debug "u_theta_k min/max: $(round(minimum(u_theta_k), digits=4)) / $(round(maximum(u_theta_k), digits=4))" # Less useful now, checked state part above


    # --- Initialize U and V for AM ---
    # Shift previous solution if available and non-zero, otherwise initialize with small constant non-zero values
    init_val = 0.01 # Small constant value for initialization
    warm_start_threshold = 1e-6 # Threshold to consider previous solution non-zero

    # --- Improved Initialization Logic ---
    if controller.last_U_opt !== nothing && length(controller.last_U_opt) == m * Nc && norm(controller.last_U_opt) > warm_start_threshold
        U_current = [controller.last_U_opt[m+1 : m*Nc]; zeros(m)] # Shift and pad
        #@debug "Initialized U_current from previous non-zero solution (shifted)."
    else
        #if controller.last_U_opt !== nothing
        #     @debug "Previous U_current was zero or invalid, re-initializing."
        #end
        U_current = fill(init_val, m * Nc) # Use fill for constant value
        #@debug "Initialized U_current with small constant values (val=$init_val)."
    end

    if controller.last_V_opt !== nothing && length(controller.last_V_opt) == n * Nc && norm(controller.last_V_opt) > warm_start_threshold
        V_current = [controller.last_V_opt[n+1 : n*Nc]; zeros(n)] # Shift and pad
         #@debug "Initialized V_current from previous non-zero solution (shifted)."
    else
         #if controller.last_V_opt !== nothing
         #    @debug "Previous V_current was zero or invalid, re-initializing."
         #end
        V_current = fill(init_val, n * Nc) # Use fill for constant value
         #@debug "Initialized V_current with small constant values (val=$init_val)."
    end
    # --- End Initialization Logic ---

    #@debug "Initial U_current norm: $(round(norm(U_current), digits=4))"
    #@debug "Initial V_current norm: $(round(norm(V_current), digits=4))"


    # --- Alternating Minimization Loop ---
    solve_success = true
    # MODIFIED: Increase regularization factor significantly
    qp_regularization = 1e1 # Increased from 1e-4, was originally 1e-6
    @info "Using QP Regularization: $qp_regularization" # Changed to info

    for iter = 1:N_iter_am
        #@debug "AM Iteration $iter/$N_iter_am"

        # --- 1. Optimize for U (given V_current) ---
        local results_U
        try
            M_V = build_M_V(V_current, Nc, m, n) # (nu*Nc) x (m*Nc)
            # Check norm AFTER building the matrix
            if norm(V_current) < 1e-9 || norm(M_V) < 1e-9 # Check if V_current or M_V is essentially zero
                 @warn "AM Iter $iter: V_current or M_V matrix is near zero. Skipping U optimization."
                 # U_current remains unchanged
            else
                # Calculate U-QP Hessian and Gradient
                P_U_dense = M_V' * P_sparse * M_V   # (m*Nc) x (m*Nc) - Can be dense
                P_U = sparse(0.5 * (P_U_dense + P_U_dense')) # Ensure symmetry
                P_U += (qp_regularization * I) # Add regularization

                q_U = M_V' * q_theta_k              # (m*Nc)

                # Map original theta constraints to U space
                A_con_U = A_con * M_V               # (num_theta_cons x m*Nc)

                # --- NEW: Build Augmented Constraints for U-Step ---
                # Identity matrix for U physical bounds (size: num_u_phys_cons x m*Nc)
                I_u_bounds = sparse(I, num_u_phys_cons, m * Nc)
                # Zero block placeholder for V bounds (size: num_v_phys_cons x m*Nc)
                Z_v_bounds_for_u_step = spzeros(num_v_phys_cons, m * Nc)

                # Stack vertically: [Theta_Cons_Mapped_to_U; U_Bounds_Identity; V_Bounds_Zero]
                A_U_aug = sparse([A_con_U; I_u_bounds; Z_v_bounds_for_u_step]) # Size: (num_total_cons, m*Nc)

                # Stack bounds vertically: [Theta_Bounds; U_Phys_Bounds; V_Bounds_Placeholder]
                l_U_aug = [l_theta_k; U_phys_min; fill(-Inf, num_v_phys_cons)] # Size: (num_total_cons,)
                u_U_aug = [u_theta_k; U_phys_max; fill(+Inf, num_v_phys_cons)] # Size: (num_total_cons,)
                # --- End Augmented Constraints ---


                # --- Minimal Debugging ---
                # DEBUG: Check final augmented bounds before solving
                @debug "Augmented U bounds check: any(l_U_aug .> u_U_aug) = $(any(l_U_aug .> u_U_aug))"
                 if any(isfinite.(l_U_aug)) && any(isfinite.(u_U_aug)) # Avoid errors if all are Inf
                    @debug "l_U_aug finite min/max: $(minimum(filter(isfinite, l_U_aug))) / $(maximum(filter(isfinite, l_U_aug)))"
                    @debug "u_U_aug finite min/max: $(minimum(filter(isfinite, u_U_aug))) / $(maximum(filter(isfinite, u_U_aug)))"
                else
                    @debug "l_U_aug or u_U_aug contains no finite values."
                end
                #@debug "     P_U size=$(size(P_U)), nnz=$(nnz(P_U))"
                #@debug "     q_U (used) size=$(size(q_U)), norm=$(round(norm(q_U), digits=4))"
                #@debug "     A_U_aug size=$(size(A_U_aug)), nnz=$(nnz(A_U_aug))"
                #@debug "     Initial U (warm start): size=$(size(U_current)), norm=$(round(norm(U_current), digits=4))"
                # --- End Minimal Debugging ---


                # Update and solve U-QP
                # MODIFIED: Use augmented constraints
                OSQP.update!(model_U; Px=P_U.nzval, q=q_U, Ax=A_U_aug.nzval, l=l_U_aug, u=u_U_aug)
                OSQP.warm_start!(model_U; x=U_current) # Warm start with previous U
                results_U = OSQP.solve!(model_U)
                #@debug "  U-Step OSQP Info: $(results_U.info)"


                if results_U.info.status_val in (1, 2) # Solved or solved inaccurately
                    U_current = results_U.x
                    #@debug "  U-Step solved. New U_current norm: $(round(norm(U_current), digits=4))"
                    if results_U.info.status_val == 2
                         @warn "AM Iter $iter: OSQP for U solved inaccurately."
                    end
                else
                    @warn "AM Iter $iter: OSQP for U failed with status code: $(results_U.info.status_val), status string: $(results_U.info.status). Keeping previous U."
                    # Log more details on failure
                    @debug "  U-Step Failed QP Data: P_U nnz=$(nnz(P_U)), q_U norm=$(norm(q_U)), A_U_aug nnz=$(nnz(A_U_aug)), l_U_aug range=($(minimum(l_U_aug)), $(maximum(l_U_aug))), u_U_aug range=($(minimum(u_U_aug)), $(maximum(u_U_aug)))"
                    solve_success = false
                    break # Exit AM loop on failure
                end
            end
        catch e
             @error "AM Iter $iter: Error during U optimization: $e" stacktrace=stacktrace(catch_backtrace())
             solve_success = false
             break # Exit AM loop on error
        end

        # --- 2. Optimize for V (given U_current) ---
        # Only proceed if U-step was successful and AM loop shouldn't break
        if !solve_success
            break
        end

        local results_V
        try
            M_U = build_M_U(U_current, Nc, m, n) # (nu*Nc) x (n*Nc)
             # Check norm AFTER building the matrix
             if norm(U_current) < 1e-9 || norm(M_U) < 1e-9 # Check if U_current or M_U is essentially zero
                 @warn "AM Iter $iter: U_current or M_U matrix is near zero. Skipping V optimization."
                 # V_current remains unchanged
             else
                # Calculate V-QP Hessian and Gradient
                P_V_dense = M_U' * P_sparse * M_U   # (n*Nc) x (n*Nc)
                P_V = sparse(0.5 * (P_V_dense + P_V_dense')) # Ensure symmetry
                P_V += (qp_regularization * I) # Add regularization

                q_V = M_U' * q_theta_k              # (n*Nc)

                # Map original theta constraints to V space
                A_con_V = A_con * M_U               # (num_theta_cons x n*Nc)

                # --- NEW: Build Augmented Constraints for V-Step ---
                # Identity matrix for V physical bounds (size: num_v_phys_cons x n*Nc)
                I_v_bounds = sparse(I, num_v_phys_cons, n * Nc)
                # Zero block placeholder for U bounds (size: num_u_phys_cons x n*Nc)
                Z_u_bounds_for_v_step = spzeros(num_u_phys_cons, n * Nc)

                # Stack vertically: [Theta_Cons_Mapped_to_V; U_Bounds_Zero; V_Bounds_Identity]
                A_V_aug = sparse([A_con_V; Z_u_bounds_for_v_step; I_v_bounds]) # Size: (num_total_cons, n*Nc)

                # Stack bounds vertically: [Theta_Bounds; U_Bounds_Placeholder; V_Phys_Bounds]
                l_V_aug = [l_theta_k; fill(-Inf, num_u_phys_cons); V_phys_min] # Size: (num_total_cons,)
                u_V_aug = [u_theta_k; fill(+Inf, num_u_phys_cons); V_phys_max] # Size: (num_total_cons,)
                # --- End Augmented Constraints ---


                # --- Minimal Debugging ---
                # DEBUG: Check final augmented bounds before solving (Optional for V-step unless it also fails)
                # @debug "Augmented V bounds check: any(l_V_aug .> u_V_aug) = $(any(l_V_aug .> u_V_aug))"
                #  if any(isfinite.(l_V_aug)) && any(isfinite.(u_V_aug)) # Avoid errors if all are Inf
                #     @debug "l_V_aug finite min/max: $(minimum(filter(isfinite, l_V_aug))) / $(maximum(filter(isfinite, l_V_aug)))"
                #     @debug "u_V_aug finite min/max: $(minimum(filter(isfinite, u_V_aug))) / $(maximum(filter(isfinite, u_V_aug)))"
                # else
                #     @debug "l_V_aug or u_V_aug contains no finite values."
                # end
                #@debug "     P_V size=$(size(P_V)), nnz=$(nnz(P_V))"
                #@debug "     q_V (used) size=$(size(q_V)), norm=$(round(norm(q_V), digits=4))"
                #@debug "     A_V_aug size=$(size(A_V_aug)), nnz=$(nnz(A_V_aug))"
                #@debug "     Initial V (warm start): size=$(size(V_current)), norm=$(round(norm(V_current), digits=4))"
                # --- End Minimal Debugging ---


                # Update and solve V-QP
                # MODIFIED: Use augmented constraints
                OSQP.update!(model_V; Px=P_V.nzval, q=q_V, Ax=A_V_aug.nzval, l=l_V_aug, u=u_V_aug)
                OSQP.warm_start!(model_V; x=V_current) # Warm start with previous V
                results_V = OSQP.solve!(model_V)
                #@debug "  V-Step OSQP Info: $(results_V.info)"


                if results_V.info.status_val in (1, 2) # Solved or solved inaccurately
                    V_current = results_V.x
                    #@debug "  V-Step solved. New V_current norm: $(round(norm(V_current), digits=4))"
                    if results_V.info.status_val == 2
                         @warn "AM Iter $iter: OSQP for V solved inaccurately."
                    end
                else
                    @warn "AM Iter $iter: OSQP for V failed with status code: $(results_V.info.status_val), status string: $(results_V.info.status). Keeping previous V."
                     # Log more details on failure
                    @debug "  V-Step Failed QP Data: P_V nnz=$(nnz(P_V)), q_V norm=$(norm(q_V)), A_V_aug nnz=$(nnz(A_V_aug)), l_V_aug range=($(minimum(l_V_aug)), $(maximum(l_V_aug))), u_V_aug range=($(minimum(u_V_aug)), $(maximum(u_V_aug)))"
                    solve_success = false
                    break # Exit AM loop on failure
                end
             end
        catch e
             @error "AM Iter $iter: Error during V optimization: $e" stacktrace=stacktrace(catch_backtrace())
             solve_success = false
             break # Exit AM loop on error
        end

        # Optional: Check for convergence?
        # if iter > 1 && norm(U_current - U_prev) < tol && norm(V_current - V_prev) < tol
        #     @debug "AM converged early at iter $iter"
        #     break
        # end
        # U_prev = copy(U_current) # Need to store previous values for convergence check
        # V_prev = copy(V_current)

    end # End AM loop

    # --- Extract Final Control Action ---
    if !solve_success
        @error "Alternating minimization failed to converge or encountered an error. Returning zero input."
        # Store zeros to avoid using potentially bad U_current/V_current
        controller.last_xk = xk
        controller.last_U_opt = zeros(m * Nc)
        controller.last_V_opt = zeros(n * Nc)
        controller.last_uk = zeros(nu) # Store zero effective theta
        return zeros(m + n)
    end

    uk_final = U_current[1:m]
    vk_final = V_current[1:n]
    #@debug "Final uk_final: $(round.(uk_final, digits=4))"
    #@debug "Final vk_final: $(round.(vk_final, digits=4))"


    # --- Store State and Input for Next Step ---
    controller.last_xk = xk
    controller.last_U_opt = U_current # Store full sequence for warm-start
    controller.last_V_opt = V_current # Store full sequence for warm-start

    # Store the *effective* theta that results from the applied u,v for state reconstruction
    theta_k_effective = vec(uk_final * vk_final')
    if length(theta_k_effective) != nu
        @error "Internal error: Effective theta dimension mismatch. Expected $nu, got $(length(theta_k_effective))."
        # Fallback if something went wrong with dimensions
         controller.last_uk = zeros(nu)
    else
         controller.last_uk = theta_k_effective
    end
    #@debug "Stored last_uk (effective theta) norm: $(round(norm(controller.last_uk), digits=4))"


    # --- Return Physical Input ---
    physical_input = [uk_final; vk_final]
    expected_physical_dim = m + n

    # Final checks (optional but recommended)
    if any(isnan, physical_input)
        @error "NaN detected in final control input. Returning zeros."
        # Reset stored optimal sequences if NaN occurred
        controller.last_U_opt = zeros(m * Nc)
        controller.last_V_opt = zeros(n * Nc)
        controller.last_uk = zeros(nu)
        return zeros(expected_physical_dim)
    end
     if length(physical_input) != expected_physical_dim
         @warn "Final physical input length incorrect ($(length(physical_input)) vs $expected_physical_dim). Returning zeros."
         return zeros(expected_physical_dim)
     end

    # Physical inputs should now be constrained by OSQP if bounds were finite.
    # Clamping here would only be needed if OSQP solution slightly violates bounds
    # due to tolerances, or if bounds were +/- Inf.
    # Example: clamp!(uk_final, u_phys_min_vec, u_phys_max_vec) etc. might still be useful.

    @info "Sending input: $(round.(physical_input, digits=4))" # Changed to info for visibility
    return physical_input
end

# --- update_controller_state! function ---
# Currently minimal, could be used for parameter adaptation later
function update_controller_state!(controller::AlternatingMPCController, time::Float64, new_state_data::Dict{Symbol, Any})
    # Example: Update reference if provided
    # if haskey(new_state_data, :target_reference)
    #     new_ref_val = new_state_data[:target_reference]
    #     # Update controller.mpc_data.xref and controller.mpc_data.Xref accordingly
    #     # ... (implementation needed) ...
    #     @info "Updated MPC reference at time $time to $new_ref_val"
    # end
    return nothing
end
