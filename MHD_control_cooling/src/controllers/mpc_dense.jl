# src/controllers/mpc_dense.jl
using LinearAlgebra
using OSQP
using SparseArrays
using JLD2 # If loading A, B from file
using TOML # To read config inside constructor

# Assuming AbstractController is defined in interface.jl
include("interface.jl")
# Include the outer product optimization function
include("../optim/outer_product_optim.jl") # Make sure path is correct

# --- MPCData Struct Definition ---
struct MPCData
    # ... (definition remains the same, but nu will be 16) ...
    A::Matrix{Float64}
    B::Matrix{Float64} # B will be nx_aug x 16
    nx::Int
    nu::Int # This will be 16
    Np::Int
    Nc::Int
    Q::Diagonal{Float64, Vector{Float64}}
    R::Diagonal{Float64, Vector{Float64}} # R will be 16x16
    xref::Vector{Float64}
    umin::Vector{Float64} # umin will have length 16
    umax::Vector{Float64} # umax will have length 16
    xmin::Vector{Float64}
    xmax::Vector{Float64}
    Sx::Matrix{Float64}
    Su::Matrix{Float64} # Su will be (nx*Np) x (16*Nc)
    Qbar::SparseMatrixCSC{Float64, Int64}
    Rbar::SparseMatrixCSC{Float64, Int64}
    P_sparse::SparseMatrixCSC{Float64, Int64}
    A_in::SparseMatrixCSC{Float64, Int64}
    A_st::SparseMatrixCSC{Float64, Int64}
    A_con::SparseMatrixCSC{Float64, Int64}
    l_in::Vector{Float64} # length 16*Nc
    u_in::Vector{Float64} # length 16*Nc
    Xmin::Vector{Float64}
    Xmax::Vector{Float64}
    Xref::Vector{Float64}
    model::OSQP.Model
end

# --- setup_mpc function ---
# (Implementation is the same, but it will work with nu=16 based on the input B matrix)
function setup_mpc(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax; osqp_settings...)
    nx = size(A, 1)
    nu = size(B, 2) # Should be 16
    @info "Setting up MPC with nx=$nx, nu=$nu"
    # ... (rest of setup_mpc implementation as before) ...
     # --- Dimension Checks ---
     @assert size(A, 1) == size(A, 2) "Matrix A must be square ($nx x $nx)."
     @assert size(B, 1) == nx "Matrix B must have $nx rows (dimension mismatch with A)."
     @assert Np > 0 "Prediction horizon Np must be positive."
     @assert Nc > 0 "Control horizon Nc must be positive."
     @assert size(Q) == (nx, nx) "Matrix Q must be square with dimension $nx x $nx."
     @assert size(R) == (nu, nu) "Matrix R must be square with dimension $nu x $nu." # nu=16 check
     @assert length(xref) == nx "Reference vector xref must have length $nx."
     @assert length(umin) == nu "Control lower bound vector umin must have length $nu." # nu=16 check
     @assert length(umax) == nu "Control upper bound vector umax must have length $nu." # nu=16 check
     @assert length(xmin) == nx "State lower bound vector xmin must have length $nx."
     @assert length(xmax) == nx "State upper bound vector xmax must have length $nx."

     # --- Dense MPC Formulation Calculations ---
     Sx = zeros(nx * Np, nx)
     Su = zeros(nx * Np, nu * Nc) # nu=16 here

     # Calculate Sx
     temp_A = A
     for i = 1:Np
         rows = (i-1)*nx+1 : i*nx
         Sx[rows, :] = temp_A
         temp_A = temp_A * A
     end

     # Calculate Su using the input Matrix B directly (B is nx x 16)
     for i = 1:Np
         rows_x = (i-1)*nx+1 : i*nx
         for j = 1:min(i, Nc)
             cols_u = (j-1)*nu+1 : j*nu # nu=16 here
             if i - j >= 0
                 Su[rows_x, cols_u] = (A^(i - j)) * B
             end
         end
     end

     # Build QP matrices
     Qbar = kron(sparse(I(Np)), sparse(Q))
     Rbar = kron(sparse(I(Nc)), sparse(R)) # R is 16x16

     Hessian_calc = Su' * Qbar * Su + Rbar
     P_sparse = sparse(Symmetric(Hessian_calc))

     # --- Constraint Formulation ---
     Umin = repeat(umin, Nc) # umin has length 16
     Umax = repeat(umax, Nc) # umax has length 16
     A_in = sparse(Matrix{Float64}(I, nu * Nc, nu * Nc)) # nu=16 here
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
     default_settings = Dict(:verbose => false, :eps_abs => 1e-4, :eps_rel => 1e-4, :max_iter => 5000)
     settings = merge(default_settings, Dict(osqp_settings))

     # Note: q vector length is nu*Nc (16*Nc)
     OSQP.setup!(model; P=P_sparse, q=zeros(nu*Nc), A=A_con, l=zeros(size(A_con,1)), u=zeros(size(A_con,1)), settings...)

     # Create and return the MPCData struct (passing the original Matrix B)
     return MPCData(A, B, nx, nu, Np, Nc, Q, R, xref, umin, umax, xmin, xmax,
                    Sx, Su, Qbar, Rbar, P_sparse, A_in, A_st, A_con,
                    l_in, u_in, Xmin, Xmax, Xref, model)

end

# --- compute_control_input function ---
# (Implementation is the same, but it returns a 16-element vector uk_vec)
function compute_control_input(xk, mpc_data::MPCData)
    (; Sx, Su, Qbar, Xref, Xmin, Xmax, l_in, u_in, A_con, model, nu, nx) = mpc_data # nu = 16 here

    # Update QP terms
    q_k = Su' * Qbar * (Sx * xk - Xref)
    l_st_k = Xmin - Sx * xk
    u_st_k = Xmax - Sx * xk
    l_k = [l_in; l_st_k]
    u_k = [u_in; u_st_k]

    # Update and solve
    OSQP.update!(model; q=q_k, l=l_k, u=u_k)
    results = OSQP.solve!(model)

    # Check status and return result
    if results.info.status_val ∉ (1, 2) # 1: solved, 2: solved inaccurate
         @warn "OSQP solver failed with status: $(results.info.status)"
         # Return zero control (16 elements) or handle appropriately
         return zeros(nu) # nu = 16
     end
    U_opt = results.x
    return U_opt[1:nu] # Return first control input vector (16 elements)
end

# --- The Controller Struct ---
mutable struct DenseMPCController <: AbstractController
    mpc_data::MPCData
    sensor_indices::Union{Vector{Int}, Nothing}
    mean_vec::Vector{Float64}
    last_xk::Union{Vector{Float64}, Nothing}
    last_uk::Union{Vector{Float64}, Nothing}
end

# --- Constructor for the Controller ---
function DenseMPCController(config::Dict)
    mpc_params = config["controller"]["params"]["DenseMPC"] # Expect DenseMPC section in config
    last_xk = nothing
    last_uk = nothing

    matrix_file = mpc_params["matrix_file"]
    local mean_vec
    try
        jldopen(matrix_file, "r") do file
            global A = file["A"]
            global B = file["B"] # B should be nx_aug x 16
            mean_vec = file["mean_vec"]
            @info "Loaded A [$(size(A))] and B [$(size(B))] from $matrix_file"
        end
    catch e
        @error "Failed to load matrices from $matrix_file" exception=(e, catch_backtrace())
        rethrow(e)
    end

    # Extract other MPC parameters from config (checking nu=16 constraints)
    Np = mpc_params["Np"]
    Nc = mpc_params["Nc"]
    nx = size(A, 1) # Augmented state dimension
    nu = size(B, 2) # Should be 16

    if nu != 16
        @warn "Loaded B matrix has nu=$(nu) columns, but expected 16 based on outer product plan. Proceeding, but check model/config."
    end

    Q_diag = ones(nx) * Float64.(mpc_params["Q"])
    R_diag = Float64.(mpc_params["R_diag"]) # Length should be 16
    @assert length(Q_diag) == nx "Length of Q_diag must match state dimension nx=$nx"
    @assert length(R_diag) == nu "Length of R_diag must match input dimension nu=$nu (expected 16)"

    xref = Float64.(mpc_params["x0ref"]) * ones(nx)
    @assert length(xref) == nx "Length of xref must match state dimension nx=$nx"

    umin = Float64.(mpc_params["umin"]) # Length should be 16
    umax = Float64.(mpc_params["umax"]) # Length should be 16
    xmin = -Inf * ones(nx)
    xmax = Inf * ones(nx)
    @assert length(umin) == nu && length(umax) == nu "Input constraints must match nu=$nu (expected 16)"
    @assert length(xmin) == nx && length(xmax) == nx "State constraints must match nx=$nx"

    osqp_settings = get(mpc_params, "osqp_settings", Dict())

    @info "Setting up Dense MPC Controller (nu=16)..."
    mpc_data = setup_mpc(A, B, Np, Nc, Diagonal(Q_diag), Diagonal(R_diag), xref, umin, umax, xmin, xmax; osqp_settings...)
    @info "MPC setup complete."

    sensor_indices = get(mpc_params, "sensor_indices", nothing)

    return DenseMPCController(mpc_data, sensor_indices, mean_vec, last_xk, last_uk)
end

# --- Implement the required interface function ---
function compute_control_action(controller::DenseMPCController, time::Float64, system_state::Dict{Symbol, Any})
    # --- STATE MAPPING ---
    # (Logic remains the same - needs to be configured correctly via sensor_indices)
    local last_uk, last_xk
    if controller.last_xk !== nothing && controller.last_uk !== nothing
        last_xk = controller.last_xk
        last_uk = controller.last_uk
    else
        last_xk = zeros(controller.mpc_data.nx)
        last_uk = zeros(controller.mpc_data.nu)
    end

    local xk::Vector{Float64}
    xk = A * last_xk + B * last_uk
    xk[1] = system_state[:wallheatflux][1] - controller.mean_vec[1]

    @info "Target reference is $(xk[1])"

    # --- Call the core MPC computation (returns 16 elements) ---
    uk_vec_16 = compute_control_input(xk, controller.mpc_data)
    controller.last_uk = uk_vec_16

    # --- Decompose 16-element vector into 8-element physical input ---
    local physical_input_8::Vector{Float64}
    if length(uk_vec_16) == 16
        try
            # Decompose uk_vec_16 (theta) into u (4x1) and v (4x1)
            # Assumes m=4, n=4 based on dmd_lib.jl preprocessing
            u_opt, v_opt = find_closest_outer_product(uk_vec_16, 4, 4)

            # Concatenate u and v to get the 8-element physical input
            physical_input_8 = [u_opt; v_opt]

        catch e
            @error "Error during outer product decomposition: $e. Falling back to zero input."
            physical_input_8 = zeros(8)
        end
    else
         @warn "MPC returned vector of length $(length(uk_vec_16)), expected 16. Cannot decompose. Falling back to zero input."
         physical_input_8 = zeros(8)
    end

    # Ensure final vector has 8 elements
    if length(physical_input_8) != 8
        @warn "Final physical input has length $(length(physical_input_8)), expected 8. Using zeros."
        physical_input_8 = zeros(8)
    end

    return physical_input_8 # Return the 8-element vector for the physical system
end

# --- update_controller_state! function remains the same (likely does nothing) ---
function update_controller_state!(controller::DenseMPCController, time::Float64, new_state_data::Dict{Symbol, Any})
    # ... (no changes needed unless adding an estimator) ...
    return nothing
end
