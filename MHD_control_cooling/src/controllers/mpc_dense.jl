# src/controllers/mpc_dense.jl
using LinearAlgebra
using OSQP
using SparseArrays
using JLD2 # If loading A, B from file
using TOML # To read config inside constructor
using Logging # Added for @debug, @warn, @error

# Assuming AbstractController is defined in interface.jl
include("interface.jl")
# Include the outer product optimization function
include("../optim/outer_product_optim.jl") # Make sure path is correct

# --- MPCData Struct Definition ---
# (No changes needed here)
struct MPCData
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
# (No changes needed here)
function setup_mpc(A, B, Np, Nc, Q, R, xref, umin, umax, xmin, xmax; osqp_settings...)
    nx = size(A, 1)
    nu = size(B, 2) # Should be 16
    @info "Setting up MPC with nx=$nx, nu=$nu"
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
# (No changes needed here)
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
    elseif results.info.status_val == 2
         @warn "OSQP solved inaccurately. Status: $(results.info.status)"
    end
    # Add solver logging

    U_opt = results.x

    if maximum(U_opt) > 1 || minimum(U_opt) < 0
        @warn "OSQP solution out of bounds. Maximum: $(maximum(U_opt)), Minimum: $(minimum(U_opt))"
        clamp!(U_opt, 0, 1)
    end

    if isnothing(U_opt)
        @error "OSQP solve returned nothing despite status being ok? Returning zeros."
        return zeros(nu)
    end
    return U_opt[1:nu] # Return first control input vector (16 elements)
end

# --- The Controller Struct ---
# (No changes needed here)
mutable struct DenseMPCController <: AbstractController
    mpc_data::MPCData
    sensor_indices::Union{Vector{Int}, Nothing} # Optional: Indices to map system_state to parts of xk
    mean_vec::Vector{Float64}
    last_xk::Union{Vector{Float64}, Nothing}    # Stores x_{k-1}
    last_uk::Union{Vector{Float64}, Nothing}    # Stores theta_{k-1}
    input_delays::Int                           # **** ADDED: Store number of input delays ****
end

# --- Constructor for the Controller ---
function DenseMPCController(config::Dict)
    mpc_params = config["controller"]["params"]["DenseMPC"]
    last_xk = nothing
    last_uk = nothing

    matrix_file = mpc_params["matrix_file"]
    local mean_vec, A, B
    try
        jldopen(matrix_file, "r") do file
            A = file["A"]
            B = file["B"] # B should be nx_aug x 16
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
    nu = size(B, 2) # Should be 16

    if nu != 16
        @warn "Loaded B matrix has nu=$(nu) columns, but expected 16 based on outer product plan. Proceeding, but check model/config."
    end

    # **** CRITICAL: Ensure this parameter exists in your config file! ****
    input_delays = get(mpc_params, "input_delays", -1) # Default to -1 to force error if missing
    if input_delays < 0
         @error "Missing or invalid 'input_delays' parameter in controller config. This is required for state reconstruction."
         throw(ArgumentError("Missing 'input_delays' in config"))
    elseif input_delays * nu >= nx
         @error "Invalid delay structure: input_delays*nu ($(input_delays*nu)) >= nx ($nx). Check config/model."
         throw(ArgumentError("Invalid delay structure"))
    end
    @info "Using input_delays (du) = $input_delays"
    # **** End of added section ****

    # Construct Q ensuring diagonal matches nx
    Q_base = Float64(get(mpc_params, "Q", 1.0)) # Get base Q value, default 1.0
    Q_diag_override = get(mpc_params, "Q_diag", nothing) # Optional override vector
    local Q_diag
    if isnothing(Q_diag_override)
        Q_diag = ones(nx) * Q_base
        # Optionally emphasize the first state (output) if desired
        # Q_diag[1] *= 10 # Example: Make first state 10x more important
    else
        Q_diag = Float64.(Q_diag_override)
    end
    @assert length(Q_diag) == nx "Length of Q_diag must match state dimension nx=$nx"

    # Construct R from R_diag
    R_diag_vec = Float64.(mpc_params["R_diag"]) # Length should be 16
    @assert length(R_diag_vec) == nu "Length of R_diag must match input dimension nu=$nu (expected 16)"

    # Construct xref - ASSUMING STEADY STATE where output = x0ref and inputs are zero
    # TODO: This might need refinement if the steady-state theta is non-zero for x0ref
    x0ref_val = Float64(mpc_params["x0ref"])
    xref = zeros(nx)
    dy = nx - input_delays * nu # Calculate number of output delays
    xref[1:dy] .= x0ref_val # Set all output delay states to the target reference
    # Keep input delay states at 0 (assuming zero steady-state input for target)
    @info "Constructed xref with first element: $(xref[1])"
    # @assert length(xref) == nx "Length of xref must match state dimension nx=$nx" # Already checked implicitly

    umin = Float64.(mpc_params["umin"]) # Length should be 16
    umax = Float64.(mpc_params["umax"]) # Length should be 16

    xmin = fill(-Inf, nx) # Use fill for clarity
    xmax = fill(Inf, nx)
    @assert length(umin) == nu && length(umax) == nu "Input constraints must match nu=$nu (expected 16)"
    @assert length(xmin) == nx && length(xmax) == nx "State constraints must match nx=$nx"

    osqp_settings = get(mpc_params, "osqp_settings", Dict())

    @info "Setting up Dense MPC Controller (nu=16)..."
    mpc_data = setup_mpc(A, B, Np, Nc, Diagonal(Q_diag), Diagonal(R_diag_vec), xref, umin, umax, xmin, xmax; osqp_settings...)
    @info "MPC setup complete."

    sensor_indices = get(mpc_params, "sensor_indices", nothing)

    # Pass input_delays to the struct
    return DenseMPCController(mpc_data, sensor_indices, mean_vec, last_xk, last_uk, input_delays)
end


# --- Implement the required interface function ---
# --- MODIFIED: Uses Direct State Assembly ---
function compute_control_action(controller::DenseMPCController, time::Float64, system_state::Dict{Symbol, Any})
    # --- Extract Data ---
    mpc_data = controller.mpc_data
    nx = mpc_data.nx
    nu = mpc_data.nu # 16
    mean_vec = controller.mean_vec
    du = controller.input_delays # Get number of input delays from struct
    dy = nx - nu * du # Number of output delays

    # --- Get Previous State (x_{k-1}) and Input (theta_{k-1}) ---
    prev_xk = controller.last_xk === nothing ? zeros(nx) : controller.last_xk
    prev_theta = controller.last_uk === nothing ? zeros(nu) : controller.last_uk

    # --- CONSTRUCT CURRENT STATE xk = [y_k, y_{k-1}, ..., theta_{k-1}, ...] ---
    xk = zeros(nx)

    # 1. Get current measurement y_k (centered)
    local yk_centered::Float64
    measurement_key = :wallheatflux # Or load from config if variable
    if haskey(system_state, measurement_key)
        # Assuming the measurement is a scalar or we take the first element
        measurement_val = isa(system_state[measurement_key], AbstractVector) ?
                          system_state[measurement_key][1] :
                          system_state[measurement_key]
        measurement_val

        if length(mean_vec) >= 1
             yk_centered = Float64(measurement_val) - mean_vec[1]
             @info "Centered measurement y_k = $yk_centered"
        else
             @error "mean_vec missing or too short, cannot center measurement."
             # Fallback: Use raw measurement? Risky. Return safe input.
             # yk_centered = Float64(measurement_val)
             return zeros(8) # Return safe physical input
        end
    else
        @error "Missing measurement key '$measurement_key' in system_state."
        return zeros(8) # Return safe physical input
    end
    yk_centered = yk_centered / 1e6
    xk[1] = yk_centered # y_k

    # 2. Past outputs: y_{k-1} down to y_{k-dy+1}
    if dy > 1
        # These come from the first dy-1 elements of prev_xk
        if length(prev_xk) >= (dy - 1)
             xk[2:dy] = prev_xk[1:dy-1]
        else
             @error "prev_xk is too short for output history reconstruction. length=$(length(prev_xk)), needed=$(dy-1)"
             return zeros(8)
        end
    end

    # 3. Past inputs: theta_{k-1} down to theta_{k-du}
    if du > 0
        # theta_{k-1} is the prev_theta calculated in the last step
        start_idx_theta_km1 = dy + 1
        end_idx_theta_km1 = dy + nu
        if length(prev_theta) == nu && end_idx_theta_km1 <= nx
            xk[start_idx_theta_km1 : end_idx_theta_km1] = prev_theta
        else
             @error "prev_theta length mismatch or index out of bounds for theta_{k-1}. length=$(length(prev_theta)), nu=$nu, end_idx=$end_idx_theta_km1, nx=$nx"
             return zeros(8)
        end


        # theta_{k-2} down to theta_{k-du}
        if du > 1
            # These come from prev_xk, starting after the outputs in prev_xk
            start_idx_prev_thetas = dy + 1
            num_older_thetas = (du - 1) * nu
            end_idx_prev_thetas = start_idx_prev_thetas + num_older_thetas - 1

            # Where they go in xk (after theta_{k-1})
            start_idx_new_thetas = end_idx_theta_km1 + 1
            end_idx_new_thetas = start_idx_new_thetas + num_older_thetas - 1

            # Check bounds before assignment
            if end_idx_prev_thetas <= length(prev_xk) && end_idx_new_thetas <= length(xk) && start_idx_prev_thetas > 0 && start_idx_new_thetas > 0
                 xk[start_idx_new_thetas : end_idx_new_thetas] = prev_xk[start_idx_prev_thetas : end_idx_prev_thetas]
            else
                 @error "Index out of bounds during state assembly for older thetas. Check dimensions/delays."
                 @debug "Indices: prev=[$start_idx_prev_thetas:$end_idx_prev_thetas], new=[$start_idx_new_thetas:$end_idx_new_thetas], len(prev_xk)=$(length(prev_xk)), len(xk)=$(length(xk))"
                 return zeros(8)
            end
        end
    end
    # --- End of State Construction ---

    @debug "Constructed state xk[1]=$(round(xk[1], digits=3)). Using prev_theta[1]=$(round(prev_theta[1], digits=3))"

    # --- Compute MPC Control Input (theta_k) ---
    current_theta_k = compute_control_input(xk, mpc_data) # Result is theta_k

    # --- Store State and Input for Next Step ---
    controller.last_xk = xk             # Store assembled x_k
    controller.last_uk = current_theta_k # Store calculated theta_k

    # --- Decompose theta_k to Physical Input [u_k; v_k] ---
    local physical_input_8::Vector{Float64}
    m_decomp, n_decomp = 4, 4 # Physical input dimensions
    expected_physical_dim = m_decomp + n_decomp # Should be 8

    if length(current_theta_k) == nu # Ensure correct length (16)
        try
            uk, vk = find_closest_outer_product(current_theta_k, m_decomp, n_decomp)
            physical_input_8 = [uk; vk]
            reconstruction_error = norm(current_theta_k - vec(uk * vk'))
            @debug "Decomposed theta_k to physical [uk; vk]: $(round.(physical_input_8, digits=3)). Recon Error: $(round(reconstruction_error, digits=4))"

        catch e
            @error "Outer product decomposition failed: $e. Using zeros."
            physical_input_8 = zeros(expected_physical_dim)
        end
    else
        @warn "MPC returned vector length $(length(current_theta_k)), expected $nu. Using zeros."
        physical_input_8 = zeros(expected_physical_dim)
    end

     # Final dimension check
     if length(physical_input_8) != expected_physical_dim
         @warn "Final physical input length incorrect ($(length(physical_input_8)) vs $expected_physical_dim). Using zeros."
         physical_input_8 = zeros(expected_physical_dim)
     end

    return physical_input_8
end

# --- update_controller_state! function ---
# (No changes needed here)
function update_controller_state!(controller::DenseMPCController, time::Float64, new_state_data::Dict{Symbol, Any})
    # This function is currently not used for state updates with direct assembly,
    # but could be used for parameter updates or other logic if needed.
    return nothing
end
