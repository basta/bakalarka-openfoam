# src/controllers/dmdc_mpc.jl
# Implements a Controller using DMDc for system identification and Dense MPC for control.

# --- Dependencies ---
using LinearAlgebra
using CircularArrayBuffers # For efficient history storage (Pkg.add("CircularArrayBuffers"))
using SparseArrays
using JuMP
# Assuming these libraries are in the same project/environment
try
    # Controller Interface
    include("interface.jl") # Defines AbstractController etc.

    # DMDc library for system ID (or loading precomputed A/B)
    include("../dmd_lib.jl") # Provides create_abc_system_from_data, create_delay_state etc.

    # Dense MPC library
    include("../mpc/mpc_lib.jl")
    using .DenseMPCLib # Provides DenseMPCController struct, setup_dense_mpc, compute_control_action

    include("../optim/outer_product_optim.jl")
catch e
    @error "Failed to include required library files (interface.jl, dmd_lib.jl, mpc_lib.jl)." exception=(e, catch_backtrace())
    rethrow(e)
end


# --- Placeholder Functions (Implement these based on your specific setup) ---

"""
Estimates the N_X-dimensional physical state vector 'x' from the raw server data.
"""
function estimate_physical_state(state_dict::Dict{Symbol, Any}, N_X::Int)::Vector{Float64}    # end
    if haskey(state_dict, :wallheatflux) && !isempty(state_dict[:wallheatflux])
        return state_dict[:wallheatflux][1:1]
    else
        return zeros(Float64, N_X) # Fallback
    end
end

"""
Preprocesses the m_phys-dimensional physical input 'u_phys' into the
N_U-dimensional input 'u_prime' used by the DMDc model's B matrix.
"""
function preprocess_input(u_phys::AbstractVector{Float64}, N_U::Int)::Vector{Float64}
    # TODO: Implement the input preprocessing used during DMDc identification.
    # Example: If u' = vec(u_e * u_m'), where u_phys = [u_e; u_m]
    m_phys = length(u_phys)
    # @info u_phys, "requested N_U is $N_U"
    if m_phys % 2 != 0 || N_U != (m_phys/2)^2
         @warn "Input preprocessing assumes u_phys = [u_e; u_m] and N_U = (m_phys/2)^2. Mismatch detected. Returning zeros."
         error("Input dim mismatch")
         return zeros(Float64, N_U)
    end
    split_idx = m_phys ÷ 2
    u_e = u_phys[1:split_idx]
    u_m = u_phys[split_idx+1:end]
    u_prime = vec(u_e * u_m') # Perform outer product and flatten
    if length(u_prime) != N_U
         @error "Preprocessed input dimension ($(length(u_prime))) does not match N_U ($N_U)."
         return zeros(Float64, N_U) # Fallback
    end
    return u_prime
end


# --- Controller Definition ---

"""
    DMDC_MPCController <: AbstractController

Controller implementing Dense MPC based on a DMDc-identified model with delays.

Handles state augmentation, normalization, and calls the underlying DenseMPCLib.
"""
mutable struct DMDC_MPCController <: AbstractController
    # --- Configuration & Parameters ---
    Np::Int                     # Prediction horizon
    N_X::Int                    # Dimension of physical state 'x'
    N_U::Int                    # Dimension of preprocessed input 'u_prime'
    m_phys::Int                 # Dimension of physical input 'u_phys'
    DELAY::Int                  # Number of delay steps in the model

    # --- Core Components ---
    dense_controller::DenseMPCLib.DenseMPCController # Pre-computed dense MPC solver setup
    mean_vec::Vector{Float64}   # Mean vector for physical state normalization (size N_X)

    # --- History Buffers ---
    # Store recent physical states and inputs to construct augmented state
    # CircularArrayBuffer is efficient for rolling windows
    x_history::CircularArrayBuffer{Float64, 2} # Stores x_k, x_{k-1}, ..., x_{k-DELAY} (size N_X x DELAY+1)
    u_phys_history::CircularArrayBuffer{Float64, 2} # Stores u_phys_{k-1}, ..., u_phys_{k-DELAY} (size m_phys x DELAY)

    # --- Last Computed Input ---
    last_computed_u_phys::Vector{Float64} # Store u_phys_k for next step's history

    # --- Constructor (Performs Offline Setup) ---
    function DMDC_MPCController(config::Dict)
        @info "Initializing DMDC_MPCController..."

        # --- Extract Parameters from Config ---
        mpc_params = get(config, "controller", Dict())
        mpc_specific_params = get(get(mpc_params, "params", Dict()), "DMDC_MPC", Dict()) # Get params specific to this controller type

        Np = get(mpc_specific_params, "Np", 5) # Prediction Horizon
        m_phys = get(mpc_specific_params, "m_phys", 8) # Physical input dimension
        data_path = get(mpc_specific_params, "data_path", "./data/dataset.jld2") # Path for DMDc ID
        abc_path = get(mpc_specific_params, "abc_path", "ABC.jld2") # Path to load/save A/B/mean/delay
        load_abc = get(mpc_specific_params, "load_abc", true) # Flag to load precomputed A/B

        # --- System Identification (DMDc) ---
        local A_full, B_full, mean_vec, DELAY, N_X, N_U
        if load_abc && isfile(abc_path)
            try
                @info "Loading pre-computed A, B, mean_vec, delay from '$abc_path'..."
                jldopen(abc_path, "r") do file
                    A_full = file["A"]
                    B_full = file["B"]
                    mean_vec = file["mean_vec"]
                    DELAY = file["delay"]
                end
                # Infer dimensions (ensure consistency)
                N_X = length(mean_vec)
                N_aug = size(A_full, 1)
                if DELAY == 0
                    N_U = size(B_full, 2)
                    if N_aug != N_X*(DELAY+1)
                         error("Loaded A_full dimension ($(size(A_full,1))) inconsistent with N_X=$N_X, DELAY=$DELAY")
                    end
                else
                    # Infer N_U assuming structure z = [x_part; u_part]
                    expected_x_part_len = N_X * (DELAY + 1)
                    u_part_len = N_aug - expected_x_part_len
                    if u_part_len < 0 || u_part_len % DELAY != 0
                        error("Cannot infer N_U from loaded A_full dimension ($N_aug), N_X ($N_X), and DELAY ($DELAY).")
                    end
                    N_U = u_part_len ÷ DELAY
                     # Check B_full columns match N_U
                    if size(B_full, 2) != N_U
                         error("Loaded B_full column count ($(size(B_full, 2))) does not match inferred N_U ($N_U).")
                    end
                end
                 @info "Loaded System: N_X=$N_X, N_U=$N_U, DELAY=$DELAY, N_aug=$N_aug"
            catch e
                @error "Failed to load ABC data from '$abc_path'. Performing DMDc identification..." exception=(e, catch_backtrace())
                load_abc = false # Force recomputation
            end
        end

        if !load_abc || !@isdefined(A_full) # If loading failed or was skipped
            @info "Performing DMDc identification from data path '$data_path'..."
            # DELAY needs to be defined before calling create_abc_system_from_data
            # Get DELAY from config or use a default
            DELAY = get(mpc_specific_params, "delay", 10) # Default delay if not loading/specified
            if DELAY < 0 error("DELAY must be non-negative.") end

            # This function needs the raw X_data and U_data (preprocessed)
            # It returns the augmented A, B, the mean of X, N_X, N_U (preprocessed dim), etc.
            A_full_id, B_full_id, _, _, _, mean_vec_id, N_X_id, N_U_id, _, _, _ = create_abc_system_from_data(data_path, DELAY)
            A_full = A_full_id
            B_full = B_full_id
            mean_vec = mean_vec_id
            N_X = N_X_id
            N_U = N_U_id # This is the dimension of the *preprocessed* U used in B_full
            @info "DMDc Identification Complete: N_X=$N_X, N_U=$N_U, DELAY=$DELAY"
            # Optionally save the computed matrices
            try
                jldsave(abc_path, A=A_full, B=B_full, delay=DELAY, mean_vec=mean_vec)
                @info "Saved computed A, B, mean_vec, delay to '$abc_path'"
            catch e
                 @warn "Could not save computed ABC data to '$abc_path'." exception=e
            end
        end

        # --- Define MPC Costs & Constraints ---
        # TODO: Load these from the config file (mpc_specific_params)
        # Example: Define Q, R, E, F vectors based on config values
        @warn "Using placeholder MPC costs and constraints. Load/define actual values from config."
        T = eltype(A_full) # Get data type
        N_aug = size(A_full, 1) # Dimension of augmented state z
        m_aug = size(B_full, 2) # Dimension of augmented input u' (should match N_U)

        Q_state_val = 1.0; R_phys_val = 0.05; Q_term_val = 10.0; target_state_val = 0
        Q_stage = spzeros(T, N_aug, N_aug); # Penalty on augmented state
        # Penalize deviation of the *physical* state part (first N_X elements) from target
        target_phys_state_deviation = zeros(N_X); # Assume target is zero mean for simplicity here
        # Example: Penalize first state element
        Q_stage[1, 1] = Q_state_val
        Q_term_mat = spzeros(T, N_aug, N_aug); Q_term_mat[1, 1] = Q_term_val

        # R penalizes the *augmented* input u' (size m_aug = N_U)
        # If you want to penalize physical input u_phys, it's more complex in this formulation
        R_stage = sparse(Diagonal(fill(R_phys_val, m_aug))) # Penalize the preprocessed input

        r_phys = -1e6*ones(N_X)
        r_prime_phys = r_phys - mean_vec # Target for normalized x_k
        z_prime_ref = zeros(T, N_aug)
        z_prime_ref[1:N_X] = r_prime_phys

        # Linear terms (q, r) - usually zero unless tracking non-zero reference
        q_stage = -2 * Q_stage * z_prime_ref
        q_term_vec = -2 * Q_term_mat * z_prime_ref
        r_stage = zeros(T, m_aug)

        Q_vec = AbstractMatrix{T}[copy(Q_stage) for _ = 1:Np]; push!(Q_vec, Q_term_mat);
        R_vec = AbstractMatrix{T}[copy(R_stage) for _ = 1:Np]
        q_vec = AbstractVector{T}[copy(q_stage) for _ = 1:Np]; push!(q_vec, q_term_vec);
        r_vec = AbstractVector{T}[copy(r_stage) for _ = 1:Np]

        # --- Dummy Constraints (Replace with config loading) ---
        # Example: Box constraints on physical inputs u_phys
        u_phys_max = 10 # Max physical input value
        # Constraint matrices E, F apply to augmented state z and input u'
        # To constrain u_phys, you typically constrain the first m_phys elements of u'
        # This assumes the first block of u' corresponds to u_phys (might not be true if u' is complex)
        # If u' = vec(u_e*u_m'), constraining u_phys directly is hard in this dense formulation.
        # We might constrain u' elements instead, or use a different MPC formulation.
        # For simplicity, let's constrain the elements of u' directly.
        u_prime_max = u_phys_max # Placeholder: Assume max value applies to elements of u'
        num_constraints_per_step = 2 * m_aug # Upper and lower bounds for each element of u'
        E_stage = spzeros(T, num_constraints_per_step, N_aug)
        F_stage = sparse(vcat(I(m_aug), -I(m_aug))) # [ I; -I ] * u' <= [max; -min]
        b_stage = ones(T, num_constraints_per_step) * u_prime_max # Assumes symmetric bounds [-max, max]

        E_term = spzeros(T, 0, N_aug); b_term = zeros(T, 0); # No terminal state constraints example

        E_vec = AbstractMatrix{T}[copy(E_stage) for _ = 1:Np]; push!(E_vec, E_term);
        F_vec = AbstractMatrix{T}[copy(F_stage) for _ = 1:Np]
        b_vec = AbstractVector{T}[copy(b_stage) for _ = 1:Np]; push!(b_vec, b_term);
        # --- End Dummy Costs & Constraints ---


        # --- Setup Dense MPC Controller (Offline Computation) ---
        @info "Setting up DenseMPCLib controller..."
        dense_controller = DenseMPCLib.setup_dense_mpc(
            A_full, B_full, Np,
            Q_vec, R_vec, q_vec, r_vec,
            E_vec, F_vec, b_vec,
            m_phys # Pass physical input dimension
        )
        @info "DenseMPCLib controller setup complete."

        # --- Initialize History Buffers ---
        # Buffer size needs to accommodate the required delay
        # x_history stores x_k down to x_{k-DELAY} -> DELAY+1 columns
        # u_phys_history stores u_{k-1} down to u_{k-DELAY} -> DELAY columns
        x_hist_buffer = CircularArrayBuffer{Float64}(N_X, DELAY + 1)
        u_phys_hist_buffer = CircularArrayBuffer{Float64}(m_phys, DELAY) # Size DELAY

        # Fill buffers with zeros initially (or load from data if available)
        fill!(x_hist_buffer, 0.0)
        fill!(u_phys_hist_buffer, 0.0)

        # Initialize last computed input
        last_u_phys = zeros(Float64, m_phys)

        # --- Create and Return Controller Instance ---
        new(Np, N_X, N_U, m_phys, DELAY,
            dense_controller,
            vec(mean_vec), # Ensure mean_vec is a vector
            x_hist_buffer,
            u_phys_hist_buffer,
            last_u_phys)
    end
end


# --- Implement Controller Interface Functions ---

"""
Computes the control action for the DMDC_MPCController.
"""
function compute_control_action(controller::DMDC_MPCController, time::Float64, state::Dict{Symbol, Any})
    # 1. Estimate Current Physical State (x_k)
    x_k = estimate_physical_state(state, controller.N_X)

    # 2. Update History Buffers
    push!(controller.x_history, x_k)
    if controller.DELAY > 0
        # Push the *physical* input computed in the *previous* step
        push!(controller.u_phys_history, controller.last_computed_u_phys)
    end

    # 3. Construct Current Augmented State (z_k_unnormalized)
    # ... (Keep existing code for constructing z_k_unnormalized) ...
    x_history_matrix = Matrix{Float64}(undef, controller.N_X, controller.DELAY + 1)
    for i in 0:controller.DELAY
        idx = controller.DELAY + 1 - i
        if idx <= length(controller.x_history)
            x_history_matrix[:, i+1] = controller.x_history[:, idx]
        else
            x_history_matrix[:, i+1] .= 0.0
            @warn "[MPC] x_history buffer not full at step k-$(i). Using zeros."
        end
    end

    u_prime_history_matrix = Matrix{Float64}(undef, controller.N_U, controller.DELAY)
    if controller.DELAY > 0
        for i in 1:controller.DELAY
            idx = controller.DELAY + 1 - i
            if idx <= length(controller.u_phys_history)
                u_phys_past = controller.u_phys_history[:, idx]
                # Preprocess past physical input to get past model input u'
                u_prime_history_matrix[:, i] = preprocess_input(u_phys_past, controller.N_U)
            else
                u_prime_history_matrix[:, i] .= 0.0
                @warn "[MPC] u_phys_history buffer not full at step k-$(i). Using zeros for u'."
            end
        end
    end
    z_k_unnormalized = vcat(vec(x_history_matrix), vec(u_prime_history_matrix))


    # 4. Normalize Augmented State
    mean_state_part = repeat(controller.mean_vec, controller.DELAY + 1)
    zero_input_part = zeros(Float64, controller.N_U * controller.DELAY)
    full_mean_vector = vcat(mean_state_part, zero_input_part)
    normalized_state_z0 = z_k_unnormalized - full_mean_vector

    # 5. Compute Optimal *Model* Control Action using DenseMPCLib
    # Calls the MODIFIED function which returns the 16-dim u'_0
    u_prime_optimal_k0, status = DenseMPCLib.compute_control_action(controller.dense_controller, normalized_state_z0)

    local u_phys_to_apply::Vector{Float64}

    # --- Logging and Status Handling ---
    if status == JuMP.OPTIMAL || status == JuMP.ALMOST_OPTIMAL
        # 6. Project the optimal 16-dim model input onto the outer product manifold
        # Use m=4, n=4 based on your preprocessing step (u_e is 4, u_m is 4)
        # Assumes find_closest_outer_product is available in the scope
        # (e.g., via include or using statement)
        m_outer = 4 # Dimension of u_e
        n_outer = 4 # Dimension of u_m
        if length(u_prime_optimal_k0) == m_outer * n_outer
             # Use the projection function you provided
             u_e_proj, u_m_proj = find_closest_outer_product(u_prime_optimal_k0, m_outer, n_outer)

             # The physical control to apply is the concatenation [u_e; u_m]
             u_phys_to_apply = [u_e_proj; u_m_proj] # Should be size 8

             @info "[MPC] Time $time: Projected u'_0 to get u_e=$(round.(u_e_proj, digits=2)), u_m=$(round.(u_m_proj, digits=2))"

             # Optional: Clamp u_e_proj, u_m_proj to physical limits if necessary
             # u_phys_to_apply = clamp.(u_phys_to_apply, min_limit, max_limit)

        else
            @error "[MPC] Optimal model input dimension ($(length(u_prime_optimal_k0))) does not match expected outer product dim ($(m_outer * n_outer)). Using fallback."
            u_phys_to_apply = controller.last_computed_u_phys # Fallback
        end

    else
        @warn "[MPC] QP solve failed at time $time with status: $status. Outputting last computed physical input."
        # Return the last successfully computed physical input as a fallback
        u_phys_to_apply = controller.last_computed_u_phys
    end

    # 7. Store the *physical* input computed for the next step's history
    controller.last_computed_u_phys = u_phys_to_apply

    # 8. Return the physical input to be applied to the system
    return u_phys_to_apply
end


"""
Updates the internal state of the controller (if any).
For this specific controller, the state is managed via history buffers
updated in `compute_control_action`. If an explicit observer were used,
it would be updated here based on `new_state_data`.
"""
function update_controller_state!(controller::DMDC_MPCController, time::Float64, new_state_data::Dict{Symbol, Any})
    # Currently, history updates happen in compute_control_action just before they are needed.
    # If you add an observer (e.g., Kalman filter) that needs updating
    # as soon as new data arrives, implement that logic here.
    # Example: update_observer!(controller.observer, time, new_state_data)
    return nothing
end
