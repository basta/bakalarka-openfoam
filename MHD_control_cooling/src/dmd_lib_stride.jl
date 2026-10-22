#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays and a specified stride.
It constructs an augmented state-space model (A_full, B_full)
and compares its simulation results against the original DMDc formulation
and the actual data.

Key Functions:
- `load_and_preprocess_data`: Loads data from JLD2, preprocesses control input U.
- `delay_embed`: Creates delay-embedded matrices with stride.
- `calculate_dmdc_matrices`: Computes DMDc system matrices (As, Bs, M) from data.
- `construct_augmented_matrices`: Builds the augmented A_full, B_full matrices.
- `create_abc_system_from_data`: Orchestrates data loading, DMDc, and ABC matrix construction.
- `create_abc_init_state`, `create_orig_init_state`: Generate initial states for simulations.
- `simulate_abc_model`, `simulate_orig_model`: Run simulations for both model types.
- `plot_comparison`: Generates and saves a comparison plot.
- `main`: Executes the full workflow.
=#
using JLD2
using LinearAlgebra
using Statistics
using Plots # Using Plots.jl for simplicity. Add Pkg.add("Plots") if needed.


# --- Core Data Processing and Model Creation Functions ---

"""
Loads X and U data, preprocesses U.
"""
function load_and_preprocess_data(data_path::String)
    println("Loading data from $data_path...")
    try
        X_data_raw = jldopen(data_path, "r")["X_data"]
        # Assuming RELEVANT_STATES is defined globally or passed as an argument
        # For this example, RELEVANT_STATES needs to be accessible
        X_data = X_data_raw[RELEVANT_STATES, 2:end]


        U_data_raw = jldopen(data_path, "r")["U_data"]

        println("Preprocessing U data...")
        U = let
            cols = []
            for u_col in eachcol(U_data_raw)
                if size(u_col, 1) == 8
                    push!(cols, vec(u_col[1:4] * u_col[5:8]'))
                else
                    # Ensure U_data_raw is defined or handle this case
                    println("Warning: Unexpected size for column in U_data_raw. Skipping. Column size: $(size(u_col))")
                end
            end
            hcat(cols...)
        end

        N_X = size(X_data, 1)
        N_U = size(U, 1)
        TOTAL_SAMPLES = size(X_data, 2) # Based on preprocessed X_data
        println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")

        return X_data, U, N_X, N_U, TOTAL_SAMPLES
    catch e
        println("Error loading data from $data_path: $e")
        println("Please ensure the JLD2 file exists, the path is correct, and RELEVANT_STATES is defined.")
        rethrow(e) # Propagate the error
    end
end

"""
Creates a delay-embedded matrix from time series data with a specified stride.
`data`: Matrix of time series data (features x samples).
`num_delays`: Number of delay steps (D). Results in D+1 blocks.
`stride`: Step size between delayed samples.
Output `embedded_data[:,t]` is `[data[:, k]; data[:, k-stride]; ...; data[:, k-num_delays*stride]]`
where `k = t + num_delays*stride`.
"""
function delay_embed(data::AbstractMatrix, num_delays::Int, stride::Int)
    n_features, n_samples = size(data)

    if num_delays < 0
        error("Number of delays must be non-negative.")
    end
    if stride <= 0
        error("Stride must be positive.")
    end

    if num_delays == 0
        # For no delays, return a copy of the data.
        # The length of this data will be handled by the calling function's slicing.
        return copy(data)
    end

    actual_lag_span = num_delays * stride
    if n_samples <= actual_lag_span
        error("Not enough samples ($n_samples) to create embedding with $num_delays delays and stride $stride (requires history of $actual_lag_span samples). Minimum $actual_lag_span + 1 samples needed.")
    end

    n_embedded_samples = n_samples - actual_lag_span
    embedded_dim = n_features * (num_delays + 1)
    embedded_data = similar(data, embedded_dim, n_embedded_samples)

    for t in 1:n_embedded_samples # t is the index for columns in embedded_data
        for d in 0:num_delays     # d is the delay block index (0 for most recent in this embedding structure, num_delays for oldest)
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features

            # source_col corresponds to data[:, t_anchor - d*stride]
            # where t_anchor is the time index of the most recent sample in the current embedded vector.
            # t_anchor = t + actual_lag_span (or t + num_delays*stride)
            source_col = t + (num_delays - d) * stride

            embedded_data[row_start:row_end, t] = data[:, source_col]
        end
    end
    return embedded_data
end

function input_lift(U::AbstractMatrix)
    # Placeholder for any specific input lifting/transformation if needed
    U_lifted = U
    return U_lifted
end

"""
Calculates DMDc matrices (As, Bs, M) from data with specified delay and stride.
`delay`: Number of delay terms (D).
`stride`: Stride for delays.
"""
function calculate_dmdc_matrices(X_data::AbstractMatrix, U::AbstractMatrix, delay::Int, stride::Int, N_X::Int, N_U::Int, TOTAL_SAMPLES::Int)
    println("Creating DMDc model with delay $delay and stride $stride...")

    # max_delay_span is the total time window covered by the delays
    max_delay_span = delay * stride

    # Prepare data indices
    # X_current is X(k+1), so k starts from max_delay_span + 1
    # Thus, k+1 starts from max_delay_span + 2
    start_idx_curr = max_delay_span + 2
    end_idx_curr = TOTAL_SAMPLES

    # Data segments for embedding: X_data up to k_max, U_data up to k_max
    # k_max is TOTAL_SAMPLES - 1 because X_current goes up to TOTAL_SAMPLES (which is X(k_max+1))
    start_idx_hist = 1
    end_idx_hist = TOTAL_SAMPLES - 1

    if start_idx_curr > end_idx_curr
        error("Not enough data points to form X_current. start_idx_curr=$start_idx_curr, end_idx_curr=$end_idx_curr. Need TOTAL_SAMPLES >= $(max_delay_span + 2). Current TOTAL_SAMPLES=$TOTAL_SAMPLES.")
    end
    if end_idx_hist < start_idx_hist + max_delay_span && delay > 0 # only check if delays are present
        error("Not enough data points in history segment for embedding. end_idx_hist=$end_idx_hist, required history span=$max_delay_span. Need end_idx_hist >= $(start_idx_hist + max_delay_span).")
    end


    U_lifted = input_lift(U) # Apply any lifting to U before embedding

    X_current = X_data[:, start_idx_curr:end_idx_curr]

    # Create segments for embedding. These segments must contain enough points for the delay embedding.
    X_hist_segment = X_data[:, start_idx_hist:end_idx_hist]
    U_hist_segment = U_lifted[:, start_idx_hist:end_idx_hist]

    X_delayed = delay_embed(X_hist_segment, delay, stride)
    U_delayed = delay_embed(U_hist_segment, delay, stride)

    # Align dimensions:
    # X_current has (end_idx_curr - start_idx_curr + 1) columns
    # X_delayed has (size(X_hist_segment,2) - max_delay_span) columns if delay > 0
    # If delay == 0, X_delayed is X_hist_segment, U_delayed is U_hist_segment

    local common_len
    if delay == 0
        # X_delayed and U_delayed are just X_hist_segment and U_hist_segment
        # X_current starts at index 2. X_hist_segment starts at index 1.
        # We need to align X_current (X(k+1)) with X_delayed (X(k))
        # Number of X(k+1) samples: size(X_current,2)
        # Number of X(k) samples available from X_hist_segment: size(X_hist_segment,2)
        # X_current[:,i] corresponds to X_data[:, i + start_idx_curr -1]
        # X_delayed[:,j] corresponds to X_data[:, j + start_idx_hist -1]
        # We need (i + start_idx_curr -1) = (j + start_idx_hist -1) + 1
        # i + max_delay_span + 2 -1 = j + 1 - 1 + 1 => i + 1 = j + 1 (since max_delay_span=0) => i=j
        common_len = min(size(X_current, 2), size(X_delayed, 2) - (start_idx_curr - (start_idx_hist + 1)))
        common_len = min(size(X_current, 2), size(X_delayed, 2) - 1) # X_delayed has one more point at the start
        common_len = min(common_len, size(U_delayed, 2) - 1)

        X_next_dmd = X_current[:, 1:common_len]
        # X_embed_dmd should correspond to X(k), so it's X_delayed shifted by 1
        X_embed_dmd = X_delayed[:, (start_idx_curr-1)-start_idx_hist+1:(start_idx_curr-1)-start_idx_hist+common_len]
        U_embed_dmd = U_delayed[:, (start_idx_curr-1)-start_idx_hist+1:(start_idx_curr-1)-start_idx_hist+common_len]

    else
        # For delay > 0, delay_embed handles the length based on max_delay_span
        # X_delayed[:,t] corresponds to embedding of X_data at time t_orig = t + start_idx_hist -1 + max_delay_span
        # X_current[:,t'] corresponds to X_data at time t'_orig = t' + start_idx_curr -1
        # We need t'_orig = t_orig + 1
        # t' + max_delay_span + 2 -1 = t + max_delay_span + 1 => t' = t
        common_len = min(size(X_current, 2), size(X_delayed, 2), size(U_delayed, 2))
        X_next_dmd = X_current[:, 1:common_len]
        X_embed_dmd = X_delayed[:, 1:common_len]
        U_embed_dmd = U_delayed[:, 1:common_len]
    end


    if common_len <= 0
        error("Resulting common length for DMDc matrices is non-positive ($common_len). Check data lengths and indices. X_current: $(size(X_current,2)), X_delayed: $(size(X_delayed,2)), U_delayed: $(size(U_delayed,2))")
    end

    #==== LIFITING EXAMPLE (applied after embedding) ====#
    # U_embed_dmd = abs.(U_embed_dmd) # Example: apply abs value lifting

    # Center X data
    local mean_X
    if NORMALIZE # NORMALIZE must be defined globally
        mean_X = mean(X_data, dims=2)
    else
        mean_X = zeros(size(X_data, 1), 1)
    end
    X_next_centered = X_next_dmd .- mean_X
    # X_embed_dmd has (delay+1) blocks of X states
    X_embed_centered = X_embed_dmd .- repeat(mean_X, delay + 1, 1)


    # Combined state for regression
    Omega = [X_embed_centered; U_embed_dmd]

    # Solve for the combined dynamics matrix M
    println("Performing DMDc regression...")
    M = X_next_centered * pinv(Omega) # pinv for robustness
    println("DMDc matrix M size: $(size(M))")

    # Extract A and B matrices
    # A_matrices multiply X_embed_centered, B_matrices multiply U_embed_dmd
    # X_embed_centered has N_X * (delay+1) rows
    A_matrices_flat = M[:, 1:N_X*(delay+1)]
    B_matrices_flat = M[:, N_X*(delay+1)+1:end]

    # Reshape into list of As and Bs
    # As[1] multiplies x(k), As[2] multiplies x(k-stride), ..., As[delay+1] multiplies x(k-delay*stride)
    # Bs[1] multiplies u(k), Bs[2] multiplies u(k-stride), ..., Bs[delay+1] multiplies u(k-delay*stride)
    As = [A_matrices_flat[:, i*N_X+1:(i+1)*N_X] for i in 0:delay]
    Bs = [B_matrices_flat[:, i*N_U+1:(i+1)*N_U] for i in 0:delay]

    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, mean_X
end


"""
Constructs augmented state-space matrices A_full, B_full.
`delay` is the number of delay terms (D). Stride is implicitly handled by As, Bs.
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, delay::Int)
    println("Constructing augmented state-space model (A_full, B_full)...")
    D = delay # Number of delay stages (e.g., D=0 means x(k) only, D=1 means x(k), x(k-stride))

    # Augmented state: [x_k; x_{k-s}; ...; x_{k-Ds}; u_{k-s}; ...; u_{k-Ds}]
    # Size of x-part: N_X * (D+1)
    # Size of u-part (delayed inputs in state): N_U * D
    N_STATES_AUG = N_X * (D + 1) + N_U * D

    # Construct T matrix (top row of A_full predicting x_{k+1})
    # x_{k+1} = As[1]x_k + As[2]x_{k-s} + ... + As[D+1]x_{k-Ds} +
    #             Bs[2]u_{k-s} + ... + Bs[D+1]u_{k-Ds}
    # (Bs[1]u_k is handled by B_full)
    T = zeros(N_X, N_STATES_AUG)
    for i in 1:(D+1) # Corresponds to As[1]...As[D+1]
        T[:, (i-1)*N_X+1:i*N_X] = As[i]
    end
    if D > 0
        u_del_offset_in_state = N_X * (D + 1) # Start index of u_{k-s} block in augmented state
        for i in 1:D # Corresponds to Bs[2]...Bs[D+1] for u_{k-s} ... u_{k-Ds}
            # Bs[i+1] multiplies u_{k - i*stride}
            T[:, u_del_offset_in_state+(i-1)*N_U+1:u_del_offset_in_state+i*N_U] = Bs[i+1]
        end
    end

    # Construct Sx matrix (state shifting for x part)
    # [x_k; x_{k-s}; ... ] -> [x_{k+1}; x_k; x_{k-s}; ...]
    Sx = zeros(N_X * D, N_X * (D + 1))
    if D > 0
        # x_k becomes the first delayed block x_{k+1 - s} (conceptually)
        # The augmented state stores x_k, x_{k-s}, ...
        # So the new x_k goes to the first x_del block, old x_k goes to second x_del block etc.
        # This means Identity shifts the blocks down.
        # x_k (from current state) becomes x_{k-s} (in next augmented state at this position)
        # x_{k-s} (from current state) becomes x_{k-2s} (...)
        for i in 1:D # For each of the D delayed x-blocks in the new state
            # The i-th delayed block is filled by the (i-1)-th block from old state
            Sx[(i-1)*N_X+1:i*N_X, (i-1)*N_X+1:i*N_X] = I(N_X)
        end
    end


    # Construct Su matrix (input shifting for u part in state)
    # [u_{k-s}; u_{k-2s}; ...] -> [u_k; u_{k-s}; ...] (u_k comes from B_full's effect on this part)
    Su = zeros(N_U * D, N_U * D)
    if D > 0
        # u_{k-s} (from current state) becomes u_{k-2s} (in next augmented state at this position)
        # This means Identity shifts the blocks down.
        # The first block u_{k-s} is filled by current u_k via B_full.
        if D > 1 # Need to shift u_{k-s} to u_{k-2s} etc.
            for i in 1:(D-1) # For u_{k-(i+1)s} block in new state
                # It's filled by u_{k-is} from old state's delayed u part
                Su[i*N_U+1:(i+1)*N_U, (i-1)*N_U+1:i*N_U] = I(N_U)
            end
        end
    end


    # Assemble A_full
    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    A_full[1:N_X, :] = T # Dynamics for x_k
    if D > 0
        # Shift x states
        A_full[N_X+1:N_X*(D+1), 1:N_X*(D+1)] = Sx[1:N_X*D, 1:N_X*(D+1)]

        # Shift u states (delayed part)
        u_del_start_idx_in_state = N_X * (D + 1) + 1
        A_full[u_del_start_idx_in_state:end, u_del_start_idx_in_state:end] = Su
    end

    # Assemble B_full (for current input u_k)
    B_full = zeros(N_STATES_AUG, N_U)
    B_full[1:N_X, :] = Bs[1] # Effect of u_k on x_{k+1}
    if D > 0
        # Current u_k becomes the first delayed input u_{k-s} (conceptually for next step)
        # This means u_k populates the block for u_{k-s} in the augmented state.
        u_first_del_start_idx = N_X * (D + 1) + 1
        B_full[u_first_del_start_idx:u_first_del_start_idx+N_U-1, :] = I(N_U)
    end

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end


"""
Main function to create the ABC system from data.
"""
function create_abc_system_from_data(data_path::String, delay::Int, stride::Int)
    X_data, U, N_X, N_U, TOTAL_SAMPLES = load_and_preprocess_data(data_path)
    As, Bs, M, mean_X = calculate_dmdc_matrices(X_data, U, delay, stride, N_X, N_U, TOTAL_SAMPLES)
    A_full, B_full = construct_augmented_matrices(As, Bs, N_X, N_U, delay)

    return A_full, B_full, As, Bs, M, vec(mean_X), N_X, N_U, X_data, U, TOTAL_SAMPLES
end

# --- Simulation and Initial State Functions ---

"""
Simulates the ABC model (augmented state-space).
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state::AbstractVector, inputs::AbstractMatrix, sim_length::Int)
    n_state_vars_aug = size(A_full, 1)
    _, n_input_steps = size(inputs) # n_inputs

    if n_input_steps < sim_length
        error("Not enough input data provided for simulation ($sim_length steps needed, $n_input_steps provided).")
    end

    states_aug = similar(initial_state, n_state_vars_aug, sim_length)
    current_state_aug = copy(initial_state)

    # Define the discrete-time update function
    abc_fn(state_aug, input_u_k) = A_full * state_aug + B_full * input_u_k

    for k in 1:sim_length
        states_aug[:, k] = current_state_aug
        current_state_aug = abc_fn(current_state_aug, inputs[:, k])
    end
    return states_aug # Returns the full augmented state trajectory
end

"""
Simulates the original DMDc model formulation using the M matrix.
`delay` is number of delay terms (D). Stride is implicit in M.
"""
function simulate_orig_model(initial_state_orig_model::AbstractVector,
    inputs::AbstractMatrix, # U_k, U_{k+1}, ...
    M_dmdc::AbstractMatrix,
    sim_length::Int,
    delay::Int, # D
    n_x::Int,
    n_u::Int)

    # initial_state_orig_model should be [x_embed_centered_k0; u_embed_k0]
    # where x_embed_centered_k0 = [x_k0_c; x_{k0-s}_c; ...; x_{k0-Ds}_c]
    # and u_embed_k0 = [u_k0; u_{k0-s}; ...; u_{k0-Ds}]
    # The simulation loop will update this structure.

    states_out_centered = zeros(n_x, sim_length) # Stores x_{k+1}_centered

    # current_full_omega_k represents [X_embed_centered(k); U_embed(k)]
    current_full_omega_k = copy(initial_state_orig_model)

    # The M_dmdc matrix expects Omega = [X_embed_centered; U_embed_dmd]
    # X_embed_centered has (delay+1)*n_x elements
    # U_embed_dmd has (delay+1)*n_u elements

    for k_sim_step in 1:sim_length
        # current_full_omega_k is [x_k_embed_c; u_k_embed]
        x_next_centered = M_dmdc * current_full_omega_k
        states_out_centered[:, k_sim_step] = x_next_centered

        # Update current_full_omega_k for the next step (k+1)
        # It needs to become [X_embed_centered(k+1); U_embed(k+1)]

        # 1. Update X_embed_centered part: [x_{k+1}_c; x_k_c; ...; x_{k+1-Ds}_c]
        new_x_embed_centered = zeros(n_x * (delay + 1))
        new_x_embed_centered[1:n_x] = x_next_centered # This is x_{k+1}_c
        if delay > 0
            # Shift old x_embed contents
            # Old x_embed was [x_k_c; x_{k-s}_c; ...; x_{k-Ds}_c]
            # We need [x_k_c; ...; x_{k-(D-1)s}_c] for the rest of new_x_embed_centered
            old_x_part_for_shift = current_full_omega_k[1:n_x*delay]
            new_x_embed_centered[n_x+1:end] = old_x_part_for_shift
        end

        # 2. Update U_embed part: [u_{k+1}; u_k; ...; u_{k+1-Ds}]
        new_u_embed = zeros(n_u * (delay + 1))
        if k_sim_step + 1 <= size(inputs, 2)
            new_u_embed[1:n_u] = inputs[:, k_sim_step+1] # This is u_{k+1}
        else
            # Handle case where we run out of inputs for u_{k+1} (e.g., predict one step beyond given inputs)
            # For simplicity, assume zero or repeat last input if needed for prediction
            println("Warning: Ran out of inputs for U_embed in simulate_orig_model. Using zeros for u_k+1.")
            new_u_embed[1:n_u] = zeros(n_u)
        end

        if delay > 0
            # Shift old u_embed contents
            # Old U_embed was [u_k; u_{k-s}; ...; u_{k-Ds}]
            # We need [u_k; ...; u_{k-(D-1)s}] for the rest of new_u_embed
            old_u_part_for_shift = current_full_omega_k[n_x*(delay+1)+1:n_x*(delay+1)+n_u*delay]
            new_u_embed[n_u+1:end] = old_u_part_for_shift
        end
        current_full_omega_k = [new_x_embed_centered; new_u_embed]
    end
    return states_out_centered
end


"""
Creates the initial state vector for the ABC model.
State: [x(k); x(k-s); ...; x(k-Ds); u(k-s); ...; u(k-Ds)] (all x terms centered)
`delay` is D. `start_time` is k.
"""
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, stride::Int, mean_x::AbstractVector, N_X::Int, N_U::Int)
    n_states_aug = N_X * (delay + 1) + N_U * delay
    init_state_aug = zeros(n_states_aug)

    min_required_time_x = start_time - delay * stride
    min_required_time_u = start_time - delay * stride # for u(k-Ds) if D > 0, or start_time-stride for u(k-s)

    if min_required_time_x < 1 && delay > 0 # only matters if there are delays
        error("Cannot create initial ABC state at time $start_time with delay $delay, stride $stride for X: requires data back to index $min_required_time_x, but data starts at 1.")
    end
    if delay > 0 && min_required_time_u < 1
        error("Cannot create initial ABC state at time $start_time with delay $delay, stride $stride for U: requires data back to index $min_required_time_u, but data starts at 1.")
    end
    if start_time > size(X, 2) || start_time > size(U, 2)
        error("start_time $start_time is out of bounds for X (max $(size(X,2))) or U (max $(size(U,2))).")
    end


    # Current state x(k) centered
    init_state_aug[1:N_X] = X[:, start_time] .- mean_x

    # Delayed states x(k-s) ... x(k-Ds) centered
    for i in 1:delay # i from 1 to D
        state_idx_start = N_X * i + 1
        state_idx_end = N_X * (i + 1)
        hist_time_idx = start_time - i * stride
        if hist_time_idx < 1
            error("Index $hist_time_idx out of bounds for X data in create_abc_init_state.")
        end
        init_state_aug[state_idx_start:state_idx_end] = X[:, hist_time_idx] .- mean_x
    end

    # Delayed inputs u(k-s) ... u(k-Ds) for the augmented state part
    if delay > 0
        u_del_aug_part_start_idx = N_X * (delay + 1) + 1
        for i in 1:delay # i from 1 to D, for u(k-i*s)
            input_idx_start = u_del_aug_part_start_idx + N_U * (i - 1)
            input_idx_end = input_idx_start + N_U - 1
            hist_time_idx = start_time - i * stride
            if hist_time_idx < 1
                error("Index $hist_time_idx out of bounds for U data in create_abc_init_state.")
            end
            init_state_aug[input_idx_start:input_idx_end] = U[:, hist_time_idx]
        end
    end
    return init_state_aug
end

"""
Creates the initial state vector for the original DMDc model simulation (M matrix).
Initial Omega: [x_embed_centered(k); u_embed(k)]
x_embed_centered(k) = [x(k)_c; x(k-s)_c; ...; x(k-Ds)_c]
u_embed(k)          = [u(k);   u(k-s);   ...; u(k-Ds)]
`delay` is D. `start_time` is k.
"""
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, stride::Int, mean_x::AbstractVector, N_X::Int, N_U::Int)

    len_x_embed = N_X * (delay + 1)
    len_u_embed = N_U * (delay + 1)
    init_omega = zeros(len_x_embed + len_u_embed)

    min_required_time = start_time - delay * stride
    if min_required_time < 1 && delay > 0
        error("Cannot create initial original model state at time $start_time with delay $delay, stride $stride: requires data back to index $min_required_time, but data starts at 1.")
    end
    if start_time > size(X, 2) || start_time > size(U, 2)
        error("start_time $start_time is out of bounds for X (max $(size(X,2))) or U (max $(size(U,2))).")
    end

    # Fill X_embed_centered part
    for i in 0:delay # i from 0 to D, for x(k-i*s)
        idx_start = i * N_X + 1
        idx_end = (i + 1) * N_X
        hist_time_idx = start_time - i * stride
        if hist_time_idx < 1
            error("Index $hist_time_idx out of bounds for X data in create_orig_init_state.")
        end
        init_omega[idx_start:idx_end] = X[:, hist_time_idx] .- mean_x
    end

    # Fill U_embed part
    u_embed_part_offset = len_x_embed
    for i in 0:delay # i from 0 to D, for u(k-i*s)
        idx_start = u_embed_part_offset + i * N_U + 1
        idx_end = u_embed_part_offset + (i + 1) * N_U
        hist_time_idx = start_time - i * stride
        if hist_time_idx < 1
            error("Index $hist_time_idx out of bounds for U data in create_orig_init_state.")
        end
        init_omega[idx_start:idx_end] = U[:, hist_time_idx]
    end
    return init_omega
end


"""
Plots the comparison between actual data and model simulations.
"""
function plot_comparison(actual_data::AbstractMatrix, states_abc::AbstractMatrix, states_orig::AbstractMatrix,
    plot_state_idx::Int, sim_start_time::Int, sim_len::Int, delay::Int, stride::Int)

    println("Plotting results...")
    time_axis = 1:sim_len

    p = plot(time_axis, actual_data[plot_state_idx, 1:sim_len],
        label="Actual Data (X$plot_state_idx)",
        linewidth=2, color=:black)
    plot!(p, time_axis, states_abc[plot_state_idx, 1:sim_len],
        label="ABC Model (X$plot_state_idx)",
        linewidth=1.5, linestyle=:dash, color=:red)
    plot!(p, time_axis, states_orig[plot_state_idx, 1:sim_len],
        label="Original DMDc Model (X$plot_state_idx)",
        linewidth=1.5, linestyle=:dot, color=:blue)

    title!("Model Comparison (Delay=$delay, Stride=$stride)")
    xlabel!("Time Step (relative to data index $sim_start_time)")
    ylabel!("State Value")
    plot!(legend=:outertopright)

    figname = "dmdc_comparison_delay$(delay)_stride$(stride)_start$(sim_start_time).png"
    try
        savefig(p, figname)
        println("Plot saved to $figname")
    catch e
        println("Error saving plot: $e. Make sure a plotting backend is working (e.g., GR).")
    end
end


# --- Configuration ---
# Define RELEVANT_STATES before it's used in load_and_preprocess_data if not passed as arg
const RELEVANT_STATES = [7] # Example: use the first state. Adjust as needed.
const DATA_PATH = "./data.jld2" # Adjust path if needed. Create a dummy data.jld2 for testing.
const DELAY = 100          # Number of delay steps (D)
const STRIDE = 4         # Stride for delay embeddings (stride=1 is standard no-stride)
const START_IDX = 700     # Starting index in X_data for simulation comparison (1-based)
const SIM_LEN = 50       # Default length of the simulation
const PLOT_STATE_IDX = 1 # Which state component to plot (1-based for X_data)
const NORMALIZE = false   # Whether to center X data around its mean

# --- Main Execution ---
function main()
    # Create dummy data if DATA_PATH doesn't exist, for self-contained testing
    if !isfile(DATA_PATH)
        println("Creating dummy data file: $DATA_PATH")
        _X_data = randn(max(RELEVANT_STATES...), 100 + maximum(RELEVANT_STATES) * 2) # ensure enough samples and states
        _U_data = randn(8, 100) # each column is 8x1 for the preprocessing U logic
        jldsave(DATA_PATH; X_data=_X_data, U_data=_U_data)
        println("Dummy data created. Rerun the script.")
        return
    end

    # 1. Create ABC System
    A_full, B_full, As, Bs, M_dmdc, mean_X_vec, N_X, N_U, X_data, U, TOTAL_SAMPLES =
        create_abc_system_from_data(DATA_PATH, DELAY, STRIDE)

    # 2. Simulation Setup
    # sim_start_time is the index in X_data (which is already X_data_raw[:, 2:end])
    # It's the 'k' for which we want to form the initial state.
    # We need DELAY*STRIDE history points before sim_start_time.
    # So, earliest valid sim_start_time is DELAY*STRIDE + 1.
    required_history_span = DELAY * STRIDE
    min_sim_start_time = required_history_span + 1

    sim_start_time = max(START_IDX, min_sim_start_time)

    actual_sim_len = SIM_LEN
    if sim_start_time + actual_sim_len - 1 > TOTAL_SAMPLES
        println("Warning: Simulation length ($actual_sim_len) starting at $sim_start_time exceeds data bounds (TOTAL_SAMPLES=$TOTAL_SAMPLES).")
        actual_sim_len = TOTAL_SAMPLES - sim_start_time + 1
        if actual_sim_len <= 0
            println("Error: Not enough data to run simulation even for 1 step. sim_start_time=$sim_start_time, TOTAL_SAMPLES=$TOTAL_SAMPLES")
            return
        end
        println("Adjusting sim_len to $actual_sim_len.")
    end

    # Inputs for simulation: U[:, k_sim], U[:, k_sim+1], ..., U[:, k_sim+L-1]
    sim_inputs = U[:, sim_start_time:sim_start_time+actual_sim_len-1]

    # 3. Create Initial States
    # initial states are for time 'sim_start_time'
    init_state_abc = create_abc_init_state(X_data, U, sim_start_time, DELAY, STRIDE, mean_X_vec, N_X, N_U)
    init_state_orig = create_orig_init_state(X_data, U, sim_start_time, DELAY, STRIDE, mean_X_vec, N_X, N_U)


    # 4. Run Simulations
    println("Running simulations...")
    # simulate_abc_model outputs augmented states, take the first N_X rows for x(k)
    states_abc_augmented_full = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, actual_sim_len)
    states_abc = states_abc_augmented_full[1:N_X, :] .+ mean_X_vec # Add mean back

    # simulate_orig_model outputs x_next_centered
    states_orig_centered = simulate_orig_model(init_state_orig, sim_inputs, M_dmdc, actual_sim_len, DELAY, N_X, N_U)
    states_orig = states_orig_centered .+ mean_X_vec # Add mean back

    # 5. Get Actual Data for Comparison
    # actual_data should be X_data[:, k_sim : k_sim+L-1]
    actual_data_comparison = X_data[:, sim_start_time:sim_start_time+actual_sim_len-1]
    println("Simulations complete.")

    # 6. Plot Results
    if PLOT_STATE_IDX > N_X || PLOT_STATE_IDX < 1
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is out of bounds for N_X=$N_X. Plotting first state.")
        plot_idx_to_use = 1
    else
        plot_idx_to_use = PLOT_STATE_IDX
    end
    plot_comparison(actual_data_comparison, states_abc, states_orig, plot_idx_to_use, sim_start_time, actual_sim_len, DELAY, STRIDE)

    println("Script finished.")

    # 7. Save results
    try
        jldsave("ABC.jld2"; compress=true,
            A=A_full,
            B=B_full,
            delay=DELAY,
            stride=STRIDE,
            mean_vec=mean_X_vec,
            M_dmdc=M_dmdc,
            N_X=N_X,
            N_U=N_U,
            TOTAL_SAMPLES=TOTAL_SAMPLES, # Total samples in preprocessed X_data
            RELEVANT_STATES=RELEVANT_STATES, # Save the const used
            sim_start_time_used=sim_start_time,
            sim_len_used=actual_sim_len
        )
        println("Results saved to ABC_stride.jld2")
    catch e
        println("Error saving results to JLD2 file: $e")
    end
end

if @isdefined(main)
    main()
end
