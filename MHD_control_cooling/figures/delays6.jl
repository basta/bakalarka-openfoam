#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays.
It constructs an augmented state-space model (A_full, B_full)
and compares its simulation results for various delay amounts against the actual data.

Key Functions:
- `load_and_preprocess_data`: Loads data from JLD2, preprocesses control input U.
- `delay_embed`: Creates delay-embedded matrices.
- `calculate_dmdc_matrices`: Computes DMDc system matrices (As, Bs, M) from data.
- `construct_augmented_matrices`: Builds the augmented A_full, B_full matrices.
- `create_abc_init_state`: Generate initial states for ABC model simulations.
- `simulate_abc_model`: Run simulations for the ABC model.
- `main_multi_delay_comparison`: Executes the workflow for multiple delays and plots comparison.
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
        # Apply RELEVANT_STATES selection and remove first time step if it's an offset
        X_data = X_data_raw[RELEVANT_STATES, 2:end]

        U_data_raw = jldopen(data_path, "r")["U_data"]

        println("Preprocessing U data...")
        U = let
            cols = []
            for u_col in eachcol(U_data_raw)
                if size(u_col, 1) == 8 # Specific to user's U_data structure
                    push!(cols, vec(u_col[1:4] * u_col[5:8]'))
                else
                    # Fallback if U_data column is not size 8 (e.g. already processed or different format)
                    # This assumes if not 8, it's already in the desired N_U x 1 format per column
                    # Adjust this logic if your U_data can have other raw formats
                    if length(RELEVANT_STATES) > 0 && size(X_data,1) > 0 && size(U_data_raw,1) == 16 # Heuristic for the example's U structure
                         println("Warning: U_data column size is not 8. Assuming it might be 16 (e.g. from lifted U). Attempting to use as is if N_U matches.")
                         push!(cols, vec(u_col)) # This might need adjustment based on actual N_U
                    else
                         push!(cols, vec(u_col))
                    end
                end
            end
            hcat(cols...)
        end

        N_X = size(X_data, 1)
        N_U = size(U, 1) # N_U is determined after preprocessing U
        TOTAL_SAMPLES = size(X_data, 2)

        # Ensure U has the same number of samples as X_data after its own preprocessing
        if size(U, 2) != TOTAL_SAMPLES
            min_samples = min(size(U,2), TOTAL_SAMPLES)
            println("Warning: X_data and U have different number of samples after preprocessing. Truncating to $min_samples samples.")
            X_data = X_data[:, 1:min_samples]
            U = U[:, 1:min_samples]
            TOTAL_SAMPLES = min_samples
        end

        println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")
        println("N_X = $N_X, N_U = $N_U, TOTAL_SAMPLES = $TOTAL_SAMPLES")

        if N_U == 0 && size(U,2) > 0 # If U became empty but had columns, N_U should be its row size
            _N_U_inferred = size(jldopen(data_path, "r")["U_data"], 1)
            if _N_U_inferred == 8 # Specific to user's U_data structure
                N_U = 16 # After outer product
            else
                N_U = _N_U_inferred # Or some other default/logic
            end
             println("Corrected N_U to $N_U based on raw U_data structure, as processed U was empty.")
        end


        return X_data, U, N_X, N_U, TOTAL_SAMPLES
    catch e
        println("Error loading data from $data_path: $e")
        println("Please ensure the JLD2 file exists, the path is correct, and RELEVANT_STATES are valid.")
        rethrow(e)
    end
end

"""
Creates a delay-embedded matrix from time series data.
"""
function delay_embed(data::AbstractMatrix, delay::Int)
    n_features, n_samples = size(data)
    if delay < 0
        error("Delay must be non-negative.")
    elseif delay == 0
        return data # No embedding needed for zero delay
    end

    if n_samples <= delay
        error("Not enough samples ($n_samples) to create embedding with delay ($delay). Need n_samples > delay.")
    end

    n_embedded_samples = n_samples - delay
    embedded_dim = n_features * (delay + 1)
    # Ensure data type is preserved, especially for Float types
    embedded_data = similar(data, eltype(data), embedded_dim, n_embedded_samples)


    for t in 1:n_embedded_samples
        for d in 0:delay
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features
            source_col = t + delay - d # Correct indexing for embedding
            embedded_data[row_start:row_end, t] = data[:, source_col]
        end
    end
    return embedded_data
end

function input_lift(U::AbstractMatrix)
    # Placeholder for any specific input lifting.
    # The original script had U_embed_dmd = abs.(U_embed_dmd) later.
    # If lifting is only on U_embed_dmd, it's fine. If U itself needs lifting before embedding:
    # U_lifted = abs.(U) # Example, if this was intended for the raw U
    U_lifted = U # Current script applies abs later
    return U_lifted
end

"""
Calculates DMDc matrices (As, Bs, M) from data.
"""
function calculate_dmdc_matrices(X_data::AbstractMatrix, U_original::AbstractMatrix, delay::Int, N_X::Int, N_U::Int, TOTAL_SAMPLES::Int)
    println("Creating DMDc model with delay $delay...")
    max_delay = delay # Using 'delay' directly for clarity

    # Prepare data indices for DMDc regression
    # X_next_dmd uses X_data from column max_delay + 2 onwards
    # X_embed_dmd and U_embed_dmd use data up to TOTAL_SAMPLES - 1 for history
    start_idx_curr = max_delay + 2 # Index for X_data to form X_next_dmd (i.e., x(k+1))
    end_idx_curr = TOTAL_SAMPLES   # Last index for X_data forming X_next_dmd

    # Data for embedding: X_data[:, 1 : TOTAL_SAMPLES-1], U_original[:, 1 : TOTAL_SAMPLES-1]
    # This provides sequence x(1)...x(T-1) and u(1)...u(T-1)
    # delay_embed will take care of constructing [x(k)...x(k-d)] from this.
    # The 'k' in x(k) for X_embed_dmd will correspond to time step of X_next_dmd's 'k+1'
    # So, if X_next_dmd starts at X_data[:, max_delay+2], its corresponding X_embed_dmd
    # should be built from X_data up to column max_delay+1.

    # X_prime (X_next_dmd) = X_data[:, (max_delay+2):TOTAL_SAMPLES]
    # X_dmd (for X_embed_dmd) = X_data[:, 1:(TOTAL_SAMPLES - (max_delay+2) + (max_delay+1))]
    # U_dmd (for U_embed_dmd) = U_original[:, 1:(TOTAL_SAMPLES - (max_delay+2) + (max_delay+1))]

    # Let's define the range for X_prime (target states)
    X_prime_range = (max_delay + 2):TOTAL_SAMPLES
    num_regression_samples = length(X_prime_range)

    if num_regression_samples <= 0
        error("Not enough data points (num_regression_samples=$num_regression_samples) to create DMDc model with delay $max_delay and data size $TOTAL_SAMPLES. Need TOTAL_SAMPLES > max_delay + 1.")
    end

    X_next_dmd = X_data[:, X_prime_range]

    # Construct Omega: [X_embed_centered; U_embed_dmd]
    # X_embed_dmd is [x(k); x(k-1); ...; x(k-D)]
    # U_embed_dmd is [u(k); u(k-1); ...; u(k-D)]
    # where 'k' corresponds to the time step *before* X_next_dmd.
    # So, if X_next_dmd(:,t) is X_data[:, X_prime_range[t]], then
    # X_embed_dmd needs history up to X_prime_range[t]-1
    # U_embed_dmd needs history up to X_prime_range[t]-1 (or current u if u(k) is used)

    # Simpler: define X and Y for regression: Y = G * Z
    # Y = X_data[:, (max_delay+1): (TOTAL_SAMPLES-1)]  (this is X_k+1, if X starts at 0)
    # Z = [ X_data[:, max_delay : (TOTAL_SAMPLES-2)] ; ... ; X_data[:, 0: (TOTAL_SAMPLES-2-max_delay)] ;
    #       U_data[:, max_delay : (TOTAL_SAMPLES-2)] ; ... ; U_data[:, 0: (TOTAL_SAMPLES-2-max_delay)] ]

    # Using the script's original indexing logic for Omega:
    # X_current = X_data[:, start_idx_curr:end_idx_curr] (This is X_next_dmd)
    # X_delayed_source = X_data[:, 1 : end_idx_curr-1]
    # U_delayed_source = U_original[:, 1 : end_idx_curr-1]

    # Corrected ranges for embedding based on X_next_dmd
    # X_next_dmd uses columns X_prime_range from X_data.
    # The state X_embed_dmd and U_embed_dmd should correspond to the time step *before* each column in X_next_dmd.
    # The last column used for history in X_embed_dmd is X_prime_range[end]-1.
    # The first column used for history in X_embed_dmd is X_prime_range[1]-1.
    # The number of samples for embedding source data is (X_prime_range[end]-1) - (X_prime_range[1]-1) + 1 + max_delay
    # This is complex. Let's use the common_len approach from original:

    # Data for X_delayed: X_data up to the time step *before* the first X_next_dmd
    # If X_next_dmd starts with X_data[:, max_delay+2], then X_delayed needs history ending at max_delay+1 for its first column.
    # This means X_delayed is built from X_data[:, 1 : (TOTAL_SAMPLES-1)]
    # And U_delayed is built from U_original[:, 1 : (TOTAL_SAMPLES-1)]

    X_embed_source = X_data[:, 1:(TOTAL_SAMPLES-1)]
    U_embed_source = U_original[:, 1:(TOTAL_SAMPLES-1)]

    X_delayed_full = delay_embed(X_embed_source, max_delay) # x(k)...x(k-D)
    U_delayed_full = delay_embed(U_embed_source, max_delay) # u(k)...u(k-D)

    # Align dimensions: X_next_dmd must align with columns of X_delayed_full and U_delayed_full
    # X_next_dmd has num_regression_samples columns.
    # X_delayed_full has size(X_embed_source,2) - max_delay columns.
    # size(X_embed_source,2) = TOTAL_SAMPLES-1.
    # So X_delayed_full has (TOTAL_SAMPLES-1) - max_delay columns.
    # num_regression_samples = TOTAL_SAMPLES - (max_delay+2) + 1 = TOTAL_SAMPLES - max_delay - 1.
    # These lengths match.

    X_embed_dmd = X_delayed_full
    U_embed_dmd = U_delayed_full

    #==== LIFITING HERE ====#
    U_embed_dmd = abs.(U_embed_dmd)


    local mean_X
    if NORMALIZE
        mean_X = mean(X_data, dims=2) # Mean of the original states (after RELEVANT_STATES)
    else
        mean_X = zeros(N_X, 1)
    end
    # Center X_next_dmd
    X_next_centered = X_next_dmd .- mean_X
    # Center X_embed_dmd: each block of N_X states needs to be centered
    X_embed_centered = similar(X_embed_dmd)
    for i in 0:max_delay
        row_idx = (i*N_X + 1) : ((i+1)*N_X)
        X_embed_centered[row_idx, :] = X_embed_dmd[row_idx, :] .- mean_X
    end

    # Combined state for regression
    Omega = [X_embed_centered; U_embed_dmd]

    # Solve for the combined dynamics matrix M
    println("Performing DMDc regression...")
    if size(Omega, 2) == 0
        error("Omega matrix is empty. Check data slicing and delay parameters. Omega columns: $(size(Omega,2))")
    end
    M = X_next_centered * pinv(Omega) # Y = M * Omega
    println("DMDc matrix M size: $(size(M))")

    # Extract A and B matrices
    # M = [A0 A1 ... AD B0 B1 ... BD]
    # Size of X_embed_centered part is N_X * (max_delay + 1)
    A_matrices_flat = M[:, 1:(N_X*(max_delay+1))]
    B_matrices_flat = M[:, (N_X*(max_delay+1)+1):end]

    As = [A_matrices_flat[:, i*N_X+1:(i+1)*N_X] for i in 0:max_delay]
    Bs = [B_matrices_flat[:, i*N_U+1:(i+1)*N_U] for i in 0:max_delay]


    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, vec(mean_X) # Return mean_X as a vector
end


"""
Constructs augmented state-space matrices A_full, B_full.
x_aug(k+1) = A_full * x_aug(k) + B_full * u(k)
x_aug = [x(k); x(k-1); ...; x(k-D); u(k-1); ...; u(k-D)]
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, delay::Int)
    println("Constructing augmented state-space model (A_full, B_full) for delay $delay...")
    D = delay # Number of delays for past states/inputs in augmented state

    # Size of augmented state vector
    # x(k), x(k-1)...x(k-D) -> N_X * (D+1) states
    # u(k-1)...u(k-D)     -> N_U * D states (u(k) is the input, not part of state)
    N_STATES_AUG = N_X * (D + 1) + N_U * D

    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    B_full = zeros(N_STATES_AUG, N_U)

    # Equation for x(k+1):
    # x(k+1) = A0*x(k) + A1*x(k-1) + ... + AD*x(k-D) +
    #          B0*u(k) + B1*u(k-1) + ... + BD*u(k-D)

    # Fill first N_X rows of A_full (dynamics of x(k+1))
    # Part from As (terms related to past x states)
    for i in 0:D # Corresponds to As[i+1] in 1-based Julia indexing
        col_start_idx = i * N_X + 1
        col_end_idx = (i + 1) * N_X
        A_full[1:N_X, col_start_idx:col_end_idx] = As[i+1] # As[1] is A0, As[D+1] is AD
    end

    # Part from Bs (terms related to past u states, u(k-1) to u(k-D))
    # These are u(k-1)...u(k-D) in the augmented state vector
    # Their part in A_full affects x(k+1)
    # Bs[1] is B0 (related to u(k) -> goes into B_full)
    # Bs[2] is B1 (related to u(k-1))
    # ...
    # Bs[D+1] is BD (related to u(k-D))
    if D > 0
        aug_state_u_offset = N_X * (D + 1) # Start index of u(k-1) in augmented state
        for i in 1:D # Corresponds to u(k-i) which uses Bs[i+1] (B_i)
            col_start_idx = aug_state_u_offset + (i - 1) * N_U + 1
            col_end_idx = aug_state_u_offset + i * N_U
            A_full[1:N_X, col_start_idx:col_end_idx] = Bs[i+1] # Bs[2] for u(k-1), ..., Bs[D+1] for u(k-D)
        end
    end

    # Fill B_full (effect of current input u(k))
    # x(k+1) gets B0*u(k)
    B_full[1:N_X, :] = Bs[1] # Bs[1] is B0

    # Fill identities for shifting states
    # x(k) at next step becomes x(k-1) at current step etc.
    if D > 0
        # Shift x(k) to x(k-1), x(k-1) to x(k-2), ..., x(k-D+1) to x(k-D)
        # Rows N_X+1 : N_X*(D+1) are for x(k-1)...x(k-D) at t+1
        # Columns 1 : N_X*D are for x(k)...x(k-D+1) at t
        for i in 0:(D-1) # for x(k-i) at t becoming x(k-i-1) at t+1
            dest_row_start = N_X * (i + 1) + 1
            dest_row_end = N_X * (i + 2)
            src_col_start = N_X * i + 1
            src_col_end = N_X * (i + 1)
            A_full[dest_row_start:dest_row_end, src_col_start:src_col_end] = I(N_X)
        end

        # Shift u(k) to u(k-1), u(k-1) to u(k-2), ..., u(k-D+1) to u(k-D)
        # u(k-1) at t+1 is u(k) at t. This is handled by B_full for u part of aug state.
        dest_u_block_start_row = N_X * (D + 1) + 1 # Start of u(k-1) block in x_aug(t+1)
        B_full[dest_u_block_start_row : dest_u_block_start_row + N_U - 1, :] = I(N_U)

        if D > 1
            # u(k-i) at t+1 is u(k-i+1) at t, for i=2..D
            # Rows for u(k-2)...u(k-D) at t+1
            # Columns for u(k-1)...u(k-D+1) at t
            src_u_block_start_col = N_X * (D + 1) + 1 # Start of u(k-1) block in x_aug(t)
            for i in 1:(D-1) # for u(k-i) at t becoming u(k-i-1) at t+1
                             # e.g. u(k-1) at t becomes u(k-2) at t+1
                dest_row_start = N_X * (D + 1) + i * N_U + 1
                dest_row_end = N_X * (D + 1) + (i + 1) * N_U

                src_col_start = N_X * (D + 1) + (i - 1) * N_U + 1
                src_col_end = N_X * (D + 1) + i * N_U
                A_full[dest_row_start:dest_row_end, src_col_start:src_col_end] = I(N_U)
            end
        end
    end
    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end


# --- Simulation and Initial State Functions ---

"""
Simulates the ABC model (augmented state-space).
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state_aug::AbstractVector, inputs_uk::AbstractMatrix, sim_length::Int)
    n_state_vars_aug = size(A_full, 1)
    _ , n_input_steps = size(inputs_uk) # inputs_uk is u(k) sequence

    if n_input_steps < sim_length
        error("Not enough input data provided for simulation ($sim_length steps needed, $n_input_steps provided).")
    end
    if length(initial_state_aug) != n_state_vars_aug
        error("Initial state dimension ($(length(initial_state_aug))) does not match A_full dimension ($n_state_vars_aug).")
    end


    # Store the history of the augmented state
    states_augmented_history = similar(initial_state_aug, n_state_vars_aug, sim_length)
    current_augmented_state = copy(initial_state_aug)

    for k in 1:sim_length
        states_augmented_history[:, k] = current_augmented_state
        current_input_uk = inputs_uk[:, k] # u(k)
        current_augmented_state = A_full * current_augmented_state + B_full * current_input_uk
    end
    return states_augmented_history # Returns the full augmented state history
end


"""
Creates the initial state vector for the ABC model.
Augmented State x_aug(t_start) = [x(t_start); x(t_start-1); ...; x(t_start-D); u(t_start-1); ...; u(t_start-D)]
All x components are centered.
"""
function create_abc_init_state(X_data_hist::AbstractMatrix, U_data_hist::AbstractMatrix,
                               sim_start_time_idx::Int, delay::Int, mean_x_vec::AbstractVector,
                               N_X::Int, N_U::Int)
    D = delay
    n_states_aug = N_X * (D + 1) + N_U * D
    init_state_aug = zeros(n_states_aug)

    # Check boundaries for historical data access
    # For x(t_start-D), need index sim_start_time_idx - D
    # For u(t_start-D), need index sim_start_time_idx - D
    min_required_data_idx = sim_start_time_idx - D
    if min_required_data_idx < 1
        error("Cannot create initial state at sim_start_time_idx $sim_start_time_idx with delay $D: requires data back to index $min_required_data_idx, which is < 1.")
    end
    if sim_start_time_idx > size(X_data_hist, 2) || sim_start_time_idx > size(U_data_hist, 2)
         error("sim_start_time_idx $sim_start_time_idx is out of bounds for X_data (max $(size(X_data_hist, 2))) or U_data (max $(size(U_data_hist, 2))).")
    end


    # Fill x(t_start), x(t_start-1), ..., x(t_start-D) - centered
    for i in 0:D # For x(t_start - i)
        data_idx = sim_start_time_idx - i
        state_block_start = i * N_X + 1
        state_block_end = (i + 1) * N_X
        init_state_aug[state_block_start:state_block_end] = X_data_hist[:, data_idx] .- mean_x_vec
    end

    # Fill u(t_start-1), ..., u(t_start-D) - original values
    if D > 0
        u_block_offset_in_aug_state = N_X * (D + 1) # Start index for u(t_start-1) part
        for i in 1:D # For u(t_start - i)
            data_idx = sim_start_time_idx - i
            input_block_start = u_block_offset_in_aug_state + (i - 1) * N_U + 1
            input_block_end = u_block_offset_in_aug_state + i * N_U
            init_state_aug[input_block_start:input_block_end] = U_data_hist[:, data_idx]
        end
    end
    return init_state_aug
end


# --- Configuration ---
const DATA_PATH = "./data.jld2" # Adjust path if needed
const START_IDX = 2000      # Preferred starting index for simulation comparison in original data
const SIM_LEN = 200         # Default length of the simulation comparison window
const PLOT_STATE_IDX = 1    # Which state component (from RELEVANT_STATES) to plot (1-based)
const RELEVANT_STATES = [7] # Relevant states indices from the raw X_data to be used
const NORMALIZE = true      # Whether to center the X data during DMDc

# --- Main Execution ---
function main_multi_delay_comparison()
    delays_to_compare = [0, 5, 20, 50, 100]
    println("Starting DMDc multi-delay comparison for delays: $delays_to_compare")

    # 1. Load and Preprocess Data (once)
    # N_X_loaded, N_U_loaded are dimensions *after* processing (e.g. RELEVANT_STATES, U preprocessing)
    X_data_loaded, U_data_loaded, N_X_loaded, N_U_loaded, TOTAL_SAMPLES_loaded = load_and_preprocess_data(DATA_PATH)

    if N_X_loaded == 0
        error("N_X_loaded is 0. Check RELEVANT_STATES and data loading.")
    end
     if N_U_loaded == 0 && size(U_data_loaded,1) > 0 # If U_data_loaded has rows but N_U_loaded is 0
        N_U_loaded = size(U_data_loaded,1) # Infer N_U from processed U_data_loaded
        println("Warning: N_U_loaded was 0, inferred as $(N_U_loaded) from U_data_loaded rows.")
    end
    if N_U_loaded == 0 # Still zero, this is an issue for DMDc with control
        println("Warning: N_U_loaded is 0 after attempting inference. DMDc might behave like DMD if B matrices are effectively zero.")
        # For safety, ensure Bs can be formed, even if it's zero dimensional inputs
        # This might mean calculate_dmdc_matrices needs to handle N_U=0 gracefully or error.
        # For now, we proceed, but Bs will be empty or cause errors if N_U truly is 0.
        # A robust way is to check N_U before calculate_dmdc_matrices if control is essential.
    end


    # 2. Setup Simulation Parameters for Comparison Plot
    # Determine a common simulation window [global_sim_start_time, global_sim_end_time]
    # Initial state requires history up to max_delay. So, sim_start_time >= max(delays_to_compare) + 1
    max_hist_delay_needed = isempty(delays_to_compare) ? 0 : maximum(delays_to_compare)

    # global_sim_start_time is the index in X_data_loaded/U_data_loaded
    global_sim_start_time = max(START_IDX, max_hist_delay_needed + 1)

    # Ensure sim_len is valid
    effective_sim_len = SIM_LEN
    if global_sim_start_time + effective_sim_len - 1 > TOTAL_SAMPLES_loaded
        println("Warning: Simulation window (start: $global_sim_start_time, len: $effective_sim_len) exceeds data bounds (total: $TOTAL_SAMPLES_loaded). Adjusting sim_len.")
        effective_sim_len = TOTAL_SAMPLES_loaded - global_sim_start_time + 1
    end

    if effective_sim_len <= 0
        error("Cannot determine a valid simulation length (effective_sim_len=$effective_sim_len). Check START_IDX, SIM_LEN, data size, and delay values. global_sim_start_time=$global_sim_start_time, TOTAL_SAMPLES_loaded=$TOTAL_SAMPLES_loaded")
    end
    println("Global simulation window: Start Index = $global_sim_start_time, Length = $effective_sim_len")

    # Extract actual data and inputs for this common window
    actual_data_to_plot_segment = X_data_loaded[:, global_sim_start_time : global_sim_start_time + effective_sim_len - 1]
    # PLOT_STATE_IDX is 1-based for X_data_loaded (which already incorporates RELEVANT_STATES)
    actual_data_trace = actual_data_to_plot_segment[PLOT_STATE_IDX, :]

    sim_inputs_global = U_data_loaded[:, global_sim_start_time : global_sim_start_time + effective_sim_len - 1]

    # 3. Initialize Plot
    plot_obj = plot(1:effective_sim_len, actual_data_trace,
                    label="Actual Data (X$(RELEVANT_STATES[PLOT_STATE_IDX]))", # Show original state index
                    linewidth=2.5, color=:black)

    # 4. Loop through delays, build model, simulate, and add to plot
    plot_colors = [:red, :blue, :green, :purple, :orange, :brown, :cyan]
    plot_styles = [:solid, :dash, :dot, :dashdot, :solid, :dash, :dot] # Cycle if more delays than styles

    for (idx, current_delay) in enumerate(delays_to_compare)
        println("\n--- Processing for DELAY = $current_delay ---")
        current_delay *= 2

        # Create DMDc model (As, Bs, M) and mean_X_vec for this delay
        # N_X_loaded, N_U_loaded are the dimensions of the data passed (X_data_loaded, U_data_loaded)
        As, Bs, M_dmdc, mean_X_vector = calculate_dmdc_matrices(X_data_loaded, U_data_loaded, current_delay, N_X_loaded, N_U_loaded, TOTAL_SAMPLES_loaded)

        # Construct augmented A_full, B_full matrices
        A_full, B_full = construct_augmented_matrices(As, Bs, N_X_loaded, N_U_loaded, current_delay)

        # Create Initial State for ABC model simulation
        # Uses global_sim_start_time as the reference point 'k' for the initial augmented state
        initial_aug_state = create_abc_init_state(X_data_loaded, U_data_loaded, global_sim_start_time, current_delay, mean_X_vector, N_X_loaded, N_U_loaded)

        # Run ABC Simulation
        println("Running ABC simulation for delay = $current_delay...")
        simulated_aug_states_history = simulate_abc_model(A_full, B_full, initial_aug_state, sim_inputs_global, effective_sim_len)

        # Extract the x(k) part (first N_X_loaded rows) and add mean back
        simulated_x_states = simulated_aug_states_history[1:N_X_loaded, :] .+ mean_X_vector

        # Get the specific trace to plot
        model_trace_to_plot = simulated_x_states[PLOT_STATE_IDX, :]

        # Add to plot
        current_color = plot_colors[mod1(idx, length(plot_colors))]
        current_style = plot_styles[mod1(idx, length(plot_styles))]
        curr_label = "Model (Nd=$(Int(current_delay/2)))"
        plot!(plot_obj, 1:effective_sim_len, model_trace_to_plot,
              label=curr_label,
              linewidth=1.5, color=current_color, linestyle=current_style)
        println("Finished simulation and plotting for delay = $current_delay.")
    end

    # 5. Finalize and Save Plot
    title!(plot_obj, "DMDc Model Comparison vs Actual Data")
    xlabel!(plot_obj, "Time (s) (from original data index $global_sim_start_time)")
    ylabel!(plot_obj, "Heatflux")
    plot!(plot_obj, legend=:topleft, legendfontsize=8) # Adjust legend position and size

    figname = "dmdc_multidelay_comp_state$(RELEVANT_STATES[PLOT_STATE_IDX])_start$(global_sim_start_time).svg"
    try
        savefig(plot_obj, figname)
        println("\nMulti-delay comparison plot saved to $figname")
    catch e
        println("\nError saving plot: $e")
        println("Make sure you have a compatible Plots.jl backend installed and working (e.g., GR, PyPlot).")
    end

    println("--- Script finished ---")
end

# --- Run the main function ---
main_multi_delay_comparison()

# Note: The original `simulate_orig_model` and its initial state creation
# `create_orig_init_state` are not used in this multi-delay comparison,
# which focuses on the ABC augmented model. They can be added if needed.
# The original `plot_comparison` function is also not used.
# The JLD2 saving of a single model from the old `main` is removed;
# if you need to save models for each delay, that logic would need to be
# added within the loop in `main_multi_delay_comparison`.
