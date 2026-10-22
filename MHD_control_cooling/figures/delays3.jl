#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays.
It constructs an augmented state-space model (A_full, B_full)
and compares its simulation results against the original DMDc formulation
and the actual data.

It also includes a function to calculate and plot the Autocorrelation Function (ACF)
for a selected state.

Key Functions:
- `load_and_preprocess_data`: Loads data from JLD2, preprocesses control input U.
- `calculate_and_plot_acf`: Calculates and plots the ACF for a given time series.
- `delay_embed`: Creates delay-embedded matrices.
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


# --- Configuration ---
const DATA_PATH = "./data.jld2" # Adjust path if needed
const DELAY = 580      # Number of delay steps for states and inputs, also max lag for ACF
const START_IDX = 2000   # Starting index for simulation comparison
const SIM_LEN = 700      # Default length of the simulation
const PLOT_STATE_IDX = 1 # Which state component to plot (and for ACF)
const RELEVANT_STATES = [6] # Relevant states for analysis
const NORMALIZE = false


# --- Core Data Processing and Model Creation Functions ---

"""
Loads X and U data, preprocesses U.
"""
function load_and_preprocess_data(data_path::String)
    println("Loading data from $data_path...")
    try
        X_data_raw = jldopen(data_path, "r")["X_data"]
        # Select relevant states and remove the first time step if it's an initial condition
        # The selection of [:, 2:end] depends on whether the first column is t=0 or t=1
        # and if it's used differently. Assuming it's consistent with original script's intent.
        X_data = X_data_raw[RELEVANT_STATES, 2:end]


        U_data_raw = jldopen(data_path, "r")["U_data"]

        println("Preprocessing U data...")
        U = let
            cols = []
            for u_col in eachcol(U_data_raw)
                if size(u_col, 1) == 8
                    push!(cols, vec(u_col[1:4] * u_col[5:8]'))
                else
                    println("Warning: Unexpected size for column in U_data_raw. Skipping.")
                end
            end
            hcat(cols...)
        end

        N_X = size(X_data, 1) # Number of features in X_data (after selecting RELEVANT_STATES)
        N_U = size(U, 1)
        TOTAL_SAMPLES = size(X_data, 2) # Number of samples in X_data (after column selection)
        println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")

        return X_data, U, N_X, N_U, TOTAL_SAMPLES
    catch e
        println("Error loading data from $data_path: $e")
        println("Please ensure the JLD2 file exists, the path is correct, and contains 'X_data' and 'U_data'.")
        rethrow(e) # Propagate the error
    end
end

"""
Calculates and plots the Autocorrelation Function (ACF) for a given time series.
"""
function calculate_and_plot_acf(time_series::AbstractVector, max_lag::Int, state_label::String, output_dir::String=".")
    println("Calculating ACF for state heatflux up to τ=$max_lag...")
    n = length(time_series)
    if n <= max_lag
        println("Warning: Time series length ($n) is not sufficiently greater than max_lag ($max_lag) for ACF calculation. Skipping ACF plot.")
        return
    end

    mean_ts = mean(time_series)
    var_ts = var(time_series, corrected=false) # Use uncorrected variance for ACF definition

    if var_ts == 0
        println("Warning: Variance of the time series is zero. Cannot compute ACF. Skipping ACF plot.")
        return
    end

    acf_values = zeros(Float64, max_lag + 1)

    for lag in 0:max_lag
        autocov_sum = 0.0
        for i in 1:(n-lag)
            autocov_sum += (time_series[i] - mean_ts) * (time_series[i+lag] - mean_ts)
        end
        acf_values[lag+1] = (autocov_sum / n) / var_ts # Some definitions use (n-lag) in autocov sum, others use n. Using n for consistency.
        # Normalizing by total variance.
    end
    # ACF at lag 0 should be 1
    acf_values[1] = 1.0


    println("Plotting ACF...")
    lags = 0:max_lag
    plot(lags, acf_values, seriestype=:bar,
        title="ACF for heatflux up to τ=$max_lag",
        xlabel="Lag (\$\\tau\$)",
        ylabel="Autocorrelation",
        legend=false,
        xticks=0:round(Int, max_lag / 10):max_lag, # Adjust tick frequency as needed
        ylims=(-1, 1.1)) # Typical ACF range

    # Add confidence bounds (e.g., +/- 1.96/sqrt(n))
    # These are approximate for large n and if the series were white noise
    # conf_bound = 1.96 / sqrt(n)
    # hline!([conf_bound, -conf_bound], linestyle=:dash, color=:red, label="95% Conf. Bounds")


    figname = joinpath(output_dir, "acf_plot_state_$(state_label)_delay$(max_lag).png")
    try
        savefig(figname)
        println("ACF plot saved to $figname")
    catch e
        println("Error saving ACF plot: $e")
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
        return data
    end

    if n_samples <= delay
        error("Not enough samples ($n_samples) to create embedding with delay ($delay).")
    end

    n_embedded_samples = n_samples - delay
    embedded_dim = n_features * (delay + 1)
    # Initialize with zeros or appropriate type
    embedded_data = zeros(eltype(data), embedded_dim, n_embedded_samples)


    for t in 1:n_embedded_samples
        for d in 0:delay
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features
            source_col = t + delay - d # This correctly maps to x(t), x(t-1)...x(t-delay) for column t of embedded_data
            # if embedded_data[:,t] should represent [x_t; x_{t-1}; ...; x_{t-delay}]
            # The original script's `source_col = t + delay - d` implies that
            # embedded_data[:,t] contains [data[:, t+delay]; data[:,t+delay-1]; ...; data[:,t]]
            # Let's stick to the original script's logic for `delay_embed`
            embedded_data[row_start:row_end, t] = data[:, source_col]
        end
    end
    return embedded_data
end

function input_lift(U::AbstractMatrix)
    U_lifted = U # Placeholder for potential future lifting functions
    return U_lifted
end

"""
Calculates DMDc matrices (As, Bs, M) from data.
"""
function calculate_dmdc_matrices(X_data::AbstractMatrix, U::AbstractMatrix, delay::Int, N_X::Int, N_U::Int, TOTAL_SAMPLES::Int)
    println("Creating DMDc model with delay $delay...")
    max_delay = delay # This is D in many notations

    # Prepare data indices for DMDc
    # X_current (X') will be X_data[:, max_delay+2 : TOTAL_SAMPLES]
    # X_delayed (X) will be delay_embed(X_data[:, 1 : TOTAL_SAMPLES-1], max_delay)
    # U_delayed (Gamma) will be delay_embed(U_data[:, 1 : TOTAL_SAMPLES-1], max_delay)

    # Indices for X_data to form X_current (X_k+1)
    # X_current needs to align with the output of the delay embedding
    # If X_delayed ends at sample (TOTAL_SAMPLES-1-max_delay) + max_delay = TOTAL_SAMPLES-1 (original time index)
    # then X_current should start from the next time step.
    # The number of columns in delay_embed(data, delay) is size(data,2) - delay.
    # So, X_delayed has (TOTAL_SAMPLES-1) - max_delay columns.
    # The time steps covered by X_delayed[:,j] are [ (j+max_delay)-th original sample, ..., j-th original sample ]
    # The last column of X_delayed, X_delayed[:, end], uses original samples up to TOTAL_SAMPLES-1.
    # X_current should be the state *after* the latest state in the corresponding X_delayed column.
    # The latest state in X_delayed[:,j] is data[:, j+max_delay]. So X_current[:,j] should be data[:, j+max_delay+1].

    # Number of samples available for the "history" part of DMDc
    num_hist_samples = TOTAL_SAMPLES - 1
    if num_hist_samples <= max_delay
        error("Not enough samples in X_data (after initial processing) to form X_delayed for DMDc with delay $max_delay. Have $num_hist_samples, need > $max_delay.")
    end

    X_for_delay_embed = X_data[:, 1:num_hist_samples]
    U_for_delay_embed = U[:, 1:num_hist_samples]

    X_delayed = delay_embed(X_for_delay_embed, max_delay) # Represents [x_k, x_{k-1}, ..., x_{k-D}]
    U_delayed = delay_embed(U_for_delay_embed, max_delay) # Represents [u_k, u_{k-1}, ..., u_{k-D}]

    # X_current (X_k+1)
    # The first column of X_delayed corresponds to original time index max_delay+1 (using 1-based indexing for original data before slicing for DMDc)
    # So, X_current should start from original time index max_delay+2
    # The number of columns in X_delayed is size(X_for_delay_embed, 2) - max_delay
    num_cols_embedded = size(X_for_delay_embed, 2) - max_delay

    # X_current should be X_data[:, (max_delay+1)+1 : (max_delay+1) + num_cols_embedded ]
    # which is X_data[:, max_delay+2 : max_delay+1 + num_cols_embedded]
    # Example: TOTAL_SAMPLES (in X_data) = 100, max_delay = 5
    # X_for_delay_embed = X_data[:, 1:99]
    # X_delayed = delay_embed(X_data[:, 1:99], 5) -> 99-5 = 94 columns
    # First column of X_delayed uses X_data[:, 1] to X_data[:, 6]
    # Last column of X_delayed uses X_data[:, 94] to X_data[:, 99]
    # X_current should be X_data[:, 7 : 94+6] = X_data[:, 7:100]
    # So, X_current starts at index max_delay + 1 + 1 (from X_data)
    # and has num_cols_embedded columns.
    start_idx_X_current = max_delay + 1 + 1 # This is the starting column index in the *original X_data before any slicing*
    # But X_data here is already X_data_raw[RELEVANT_STATES, 2:end]
    # So, if X_data has TOTAL_SAMPLES columns, indices are 1..TOTAL_SAMPLES
    start_idx_X_current_local = max_delay + 1 # This is the index for X_current in the *current* X_data
    # e.g. if X_delayed[:,1] uses X_data[:,1]...X_data[:,max_delay+1]
    # then X_current[:,1] should be X_data[:,max_delay+2]
    # Let's re-evaluate based on typical DMDc formulation:
    # Omega = [X_k_embed; U_k_embed]
    # X_k+1_target = X_{k+1}
    # X_k_embed contains [x_k, x_{k-1}, ..., x_{k-D}]
    # U_k_embed contains [u_k, u_{k-1}, ..., u_{k-D}]
    # The time index 'k' for the embedding runs from D+1 to N-1 (if N is total original samples)
    # X_{k+1}_target then runs from D+2 to N

    # Let's use the script's original indexing logic if it was working.
    # The original script had:
    # start_idx_curr = max_delay + 2
    # end_idx_curr = TOTAL_SAMPLES
    # start_idx_hist = 1
    # end_idx_hist = TOTAL_SAMPLES - 1
    # X_current = X_data[:, start_idx_curr:end_idx_curr]
    # X_delayed = delay_embed(X_data[:, start_idx_hist:end_idx_hist], max_delay)
    # U_delayed = delay_embed(U[:, start_idx_hist:end_idx_hist], max_delay)

    # This implies X_data and U are 1-indexed from the start of the relevant segment.
    # TOTAL_SAMPLES is size(X_data, 2) which is after X_data_raw[:,2:end]

    _X_hist = X_data[:, 1:(TOTAL_SAMPLES-1)]
    _U_hist = U[:, 1:(TOTAL_SAMPLES-1)]

    if size(_X_hist, 2) <= max_delay || size(_U_hist, 2) <= max_delay
        error("Not enough samples for history embedding. X_hist size: $(size(_X_hist,2)), U_hist size: $(size(_U_hist,2)), max_delay: $max_delay")
    end

    X_embed_dmd = delay_embed(_X_hist, max_delay) # Columns: (TOTAL_SAMPLES-1) - max_delay
    U_embed_dmd = delay_embed(_U_hist, max_delay) # Columns: (TOTAL_SAMPLES-1) - max_delay

    # X_next_dmd should be X_data shifted.
    # If X_embed_dmd[:,i] uses _X_hist[:, i]..._X_hist[:,i+max_delay] (from delay_embed logic)
    # then X_next_dmd[:,i] should be X_data[:, i+max_delay+1]
    num_embedded_cols = size(X_embed_dmd, 2)
    X_next_dmd = X_data[:, (max_delay+1)+1:(max_delay+1)+num_embedded_cols]


    # Align dimensions (already done by careful slicing)
    common_len = size(X_next_dmd, 2)
    if common_len != size(X_embed_dmd, 2) || common_len != size(U_embed_dmd, 2)
        error("Dimension mismatch after slicing for DMDc: X_next $(size(X_next_dmd,2)), X_embed $(size(X_embed_dmd,2)), U_embed $(size(U_embed_dmd,2))")
    end
    if common_len <= 0
        error("Resulting common length for DMDc matrices is non-positive ($common_len).")
    end

    #==== LIFITING HERE ====#
    U_embed_dmd = abs.(U_embed_dmd) # As per original script

    # Center X data
    local mean_X_val # Use a different name to avoid conflict with function argument
    if NORMALIZE
        # Calculate mean from the portion of X_data that corresponds to X_embed_dmd states
        # This is complex. Original script used mean(X_data, dims=2) which is mean over ALL preprocessed X_data.
        mean_X_val = mean(X_data, dims=2)
    else
        mean_X_val = zeros(N_X, 1)
    end
    X_next_centered = X_next_dmd .- mean_X_val

    # Centering X_embed_dmd: each block of N_X rows corresponds to a time snapshot
    X_embed_centered = similar(X_embed_dmd)
    for i in 0:max_delay
        row_slice = (i*N_X+1):((i+1)*N_X)
        X_embed_centered[row_slice, :] = X_embed_dmd[row_slice, :] .- mean_X_val
    end

    # Combined state for regression
    Omega = [X_embed_centered; U_embed_dmd]

    # Solve for the combined dynamics matrix M
    println("Performing DMDc regression...")
    M = X_next_centered * pinv(Omega)
    println("DMDc matrix M size: $(size(M))")

    # Extract A and B matrices
    A_matrices_flat = M[:, 1:N_X*(max_delay+1)]
    B_matrices_flat = M[:, N_X*(max_delay+1)+1:end]

    # Reshape into list of As and Bs
    # As[1] is A0, As[2] is A1 ... As[delay+1] is AD
    # Bs[1] is B0, Bs[2] is B1 ... Bs[delay+1] is BD
    As = [A_matrices_flat[:, i*N_X+1:(i+1)*N_X] for i in 0:max_delay]
    Bs = [B_matrices_flat[:, i*N_U+1:(i+1)*N_U] for i in 0:max_delay]


    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, vec(mean_X_val) # Return the actual mean used
end


"""
Constructs augmented state-space matrices A_full, B_full.
A_full maps [x_k; x_{k-1};...;x_{k-D}; u_{k-1};...;u_{k-D}] to next state
B_full maps u_k to next state
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, delay::Int)
    println("Constructing augmented state-space model (A_full, B_full)...")
    D = delay # Number of delays (e.g., D=1 means x_k, x_{k-1})

    # Augmented state: [x_k; x_{k-1}; ...; x_{k-D}; u_{k-1}; ...; u_{k-D}]
    # Size of x part: N_X * (D+1)
    # Size of u part: N_U * D
    N_STATES_AUG = N_X * (D + 1) + N_U * D

    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    B_full = zeros(N_STATES_AUG, N_U)

    # Equation for x_{k+1}:
    # x_{k+1} = A0*x_k + A1*x_{k-1} + ... + AD*x_{k-D} + B0*u_k + B1*u_{k-1} + ... + BD*u_{k-D}

    # Row 1: Dynamics of x_k -> x_{k+1}
    # Part affecting x_{k+1} from past x states
    for i in 0:D # Corresponds to As[1] (A0) to As[D+1] (AD)
        A_full[1:N_X, (i*N_X+1):((i+1)*N_X)] = As[i+1]
    end
    # Part affecting x_{k+1} from past u states (u_{k-1} to u_{k-D})
    # Bs[2] (B1) relates to u_{k-1}, Bs[D+1] (BD) relates to u_{k-D}
    # Start index for u_{k-1} in augmented state: N_X*(D+1) + 1
    if D > 0
        for i in 1:D # Corresponds to Bs[2] (B1) to Bs[D+1] (BD)
            u_col_offset_in_A = N_X * (D + 1) + (i - 1) * N_U
            A_full[1:N_X, (u_col_offset_in_A+1):(u_col_offset_in_A+N_U)] = Bs[i+1]
        end
    end
    # Part affecting x_{k+1} from current u_k (this goes into B_full)
    B_full[1:N_X, :] = Bs[1] # B0

    # Row 2: Shifting x_k to x_{k-1} position for next step, etc.
    # x_{k} becomes x_{(k+1)-1}
    # x_{k-1} becomes x_{(k+1)-2}
    # ...
    # x_{k-D+1} becomes x_{(k+1)-D}
    if D > 0
        for i in 0:(D-1) # From x_k down to x_{k-D+1}
            # Current block: x_{k-i} located at rows i*N_X+1 : (i+1)*N_X
            # Next step, this block moves to x_{(k+1)-(i+1)}
            # which is located at rows (i+1)*N_X+1 : (i+2)*N_X
            dest_row_start = (i + 1) * N_X + 1
            dest_row_end = (i + 2) * N_X

            src_col_start = i * N_X + 1
            src_col_end = (i + 1) * N_X

            A_full[dest_row_start:dest_row_end, src_col_start:src_col_end] = I(N_X)
        end
    end

    # Row 3: Shifting u_{k-1} to u_{k-2} position, etc. And u_k to u_{k-1}
    # u_k (from B_full) becomes u_{(k+1)-1}
    # u_{k-1} becomes u_{(k+1)-2}
    # ...
    # u_{k-D+1} becomes u_{(k+1)-D}

    # u_k input goes into B_full and then shifts to first delayed u slot
    u_target_first_delay_slot_row_start = N_X * (D + 1) + 1
    if D > 0
        B_full[u_target_first_delay_slot_row_start:u_target_first_delay_slot_row_start+N_U-1, :] = I(N_U)
    end

    # Shifting existing delayed inputs
    if D > 1 # Only if there's u_{k-2} or further delays
        for i in 1:(D-1) # From u_{k-i} to u_{k-(i+1)}
            # Current block u_{k-i} located at N_X*(D+1) + (i-1)*N_U + 1
            src_col_start_u = N_X * (D + 1) + (i - 1) * N_U + 1
            src_col_end_u = src_col_start_u + N_U - 1

            # Next step, this block moves to u_{(k+1)-(i+1)}
            # which is located at N_X*(D+1) + i*N_U + 1
            dest_row_start_u = N_X * (D + 1) + i * N_U + 1
            dest_row_end_u = dest_row_start_u + N_U - 1

            A_full[dest_row_start_u:dest_row_end_u, src_col_start_u:src_col_end_u] = I(N_U)
        end
    end

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end


"""
Main function to create the ABC system from data.
"""
function create_abc_system_from_data(data_path::String, delay::Int)
    X_data, U, N_X, N_U, TOTAL_SAMPLES = load_and_preprocess_data(data_path)

    # Note: N_X here is the number of rows in the *processed* X_data
    # (i.e., after selecting RELEVANT_STATES).
    # PLOT_STATE_IDX refers to a row in this processed X_data.

    # Call ACF plotting function
    # Ensure PLOT_STATE_IDX is valid for the processed X_data
    if 1 <= PLOT_STATE_IDX <= N_X
        # The label for the state in ACF plot should reflect its original index if RELEVANT_STATES is used
        # If RELEVANT_STATES = [6] and PLOT_STATE_IDX = 1, the original state index is 6.
        original_state_label = string(RELEVANT_STATES[PLOT_STATE_IDX])
        calculate_and_plot_acf(X_data[PLOT_STATE_IDX, :], DELAY, original_state_label) # DELAY is max_acf_lag
    else
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is out of bounds for processed X_data with $N_X states. Skipping ACF plot.")
    end

    As, Bs, M, mean_X = calculate_dmdc_matrices(X_data, U, delay, N_X, N_U, TOTAL_SAMPLES)
    A_full, B_full = construct_augmented_matrices(As, Bs, N_X, N_U, delay)

    return A_full, B_full, As, Bs, M, vec(mean_X), N_X, N_U, X_data, U, TOTAL_SAMPLES
end

# --- Simulation and Reshape Functions ---

"""
Simulates the ABC model (augmented state-space).
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state::AbstractVector, inputs::AbstractMatrix, sim_length::Int)
    n_state_vars_aug = size(A_full, 1)
    n_inputs_u, n_input_steps = size(inputs)

    if n_input_steps < sim_length
        error("Not enough input data provided for simulation ($sim_length steps needed, $n_input_steps provided).")
    end
    if size(B_full, 2) != n_inputs_u
        error("B_full input dimension ($(size(B_full,2))) does not match input data dimension ($n_inputs_u).")
    end


    states_augmented = zeros(eltype(initial_state), n_state_vars_aug, sim_length)
    current_state_aug = copy(initial_state)
    abc_fn(state, input) = A_full * state + B_full * input

    for k in 1:sim_length
        states_augmented[:, k] = current_state_aug
        current_state_aug = abc_fn(current_state_aug, inputs[:, k])
    end
    return states_augmented
end

"""
Simulates the original DMDc model formulation using the M matrix.
Omega_k = [X_embed_centered_k; U_embed_k]
X_next_centered = M * Omega_k
"""
function simulate_orig_model(initial_state_orig_model_centered::AbstractVector, # [x_embed_centered; u_embed_delayed]
    inputs::AbstractMatrix, # u_k, u_{k+1}, ...
    M_dmdc::AbstractMatrix,
    sim_length::Int,
    delay::Int, # D
    n_x::Int,
    n_u::Int)

    states_out_centered = zeros(n_x, sim_length) # Stores x_{k+1}_centered

    # initial_state_orig_model_centered should contain:
    # [x_k_centered; x_{k-1}_centered; ...; x_{k-D}_centered;  (N_X*(D+1))
    #  u_{k-1}; u_{k-2}; ...; u_{k-D}]                     (N_U*D)
    # The M_dmdc matrix expects Omega = [X_embed_centered; U_embed], where U_embed is [u_k; u_{k-1}; ...; u_{k-D}]

    current_x_embed_centered = initial_state_orig_model_centered[1:n_x*(delay+1)]
    current_u_delayed_embed = delay > 0 ? initial_state_orig_model_centered[n_x*(delay+1)+1:end] : zeros(0)


    for k in 1:sim_length
        current_input_uk = inputs[:, k] # This is u_k for the current step

        # Construct U_embed_k = [u_k; u_{k-1}; ...; u_{k-D}]
        u_embed_k = zeros(n_u * (delay + 1))
        u_embed_k[1:n_u] = current_input_uk
        if delay > 0
            u_embed_k[n_u+1:end] = current_u_delayed_embed
        end

        # Lift U_embed_k if needed (abs was applied during training)
        u_embed_k_lifted = abs.(u_embed_k)


        omega_k = [current_x_embed_centered; u_embed_k_lifted]
        x_next_centered = M_dmdc * omega_k
        states_out_centered[:, k] = x_next_centered

        # Update current_x_embed_centered for the next iteration:
        # [x_{next}_centered; x_k_centered; ...; x_{k-D+1}_centered]
        if delay > 0
            current_x_embed_centered = [x_next_centered; current_x_embed_centered[1:n_x*delay]]
        else
            current_x_embed_centered = x_next_centered
        end

        # Update current_u_delayed_embed for the next iteration:
        # [u_k; u_{k-1}; ...; u_{k-D+1}]
        if delay > 0
            if delay == 1
                current_u_delayed_embed = current_input_uk # u_k becomes the new u_{k-1}
            else # delay > 1
                current_u_delayed_embed = [current_input_uk; current_u_delayed_embed[1:n_u*(delay-1)]]
            end
        end
    end
    return states_out_centered
end


"""
Creates the initial state vector for the ABC model.
Augmented state: [x_k_centered; x_{k-1}_centered; ...; x_{k-D}_centered; u_{k-1}; ...; u_{k-D}]
All x components are centered. U components are not.
"""
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time_idx::Int, delay::Int, mean_x_vec::AbstractVector, N_X::Int, N_U::Int)
    D = delay
    n_states_aug = N_X * (D + 1) + N_U * D
    init_state_aug = zeros(n_states_aug)

    # Check boundaries for X and U based on start_time_idx and delay D
    # We need X from start_time_idx down to start_time_idx - D
    # We need U from start_time_idx - 1 down to start_time_idx - D
    min_required_X_idx = start_time_idx - D
    min_required_U_idx = start_time_idx - D # For u_{k-D}, if k is start_time_idx

    if min_required_X_idx < 1
        error("Cannot create initial state for ABC: start_time_idx $start_time_idx with delay $D requires X data back to index $min_required_X_idx, which is < 1.")
    end
    if D > 0 && min_required_U_idx < 1
        error("Cannot create initial state for ABC: start_time_idx $start_time_idx with delay $D requires U data back to index $min_required_U_idx, which is < 1.")
    end


    # Fill x part: [x_k_centered; x_{k-1}_centered; ...; x_{k-D}_centered]
    for i in 0:D # i=0 is x_k, i=1 is x_{k-1}, ..., i=D is x_{k-D}
        current_x_original_time_idx = start_time_idx - i
        dest_idx_start = i * N_X + 1
        dest_idx_end = (i + 1) * N_X
        init_state_aug[dest_idx_start:dest_idx_end] = X[:, current_x_original_time_idx] .- mean_x_vec
    end

    # Fill u part: [u_{k-1}; u_{k-2}; ...; u_{k-D}]
    if D > 0
        u_part_start_in_aug_state = N_X * (D + 1) + 1
        for i in 1:D # i=1 is u_{k-1}, ..., i=D is u_{k-D}
            current_u_original_time_idx = start_time_idx - i
            dest_idx_start = u_part_start_in_aug_state + (i - 1) * N_U
            dest_idx_end = dest_idx_start + N_U - 1
            init_state_aug[dest_idx_start:dest_idx_end] = U[:, current_u_original_time_idx]
        end
    end
    return init_state_aug
end


"""
Creates the initial state vector for the original DMDc model simulation.
State: [x_k_centered; ... x_{k-D}_centered; u_{k-1}; ... u_{k-D}]
This is the initial Omega_k, but without the u_k part (which comes from inputs).
The state vector for simulation loop should be [current_x_embed_centered; current_u_delayed_embed]
"""
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time_idx::Int, delay::Int, mean_x_vec::AbstractVector, N_X::Int, N_U::Int)
    D = delay
    # Length of x_embed_centered part + u_delayed_embed part
    init_state_len = N_X * (D + 1) + N_U * D
    init_state_orig = zeros(init_state_len)

    # Create the same structure as create_abc_init_state, as this is what the simulation loop expects to deconstruct
    # Fill x part: [x_k_centered; x_{k-1}_centered; ...; x_{k-D}_centered]
    min_required_X_idx = start_time_idx - D
    if min_required_X_idx < 1
        error("Cannot create initial state for Orig: start_time_idx $start_time_idx with delay $D requires X data back to index $min_required_X_idx, which is < 1.")
    end

    for i in 0:D # i=0 is x_k, i=1 is x_{k-1}, ..., i=D is x_{k-D}
        current_x_original_time_idx = start_time_idx - i
        dest_idx_start = i * N_X + 1
        dest_idx_end = (i + 1) * N_X
        init_state_orig[dest_idx_start:dest_idx_end] = X[:, current_x_original_time_idx] .- mean_x_vec
    end

    # Fill u part (delayed inputs): [u_{k-1}; u_{k-2}; ...; u_{k-D}]
    if D > 0
        min_required_U_idx = start_time_idx - D
        if min_required_U_idx < 1
            error("Cannot create initial state for Orig: start_time_idx $start_time_idx with delay $D requires U data back to index $min_required_U_idx, which is < 1.")
        end
        u_part_start_in_state = N_X * (D + 1) + 1
        for i in 1:D # i=1 is u_{k-1}, ..., i=D is u_{k-D}
            current_u_original_time_idx = start_time_idx - i
            dest_idx_start = u_part_start_in_state + (i - 1) * N_U
            dest_idx_end = dest_idx_start + N_U - 1
            init_state_orig[dest_idx_start:dest_idx_end] = U[:, current_u_original_time_idx]
        end
    end
    return init_state_orig
end


"""
Plots the comparison between actual data and model simulations.
"""
function plot_comparison(actual_data_segment::AbstractMatrix, states_abc::AbstractMatrix, states_orig::AbstractMatrix,
    plot_state_row_idx::Int, sim_start_abs_time::Int, sim_len::Int, delay_val::Int, output_dir::String=".")

    println("Plotting results...")
    # actual_data_segment is already the segment X_data[:, sim_start_time:sim_start_time+sim_len-1]
    # states_abc and states_orig are the simulated outputs for these steps, also sim_len columns.
    # plot_state_row_idx is the row index within these matrices (e.g., 1 if PLOT_STATE_IDX=1 and RELEVANT_STATES selected one state)

    time_axis = 1:sim_len

    # Determine the original state index for labeling
    # If RELEVANT_STATES = [6] and plot_state_row_idx = 1, original label is "X6"
    original_state_label = "X?"
    if 1 <= plot_state_row_idx <= length(RELEVANT_STATES)
        original_state_label = "X$(RELEVANT_STATES[plot_state_row_idx])"
    else # Should not happen if PLOT_STATE_IDX is validated
        original_state_label = "X (Row $plot_state_row_idx)"
    end


    plot(time_axis, actual_data_segment[plot_state_row_idx, :],
        label="Actual Data ($original_state_label)",
        linewidth=2, color=:black)
    plot!(time_axis, states_abc[plot_state_row_idx, :],
        label="ABC Model ($original_state_label)",
        linewidth=1.5, linestyle=:dash, color=:red)
    plot!(time_axis, states_orig[plot_state_row_idx, :],
        label="Original DMDc Model ($original_state_label)",
        linewidth=1.5, linestyle=:dot, color=:blue)

    title!("Model Comparison (Delay = $delay_val)")
    xlabel!("Time Step (relative to original data index $sim_start_abs_time)")
    ylabel!("State Value")
    plot!(legend=:outertopright)

    figname = joinpath(output_dir, "dmdc_comparison_delay$(delay_val)_start$(sim_start_abs_time).png")
    try
        savefig(figname)
        println("Plot saved to $figname")
    catch e
        println("Error saving comparison plot: $e")
    end
end


# --- Main Execution ---
function main()
    println("Starting main execution...")
    # 1. Create ABC System (this now also plots ACF)
    # N_X from create_abc_system_from_data is the number of rows in the *processed* X_data
    A_full, B_full, As, Bs, M_dmdc, mean_X_vec, N_X_processed, N_U_processed, X_data_processed, U_processed, TOTAL_SAMPLES_processed =
        create_abc_system_from_data(DATA_PATH, DELAY)

    # 2. Simulation Setup
    # sim_start_time is the column index in X_data_processed / U_processed
    sim_start_time = START_IDX
    # Ensure sim_start_time is valid considering the delay needed for initial state
    # The initial state creation needs data from (sim_start_time - DELAY)
    # So, sim_start_time must be >= DELAY + 1 (for 1-based indexing)
    if sim_start_time < DELAY + 1
        println("Warning: START_IDX ($START_IDX) is too early for DELAY ($DELAY). Adjusting simulation start time.")
        sim_start_time = DELAY + 1
    end

    actual_sim_len = SIM_LEN
    if sim_start_time + actual_sim_len - 1 > TOTAL_SAMPLES_processed
        println("Warning: Simulation length exceeds data bounds. Adjusting sim_len.")
        actual_sim_len = TOTAL_SAMPLES_processed - sim_start_time + 1
    end

    if actual_sim_len <= 0
        println("Error: Calculated simulation length is non-positive ($actual_sim_len). Cannot simulate. Check START_IDX and data length.")
        return
    end

    sim_inputs = U_processed[:, sim_start_time:sim_start_time+actual_sim_len-1]

    # 3. Create Initial States
    println("Creating initial states for simulation starting at processed data index $sim_start_time...")
    init_state_abc = create_abc_init_state(X_data_processed, U_processed, sim_start_time, DELAY, mean_X_vec, N_X_processed, N_U_processed)
    init_state_orig = create_orig_init_state(X_data_processed, U_processed, sim_start_time, DELAY, mean_X_vec, N_X_processed, N_U_processed)

    # 4. Run Simulations
    println("Running simulations for $actual_sim_len steps...")
    states_abc_full_augmented = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, actual_sim_len)
    # Extract the x-part and add mean back
    states_abc_simulated = states_abc_full_augmented[1:N_X_processed, :] .+ mean_X_vec

    states_orig_centered_simulated = simulate_orig_model(init_state_orig, sim_inputs, M_dmdc, actual_sim_len, DELAY, N_X_processed, N_U_processed)
    states_orig_simulated = states_orig_centered_simulated .+ mean_X_vec

    # 5. Get Actual Data for Comparison
    actual_data_for_comparison = X_data_processed[:, sim_start_time:sim_start_time+actual_sim_len-1]
    println("Simulations complete.")

    # 6. Plot Results
    # PLOT_STATE_IDX is the 1-based row index into X_data_processed (after RELEVANT_STATES selection)
    if 1 <= PLOT_STATE_IDX <= N_X_processed
        plot_comparison(actual_data_for_comparison, states_abc_simulated, states_orig_simulated,
            PLOT_STATE_IDX, sim_start_time, actual_sim_len, DELAY)
    else
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is invalid for plotting $N_X_processed processed states. Skipping comparison plot.")
    end


    println("Saving model data to ABC.jld2...")
    try
        jldsave("ABC.jld2"; # Use semicolon for keyword arguments in jldsave
            A=A_full,
            B=B_full,
            delay=DELAY,
            mean_vec=mean_X_vec,
            M_dmdc=M_dmdc,
            N_X=N_X_processed, # Number of states after RELEVANT_STATES processing
            N_U=N_U_processed,
            TOTAL_SAMPLES_PROCESSED=TOTAL_SAMPLES_processed, # Samples in X_data_processed
            RELEVANT_STATES_INDICES=RELEVANT_STATES # Original indices that were selected
        )
        println("Model data saved.")
    catch e
        println("Error saving JLD2 file: $e")
    end

    println("Script finished.")
end

# --- Entry Point ---
if abspath(PROGRAM_FILE) == @__FILE__
    # Check if data file exists
    if !isfile(DATA_PATH)
        println("Error: Data file not found at $DATA_PATH")
        println("Please create a dummy 'data.jld2' file or provide the correct path.")
        # Create a dummy data file for testing if it doesn't exist
        println("Creating a dummy data.jld2 for testing purposes...")
        try
            num_original_states = 10
            num_original_inputs_raw = 8 # before processing to 16
            num_total_raw_samples = START_IDX + SIM_LEN + DELAY + 100 # ensure enough data

            _X_data_raw = rand(num_original_states, num_total_raw_samples)
            _U_data_raw = rand(num_original_inputs_raw, num_total_raw_samples)
            jldsave(DATA_PATH; X_data=_X_data_raw, U_data=_U_data_raw)
            println("Dummy data.jld2 created with X_data: $(size(_X_data_raw)), U_data: $(size(_U_data_raw)).")
            println("You might need to adjust RELEVANT_STATES if it's out of bounds for $num_original_states states.")
            # Adjust RELEVANT_STATES if dummy data is created and current RELEVANT_STATES is out of bounds
            if maximum(RELEVANT_STATES) > num_original_states
                println("Warning: RELEVANT_STATES ($RELEVANT_STATES) is out of bounds for dummy data ($num_original_states states).")
                # RELEVANT_STATES = [min(num_original_states,1)] # Reset to a safe value
                # println("RELEVANT_STATES has been reset to $RELEVANT_STATES for the dummy data.")
                # This modification would require RELEVANT_STATES to not be const.
                # For now, user needs to ensure RELEVANT_STATES is valid for dummy data if used.
            end
        catch e_dummy
            println("Failed to create dummy data file: $e_dummy")
        end
    end
    main()
end
