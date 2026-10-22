#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays.
The DMDc model is trained on delay embeddings constructed from the
original resolution data, where embeddings use a specified lag.
The comparison plot shows the original non-strided actual data against
the predictions from the DMDc model.

It also includes functions to calculate and plot:
- Autocorrelation Function (ACF)
- Average Mutual Information (AMI)
for a selected state from the original processed data.

Key Functions:
- `load_and_preprocess_data`: Loads data, preprocesses U. Returns original processed X and U.
- `calculate_and_plot_acf`: Calculates and plots the ACF on original processed data.
- `calculate_ami`: Computes AMI between two series using histograms.
- `calculate_and_plot_ami`: Calculates and plots AMI vs. lag on original processed data.
- `delay_embed`: Creates delay-embedded matrices from original data using specified lags.
- `calculate_dmdc_matrices`: Computes DMDc system matrices from full embedded data.
- `construct_augmented_matrices`: Builds the augmented A_full, B_full matrices.
- `create_abc_system_from_data`: Orchestrates data loading, analysis, DMDc, and ABC matrix construction.
- `create_abc_init_state`, `create_orig_init_state`: Generate initial states for simulations.
- `simulate_abc_model`, `simulate_orig_model`: Run simulations (predicting original time steps).
- `plot_comparison`: Generates and saves a comparison plot.
- `main`: Executes the full workflow.
=#
using JLD2
using LinearAlgebra
using Statistics
using Plots # Using Plots.jl for simplicity. Add Pkg.add("Plots") if needed.


# --- Configuration ---
const DATA_PATH = "./data.jld2" # Adjust path if needed

# Embedding Parameters for DMDc
const EMBEDDING_LAG_ORIG_STEPS = 60 # Lag (tau) between components in an embedding vector (in original steps)
const NUM_EMBEDDING_DELAYS = 50     # Number of D delays (e.g., D=2 means [x_t, x_{t-tau}, x_{t-2tau}]). Total D+1 blocks.
# Total history spanned by one embedding: NUM_EMBEDDING_DELAYS * EMBEDDING_LAG_ORIG_STEPS original steps.

# ACF/AMI Analysis Parameters
const MAX_ACF_AMI_LAG_ORIG_STEPS = 500 # Max lag for ACF/AMI plots in original time steps.

# Simulation Parameters
const START_IDX = 2000   # Starting index in the *original processed data* (after initial [:,2:end]) for simulation.
const SIM_LEN = 700      # Default length of the simulation in *original* time steps.

# Plotting and Data Selection
const PLOT_STATE_IDX = 1 # Which state component to plot (1-based index into RELEVANT_STATES).
const RELEVANT_STATES = [7] # Relevant states for analysis (original indices from raw data).
const NORMALIZE = false
const NUM_BINS_AMI = 200 # Number of bins for AMI histogram estimation


# --- Core Data Processing and Model Creation Functions ---

"""
Loads X and U data, preprocesses U.
Returns X_orig_processed, U_orig_processed (full resolution after initial selection).
"""
function load_and_preprocess_data(data_path::String)
    println("Loading data from $data_path...")
    try
        X_data_raw = jldopen(data_path, "r")["X_data"]
        X_orig_processed = X_data_raw[RELEVANT_STATES, 2:end]

        U_data_raw = jldopen(data_path, "r")["U_data"]
        U_orig_processed = let
            cols = []
            # Align U with X_orig_processed (which starts from 2nd col of X_raw)
            for u_col_idx in 2:size(U_data_raw, 2)
                u_col = U_data_raw[:, u_col_idx]
                if size(u_col, 1) == 8
                    push!(cols, vec(u_col[1:4] * u_col[5:8]'))
                else
                    # This warning might appear if U_data_raw has different number of columns
                end
            end
            if isempty(cols)
                error("No valid U data columns processed. Check U_data_raw format and alignment.")
            end
            hcat(cols...)
        end

        min_cols = min(size(X_orig_processed, 2), size(U_orig_processed, 2))
        if min_cols == 0
            error("X_orig_processed or U_orig_processed has zero columns. X_raw size: $(size(X_data_raw)), U_raw size: $(size(U_data_raw))")
        end
        X_orig_processed = X_orig_processed[:, 1:min_cols]
        U_orig_processed = U_orig_processed[:, 1:min_cols]

        N_X = size(X_orig_processed, 1)
        N_U = size(U_orig_processed, 1)
        TOTAL_SAMPLES_ORIG_PROCESSED = size(X_orig_processed, 2)
        println("Data loaded and processed: X_orig_processed ($(size(X_orig_processed))), U_orig_processed ($(size(U_orig_processed)))")

        return X_orig_processed, U_orig_processed, N_X, N_U, TOTAL_SAMPLES_ORIG_PROCESSED
    catch e
        println("Error loading and preprocessing data from $data_path: $e")
        rethrow(e)
    end
end

"""
Calculates and plots the Autocorrelation Function (ACF) for a given time series (original resolution).
"""
function calculate_and_plot_acf(time_series_orig::AbstractVector, max_lag_orig::Int, state_label::String, output_dir::String=".")
    println("Calculating ACF for state '$state_label' up to $max_lag_orig original lags...")
    n = length(time_series_orig)
    if n <= max_lag_orig
        println("Warning: Time series length ($n) is not sufficiently greater than max_lag_orig ($max_lag_orig) for ACF. Skipping ACF plot.")
        return
    end

    mean_ts = mean(time_series_orig)
    var_ts = var(time_series_orig, corrected=false)

    if var_ts == 0
        println("Warning: Variance of the time series is zero. Cannot compute ACF. Skipping ACF plot.")
        return
    end

    acf_values = zeros(Float64, max_lag_orig + 1)

    for lag in 0:max_lag_orig
        autocov_sum = 0.0
        for i in 1:(n-lag)
            autocov_sum += (time_series_orig[i] - mean_ts) * (time_series_orig[i+lag] - mean_ts)
        end
        acf_values[lag+1] = (autocov_sum / n) / var_ts
    end
    acf_values[1] = 1.0

    println("Plotting ACF...")
    lags_orig = 0:max_lag_orig
    p = plot(lags_orig, acf_values, seriestype=:bar,
        title="ACF for heatflux",
        xlabel="Lag (s)",
        ylabel="Autocorrelation",
        legend=false,
        xticks=0:round(Int, max_lag_orig / 10):max_lag_orig,
        ylims=(-1, 1.1))

    conf_bound = 1.96 / sqrt(n)
    hline!(p, [conf_bound, -conf_bound], linestyle=:dash, color=:red, label="95% Conf. Bounds")

    figname = joinpath(output_dir, "acf_plot_state_$(state_label)_maxlag_orig$(max_lag_orig).svg")
    try
        savefig(p, figname)
        println("ACF plot saved to $figname")
    catch e
        println("Error saving ACF plot: $e")
    end
end

"""
Calculates Average Mutual Information (AMI) between two series using histograms.
"""
function calculate_ami(series_x::AbstractVector, series_y::AbstractVector, num_bins::Int)
    n_samples = length(series_x)
    if n_samples != length(series_y)
        error("Series_x and series_y must have the same length for AMI.")
    end
    if n_samples == 0
        return 0.0
    end

    min_x, max_x = minimum(series_x), maximum(series_x)
    min_y, max_y = minimum(series_y), maximum(series_y)
    epsilon = 1e-9
    bin_width_x = (max_x - min_x + epsilon) / num_bins
    bin_width_y = (max_y - min_y + epsilon) / num_bins

    if bin_width_x == 0 && bin_width_y == 0
        return 0.0
    elseif bin_width_x == 0
        return 0.0
    elseif bin_width_y == 0
        return 0.0
    end

    joint_counts = zeros(Int, num_bins, num_bins)
    marginal_counts_x = zeros(Int, num_bins)
    marginal_counts_y = zeros(Int, num_bins)

    for i in 1:n_samples
        bin_idx_x = min(num_bins, floor(Int, (series_x[i] - min_x) / bin_width_x) + 1)
        bin_idx_y = min(num_bins, floor(Int, (series_y[i] - min_y) / bin_width_y) + 1)
        bin_idx_x = max(1, bin_idx_x)
        bin_idx_y = max(1, bin_idx_y)
        joint_counts[bin_idx_x, bin_idx_y] += 1
        marginal_counts_x[bin_idx_x] += 1
        marginal_counts_y[bin_idx_y] += 1
    end

    joint_probs = joint_counts ./ n_samples
    marginal_probs_x = marginal_counts_x ./ n_samples
    marginal_probs_y = marginal_counts_y ./ n_samples

    h_x = 0.0
    for p_x in marginal_probs_x
        if p_x > 0
            h_x -= p_x * log2(p_x)
        end
    end
    h_y = 0.0
    for p_y in marginal_probs_y
        if p_y > 0
            h_y -= p_y * log2(p_y)
        end
    end
    h_xy = 0.0
    for i in 1:num_bins, j in 1:num_bins
        p_xy = joint_probs[i, j]
        if p_xy > 0
            h_xy -= p_xy * log2(p_xy)
        end
    end

    ami = h_x + h_y - h_xy
    return max(0.0, ami)
end

"""
Calculates and plots Average Mutual Information (AMI) for a given time series (original resolution).
"""
function calculate_and_plot_ami(time_series_orig::AbstractVector, max_lag_orig::Int, state_label::String, num_bins::Int, output_dir::String=".")
    println("Calculating AMI for heatflux")
    n = length(time_series_orig)
    if n <= max_lag_orig
        println("Warning: Time series length ($n) is not sufficiently greater than max_lag_orig ($max_lag_orig) for AMI. Skipping AMI plot.")
        return
    end

    ami_values = zeros(Float64, max_lag_orig)
    lags_orig = 1:max_lag_orig

    for lag in lags_orig
        series_t = time_series_orig[lag+1:end]
        series_t_minus_lag = time_series_orig[1:end-lag]
        ami_values[lag] = calculate_ami(series_t, series_t_minus_lag, num_bins)
    end

    println("Plotting AMI...")
    p_ami = plot(lags_orig, ami_values, seriestype=:line, marker=:circle,
        title="AMI for heatflux measurements",
        xlabel="Lag (s)",
        ylabel="Average Mutual Information (bits)",
        legend=false,
        xticks=0:round(Int, max_lag_orig / 10):max_lag_orig)

    first_min_lag = -1
    first_min_val = Inf
    if length(ami_values) >= 3
        for i in 2:(length(ami_values)-1)
            if ami_values[i] < ami_values[i-1] && ami_values[i] < ami_values[i+1]
                if ami_values[i] < first_min_val
                    first_min_lag = lags_orig[i]
                    first_min_val = ami_values[i]
                    scatter!(p_ami, [first_min_lag], [first_min_val], marker=:star, color=:red, markersize=8, label="First Local Minimum")
                    annotate!(p_ami, first_min_lag, first_min_val, text("  \$\\tau\$=$first_min_lag", :red, :left, 8))
                    break
                end
            end
        end
    end
    if first_min_lag != -1
        println("First local minimum of AMI found at original lag = $first_min_lag with value = $first_min_val")
    else
        println("No distinct first local minimum found for AMI in the plotted range.")
    end

    figname = joinpath(output_dir, "ami_plot_state_$(state_label)_maxlag_orig$(max_lag_orig)_bins$(num_bins).svg")
    savefig(p_ami, figname,)

end

"""
Creates a delay-embedded matrix from original time series data using specified lags.
`data`: original processed data (N_features x N_samples_orig)
`num_embedding_delays`: D (number of tau-delayed blocks, so D+1 total blocks)
`embedding_lag_orig_steps`: tau (lag in original samples)
Output: `(N_features * (D+1)) x (N_samples_orig - D*tau)` matrix
Each column `j` is [data[:,t_curr]; data[:,t_curr-tau]; ...; data[:,t_curr-D*tau]]
where `t_curr = D*tau + j`.
"""
function delay_embed(data::AbstractMatrix, num_embedding_delays::Int, embedding_lag_orig_steps::Int)
    D = num_embedding_delays
    tau = embedding_lag_orig_steps
    n_features, n_samples_orig = size(data)

    if D < 0
        error("Number of embedding delays must be non-negative.")
    end
    if tau <= 0 && D > 0
        error("Embedding lag must be positive if num_embedding_delays > 0.")
    end

    if D == 0 # No actual delays, just the current state
        return data
    end

    max_history_needed = D * tau # Furthest look-back from the "current" sample in an embedding
    if n_samples_orig <= max_history_needed
        error("Not enough original samples ($n_samples_orig) to create embedding with $D delays and lag $tau (requires > $max_history_needed samples).")
    end

    n_embedded_cols = n_samples_orig - max_history_needed
    embedded_dim_rows = n_features * (D + 1)
    embedded_data = zeros(eltype(data), embedded_dim_rows, n_embedded_cols)

    for j_col in 1:n_embedded_cols
        t_current_in_orig_data = max_history_needed + j_col # Index in `data` for the s(t) part of the j_col-th embedding
        for d_block in 0:D # d_block=0 is s(t), d_block=1 is s(t-tau), ...
            row_start = d_block * n_features + 1
            row_end = (d_block + 1) * n_features
            source_col_in_orig_data = t_current_in_orig_data - d_block * tau
            embedded_data[row_start:row_end, j_col] = data[:, source_col_in_orig_data]
        end
    end
    return embedded_data
end

function input_lift(U::AbstractMatrix)
    return U
end # Placeholder

"""
Calculates DMDc matrices (As, Bs, M) from original processed data using lagged embeddings.
`num_delays` is D, `lag_orig` is tau.
"""
function calculate_dmdc_matrices(X_orig_processed::AbstractMatrix, U_orig_processed::AbstractMatrix,
    num_delays::Int, lag_orig::Int,
    N_X::Int, N_U::Int, TOTAL_SAMPLES_ORIG::Int)
    println("Creating DMDc model with $num_delays delays and lag $lag_orig (original steps)...")
    D = num_delays
    tau = lag_orig

    # Data for constructing X_embed and U_embed (history part)
    # These will use data up to TOTAL_SAMPLES_ORIG - 1
    # The "current" time t in an embedding [x(t); x(t-tau); ...] will go up to TOTAL_SAMPLES_ORIG - 1
    X_history_for_embedding = X_orig_processed[:, 1:(TOTAL_SAMPLES_ORIG-1)]
    U_history_for_embedding = U_orig_processed[:, 1:(TOTAL_SAMPLES_ORIG-1)]

    if size(X_history_for_embedding, 2) <= D * tau || size(U_history_for_embedding, 2) <= D * tau
        error("Not enough samples in history data for embedding. X_hist: $(size(X_history_for_embedding,2)), U_hist: $(size(U_history_for_embedding,2)), D*tau: $(D*tau)")
    end

    X_embed_dmd = delay_embed(X_history_for_embedding, D, tau)
    U_embed_dmd = delay_embed(U_history_for_embedding, D, tau)

    # X_next_dmd (target states x(t+1))
    # If X_embed_dmd[:,j] has its "current" part from X_history_for_embedding[:, t_curr_j],
    # then X_next_dmd[:,j] should be X_orig_processed[:, t_curr_j + 1].
    # t_curr_j for X_embed_dmd[:,j] is (D*tau + j) (index in X_history_for_embedding)
    # So, X_next_dmd[:,j] comes from X_orig_processed[:, (D*tau + j) + 1]

    num_embedded_cols = size(X_embed_dmd, 2) # Should be (TOTAL_SAMPLES_ORIG-1) - D*tau

    start_idx_X_next = (D * tau) + 1 + 1 # Start index in X_orig_processed
    end_idx_X_next = start_idx_X_next + num_embedded_cols - 1

    if end_idx_X_next > TOTAL_SAMPLES_ORIG
        error("Calculated end index for X_next_dmd ($end_idx_X_next) exceeds available samples in X_orig_processed ($TOTAL_SAMPLES_ORIG).")
    end
    X_next_dmd = X_orig_processed[:, start_idx_X_next:end_idx_X_next]


    if size(X_next_dmd, 2) != num_embedded_cols
        error("Dimension mismatch for DMDc: X_next $(size(X_next_dmd,2)), expected $num_embedded_cols based on X_embed.")
    end
    if num_embedded_cols <= 0
        error("Resulting common length for DMDc matrices is non-positive ($num_embedded_cols).")
    end

    U_embed_dmd = abs.(U_embed_dmd) # As per original script's lifting choice

    mean_X_val = NORMALIZE ? mean(X_orig_processed, dims=2) : zeros(N_X, 1)
    X_next_centered = X_next_dmd .- mean_X_val

    X_embed_centered = similar(X_embed_dmd)
    for i in 0:D # For each block in the embedding (D+1 blocks total)
        row_slice = (i*N_X+1):((i+1)*N_X)
        X_embed_centered[row_slice, :] = X_embed_dmd[row_slice, :] .- mean_X_val
    end

    Omega = [X_embed_centered; U_embed_dmd] # U_embed_dmd also has D+1 blocks

    println("Performing DMDc regression (Omega size: $(size(Omega)))...")
    M = X_next_centered * pinv(Omega)
    println("DMDc matrix M size: $(size(M))")

    # M maps [x(t); x(t-tau); ...; x(t-Dtau); u(t); u(t-tau); ...; u(t-Dtau)] to x(t+1)_centered
    A_matrices_flat = M[:, 1:N_X*(D+1)]
    B_matrices_flat = M[:, N_X*(D+1)+1:end]

    As = [A_matrices_flat[:, i*N_X+1:(i+1)*N_X] for i in 0:D] # As[1]=A0, ..., As[D+1]=AD
    Bs = [B_matrices_flat[:, i*N_U+1:(i+1)*N_U] for i in 0:D] # Bs[1]=B0, ..., Bs[D+1]=BD

    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, vec(mean_X_val)
end

"""
Constructs augmented state-space matrices A_full, B_full.
`num_delays` is D. Augmented state uses D delayed X states and D delayed U states.
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, num_delays::Int)
    println("Constructing augmented state-space model (A_full, B_full)...")
    D = num_delays # Number of tau-delayed past states (e.g., D=2 means x_t, x_{t-tau}, x_{t-2tau})
    # As has D+1 elements (A0 to AD), Bs has D+1 elements (B0 to BD)

    # Augmented state: [x_k; x_{k-tau}; ...; x_{k-Dtau}; u_{k-tau}; ...; u_{k-Dtau}]
    # Size of x part: N_X * (D+1)
    # Size of u part: N_U * D
    N_STATES_AUG = N_X * (D + 1) + N_U * D
    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    B_full = zeros(N_STATES_AUG, N_U)

    # Row 1: Dynamics of x_k -> x_{k+1}
    # x_{k+1} = A0*x_k + A1*x_{k-tau} + ... + AD*x_{k-Dtau}  (from As)
    #           + B1*u_{k-tau} + ... + BD*u_{k-Dtau}        (from Bs, shifted index)
    #           + B0*u_k                                   (this goes into B_full)

    # Contribution from past X states to x_{k+1}
    for i in 0:D # As[i+1] is Ai
        A_full[1:N_X, (i*N_X+1):((i+1)*N_X)] = As[i+1]
    end
    # Contribution from past U states (u_{k-tau} to u_{k-Dtau}) to x_{k+1}
    # u_{k-tau} is the first delayed U in augmented state, corresponds to Bs[2] (B1)
    if D > 0
        for i in 1:D # Bs[i+1] is Bi, affecting x_{k+1} via u_{k-i*tau}
            u_col_offset_in_A = N_X * (D + 1) + (i - 1) * N_U # Offset for u_{k-i*tau} block in A_full
            A_full[1:N_X, (u_col_offset_in_A+1):(u_col_offset_in_A+N_U)] = Bs[i+1]
        end
    end
    # Contribution from current u_k to x_{k+1}
    B_full[1:N_X, :] = Bs[1] # B0

    # Shifting X states: x_k -> x_{k-tau} position for next step, etc.
    if D > 0
        for i in 0:(D-1) # x_{k-i*tau} moves to x_{k-(i+1)*tau} slot
            dest_row_start = (i + 1) * N_X + 1 # Slot for x_{k-(i+1)*tau}
            src_col_start = i * N_X + 1     # Slot for x_{k-i*tau}
            A_full[dest_row_start:(dest_row_start+N_X-1), src_col_start:(src_col_start+N_X-1)] = I(N_X)
        end
    end

    # Shifting U inputs: u_k (from B_full) goes to u_{k-tau} slot
    u_target_first_delay_slot_row_start = N_X * (D + 1) + 1
    if D > 0
        B_full[u_target_first_delay_slot_row_start:u_target_first_delay_slot_row_start+N_U-1, :] = I(N_U)
    end
    # Shifting existing delayed U states: u_{k-i*tau} -> u_{k-(i+1)*tau}
    if D > 1
        for i in 1:(D-1)
            src_col_start_u = N_X * (D + 1) + (i - 1) * N_U + 1 # Slot for u_{k-i*tau}
            dest_row_start_u = N_X * (D + 1) + i * N_U + 1    # Slot for u_{k-(i+1)*tau}
            A_full[dest_row_start_u:(dest_row_start_u+N_U-1), src_col_start_u:(src_col_start_u+N_U-1)] = I(N_U)
        end
    end

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end

"""
Main function to create the ABC system from data.
"""
function create_abc_system_from_data(data_path::String, num_emb_delays::Int, emb_lag_orig::Int, max_acf_ami_lag::Int)
    X_orig_processed, U_orig_processed, N_X, N_U, TOTAL_SAMPLES_ORIG = load_and_preprocess_data(data_path)

    original_state_label = "State Unknown"
    if 1 <= PLOT_STATE_IDX <= length(RELEVANT_STATES)
        original_state_label = string(RELEVANT_STATES[PLOT_STATE_IDX])
    else
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is out of bounds. Using generic label.")
        original_state_label = "State (Row $PLOT_STATE_IDX)"
    end

    if N_X > 0 && 1 <= PLOT_STATE_IDX <= N_X
        calculate_and_plot_acf(X_orig_processed[PLOT_STATE_IDX, :], max_acf_ami_lag, original_state_label)
        calculate_and_plot_ami(X_orig_processed[PLOT_STATE_IDX, :], max_acf_ami_lag, original_state_label, NUM_BINS_AMI)
    else
        println("Warning: Cannot plot ACF/AMI as PLOT_STATE_IDX ($PLOT_STATE_IDX) is invalid for X_orig_processed with $N_X states.")
    end

    As, Bs, M, mean_X = calculate_dmdc_matrices(X_orig_processed, U_orig_processed, num_emb_delays, emb_lag_orig, N_X, N_U, TOTAL_SAMPLES_ORIG)
    A_full, B_full = construct_augmented_matrices(As, Bs, N_X, N_U, num_emb_delays)

    # Return X_orig_processed and U_orig_processed for simulation
    return A_full, B_full, As, Bs, M, vec(mean_X), N_X, N_U, X_orig_processed, U_orig_processed, TOTAL_SAMPLES_ORIG
end

# --- Simulation and Reshape Functions ---
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state_aug::AbstractVector, inputs_orig::AbstractMatrix, sim_len_orig::Int)
    n_state_vars_aug = size(A_full, 1)
    n_inputs_u, n_input_steps_orig = size(inputs_orig)
    if n_input_steps_orig < sim_len_orig
        error("Not enough input data for simulation ($sim_len_orig needed, $n_input_steps_orig provided).")
    end
    if size(B_full, 2) != n_inputs_u
        error("B_full input dim ($(size(B_full,2))) mismatch input data dim ($n_inputs_u).")
    end

    states_augmented_orig = zeros(eltype(initial_state_aug), n_state_vars_aug, sim_len_orig)
    current_state_aug = copy(initial_state_aug)
    abc_fn(state, input) = A_full * state + B_full * input

    for k_orig in 1:sim_len_orig
        states_augmented_orig[:, k_orig] = current_state_aug
        current_state_aug = abc_fn(current_state_aug, inputs_orig[:, k_orig])
    end
    return states_augmented_orig
end

function simulate_orig_model(initial_state_orig_model_centered::AbstractVector, inputs_orig::AbstractMatrix, M_dmdc::AbstractMatrix,
    sim_len_orig::Int, num_emb_delays::Int, emb_lag_orig::Int, n_x::Int, n_u::Int)
    # initial_state_orig_model_centered is [X_embed_centered_at_k0; U_embed_delayed_at_k0 (without u_k0)]
    # It should be [x_k0_c; x_{k0-tau}_c; ...; x_{k0-Dtau}_c;  u_{k0-tau}; ...; u_{k0-Dtau}]
    # M_dmdc expects Omega = [x_k_c; ...; x_{k-Dtau}_c; u_k; ...; u_{k-Dtau}]

    states_out_centered_orig = zeros(n_x, sim_len_orig)
    D = num_emb_delays
    tau = emb_lag_orig

    # Deconstruct initial state (which was built like an augmented state for ABC but without u_k)
    current_x_embed_centered = initial_state_orig_model_centered[1:n_x*(D+1)]
    # The U part in initial_state_orig_model_centered is [u_{k0-tau}; ...; u_{k0-Dtau}]
    current_u_past_delayed_embed = D > 0 ? initial_state_orig_model_centered[n_x*(D+1)+1:end] : zeros(0)


    for k_orig_sim_step in 1:sim_len_orig # k_orig_sim_step is the simulation step, not absolute time
        current_input_uk = inputs_orig[:, k_orig_sim_step] # This is u_k for the current step k

        # Construct U_embed_k = [u_k; u_{k-tau}; ...; u_{k-Dtau}]
        u_embed_k = zeros(n_u * (D + 1))
        u_embed_k[1:n_u] = current_input_uk
        if D > 0
            # current_u_past_delayed_embed holds [u_prev_k-tau, ..., u_prev_k-Dtau]
            # For current step k, these become [u_{k-tau}, ..., u_{k-Dtau}]
            u_embed_k[n_u+1:end] = current_u_past_delayed_embed
        end

        u_embed_k_lifted = abs.(u_embed_k) # Apply lifting

        omega_k = [current_x_embed_centered; u_embed_k_lifted]
        x_next_centered = M_dmdc * omega_k
        states_out_centered_orig[:, k_orig_sim_step] = x_next_centered

        # Update current_x_embed_centered for the next iteration:
        # [x_{next}_centered; x_k_centered; ...; x_{k-(D-1)tau}_centered]
        if D > 0
            current_x_embed_centered = [x_next_centered; current_x_embed_centered[1:n_x*D]]
        else # D=0, embedding is just x_k
            current_x_embed_centered = x_next_centered
        end

        # Update current_u_past_delayed_embed for the next iteration:
        # This should become [u_k; u_{k-tau}; ...; u_{k-(D-1)tau}] to be used for next step's u_{k+1-tau} etc.
        if D > 0
            if D == 1
                current_u_past_delayed_embed = current_input_uk # u_k becomes the new u_{k-tau}
            else # D > 1
                current_u_past_delayed_embed = [current_input_uk; current_u_past_delayed_embed[1:n_u*(D-1)]]
            end
        end
    end
    return states_out_centered_orig
end

"""
Creates initial state for ABC model.
Augmented state: [x_k_c; x_{k-tau}_c; ...; x_{k-Dtau}_c; u_{k-tau}; ...; u_{k-Dtau}]
`start_time_idx_orig` is an index in X_orig_processed/U_orig_processed.
"""
function create_abc_init_state(X_orig::AbstractMatrix, U_orig::AbstractMatrix, start_time_idx_orig::Int,
    num_delays::Int, lag_orig::Int, mean_x_vec::AbstractVector, N_X::Int, N_U::Int)
    D = num_delays
    tau = lag_orig
    n_states_aug = N_X * (D + 1) + N_U * D
    init_state_aug = zeros(n_states_aug)

    min_required_hist_idx = start_time_idx_orig - D * tau
    if min_required_hist_idx < 1
        error("Cannot create init state ABC: start_idx $start_time_idx_orig with $D delays and lag $tau needs data to index $min_required_hist_idx < 1.")
    end

    for i in 0:D # x_{k - i*tau}
        current_x_time_idx = start_time_idx_orig - i * tau
        dest_idx_start = i * N_X + 1
        dest_idx_end = (i + 1) * N_X
        init_state_aug[dest_idx_start:dest_idx_end] = X_orig[:, current_x_time_idx] .- mean_x_vec
    end

    if D > 0 # u_{k-tau} to u_{k-Dtau}
        u_part_start_in_aug = N_X * (D + 1) + 1
        for i in 1:D # u_{k - i*tau}
            current_u_time_idx = start_time_idx_orig - i * tau
            dest_idx_start = u_part_start_in_aug + (i - 1) * N_U
            dest_idx_end = dest_idx_start + N_U - 1
            init_state_aug[dest_idx_start:dest_idx_end] = U_orig[:, current_u_time_idx]
        end
    end
    return init_state_aug
end

"""
Creates initial state for Original DMDc model simulation.
State: [x_k_c; x_{k-tau}_c; ...; x_{k-Dtau}_c; u_{k-tau}; ...; u_{k-Dtau}]
This is the state that `simulate_orig_model` deconstructs.
"""
function create_orig_init_state(X_orig::AbstractMatrix, U_orig::AbstractMatrix, start_time_idx_orig::Int,
    num_delays::Int, lag_orig::Int, mean_x_vec::AbstractVector, N_X::Int, N_U::Int)
    # This initial state structure is identical to the ABC one for this setup
    return create_abc_init_state(X_orig, U_orig, start_time_idx_orig, num_delays, lag_orig, mean_x_vec, N_X, N_U)
end

"""
Plots comparison: original non-strided actual data vs. DMDc model predictions (also on original timescale).
"""
function plot_comparison(
    X_orig_processed::AbstractMatrix,
    states_orig_simulated::AbstractMatrix, # DMDc predictions on original timescale
    plot_state_row_idx::Int,
    sim_start_idx_orig::Int, # Starting index of simulation in X_orig_processed
    sim_len_orig::Int,
    num_emb_delays::Int,
    emb_lag_orig::Int,
    output_dir::String="."
)
    println("Plotting results: Original Actual Data vs. DMDc Predictions...")
    original_state_label = (1 <= plot_state_row_idx <= length(RELEVANT_STATES)) ?
                           "X$(RELEVANT_STATES[plot_state_row_idx])" : "X (Row $plot_state_row_idx)"

    # Actual data to plot
    actual_data_time_axis = sim_start_idx_orig:(sim_start_idx_orig+sim_len_orig-1)
    if actual_data_time_axis[end] > size(X_orig_processed, 2)
        println("Warning: Simulation plot range exceeds actual data length. Truncating plot.")
        actual_data_time_axis = sim_start_idx_orig:size(X_orig_processed, 2)
        # sim_len_orig might need adjustment if states_orig_simulated is shorter
    end
    actual_data_to_plot = X_orig_processed[plot_state_row_idx, actual_data_time_axis]

    p_comp = plot(actual_data_time_axis, actual_data_to_plot,
        label="Actual Original Data ($original_state_label)",
        linewidth=2, color=:black)

    # DMDc predictions are for sim_len_orig steps, starting from sim_start_idx_orig
    # The first point in states_orig_simulated corresponds to X_orig_processed[sim_start_idx_orig]
    # No, states_orig_simulated[:,k] is prediction for x(t_start+k).
    # So its time axis is (sim_start_idx_orig+1) : (sim_start_idx_orig+sim_len_orig)
    # Let's adjust: the simulation predicts x(k+1) based on state at k.
    # If init_state is for time t_start, then first prediction is for t_start+1.
    prediction_time_axis = (sim_start_idx_orig+1):(sim_start_idx_orig+sim_len_orig)

    # Ensure prediction_time_axis does not exceed actual_data_time_axis for plotting
    valid_pred_plot_len = min(length(prediction_time_axis), length(actual_data_time_axis) - 1, size(states_orig_simulated, 2))

    if valid_pred_plot_len > 0
        plot!(p_comp, prediction_time_axis[1:valid_pred_plot_len], states_orig_simulated[plot_state_row_idx, 1:valid_pred_plot_len],
            label="Original DMDc Model (Prediction, $original_state_label)",
            linewidth=1.5, linestyle=:dot, color=:blue, marker=:circle, markersize=2)
    else
        println("Warning: No valid prediction points to plot.")
    end

    title!(p_comp, "Model Comparison (Emb Delays=$num_emb_delays, Emb Lag=$emb_lag_orig)")
    xlabel!(p_comp, "Original Time Step (from processed data)")
    ylabel!(p_comp, "State Value")
    plot!(p_comp, legend=:outertopright)

    figname = joinpath(output_dir, "dmdc_comparison_fulldata_D$(num_emb_delays)_lag$(emb_lag_orig)_start$(sim_start_idx_orig).svg")
    try
        savefig(p_comp, figname)
        println("Plot saved to $figname")
    catch e
        println("Error saving comparison plot: $e")
    end
end

# --- Main Execution ---
function main()
    println("Starting main execution with EMBEDDING_LAG_ORIG_STEPS = $EMBEDDING_LAG_ORIG_STEPS, NUM_EMBEDDING_DELAYS = $NUM_EMBEDDING_DELAYS...")
    A_full, B_full, As, Bs, M_dmdc, mean_X_vec, N_X_processed, N_U_processed,
    X_orig_processed, U_orig_processed, TOTAL_SAMPLES_ORIG_PROCESSED =
        create_abc_system_from_data(DATA_PATH, NUM_EMBEDDING_DELAYS, EMBEDDING_LAG_ORIG_STEPS, MAX_ACF_AMI_LAG_ORIG_STEPS)

    sim_start_idx_orig = START_IDX

    # Max history needed for initial state: NUM_EMBEDDING_DELAYS * EMBEDDING_LAG_ORIG_STEPS
    min_allowable_start_idx = NUM_EMBEDDING_DELAYS * EMBEDDING_LAG_ORIG_STEPS + 1
    if sim_start_idx_orig < min_allowable_start_idx
        println("Warning: Original START_IDX ($START_IDX) is too early for embedding history. Adjusting simulation start.")
        sim_start_idx_orig = min_allowable_start_idx
    end

    actual_sim_len_orig = SIM_LEN
    # Simulation predicts one step ahead, so it needs inputs for t_start to t_start + sim_len -1
    # and produces outputs for t_start+1 to t_start + sim_len
    # The actual data for comparison should go up to t_start + sim_len
    if sim_start_idx_orig + actual_sim_len_orig > TOTAL_SAMPLES_ORIG_PROCESSED
        println("Warning: Simulation length makes comparison exceed original data bounds. Adjusting sim_len_orig.")
        actual_sim_len_orig = TOTAL_SAMPLES_ORIG_PROCESSED - sim_start_idx_orig
    end

    if actual_sim_len_orig <= 0
        println("Error: Calculated simulation length is non-positive ($actual_sim_len_orig). Cannot simulate.")
        return
    end

    # Inputs for simulation: from sim_start_idx_orig up to sim_start_idx_orig + actual_sim_len_orig - 1
    sim_inputs_orig = U_orig_processed[:, sim_start_idx_orig:sim_start_idx_orig+actual_sim_len_orig-1]

    println("Creating initial states for simulation starting at original processed data index $sim_start_idx_orig...")
    init_state_abc = create_abc_init_state(X_orig_processed, U_orig_processed, sim_start_idx_orig, NUM_EMBEDDING_DELAYS, EMBEDDING_LAG_ORIG_STEPS, mean_X_vec, N_X_processed, N_U_processed)
    init_state_orig = create_orig_init_state(X_orig_processed, U_orig_processed, sim_start_idx_orig, NUM_EMBEDDING_DELAYS, EMBEDDING_LAG_ORIG_STEPS, mean_X_vec, N_X_processed, N_U_processed)

    println("Running simulations for $actual_sim_len_orig original time steps...")
    # ABC model simulation (still calculated but not plotted by default)
    states_abc_full_augmented_orig = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs_orig, actual_sim_len_orig)
    # states_abc_simulated_orig = states_abc_full_augmented_orig[1:N_X_processed, :] .+ mean_X_vec

    # Original DMDc model simulation
    states_orig_centered_simulated = simulate_orig_model(init_state_orig, sim_inputs_orig, M_dmdc, actual_sim_len_orig, NUM_EMBEDDING_DELAYS, EMBEDDING_LAG_ORIG_STEPS, N_X_processed, N_U_processed)
    states_orig_simulated = states_orig_centered_simulated .+ mean_X_vec

    println("Simulations complete.")

    if 1 <= PLOT_STATE_IDX <= N_X_processed
        plot_comparison(X_orig_processed, states_orig_simulated,
            PLOT_STATE_IDX, sim_start_idx_orig, actual_sim_len_orig, NUM_EMBEDDING_DELAYS, EMBEDDING_LAG_ORIG_STEPS)
    else
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is invalid. Skipping comparison plot.")
    end

    println("Saving model data to ABC.jld2...")
    try
        jldsave("ABC.jld2";
            A=A_full, B=B_full,
            embedding_lag_orig_steps=EMBEDDING_LAG_ORIG_STEPS,
            num_embedding_delays=NUM_EMBEDDING_DELAYS,
            mean_vec=mean_X_vec, M_dmdc=M_dmdc,
            N_X_processed=N_X_processed, N_U_processed=N_U_processed,
            TOTAL_SAMPLES_ORIG_PROCESSED=TOTAL_SAMPLES_ORIG_PROCESSED,
            RELEVANT_STATES_INDICES=RELEVANT_STATES
        )
        println("Model data saved.")
    catch e
        println("Error saving JLD2 file: $e")
    end
    println("Script finished.")
end

# --- Entry Point ---
if abspath(PROGRAM_FILE) == @__FILE__
    if !isfile(DATA_PATH)
        println("Error: Data file not found at $DATA_PATH. Creating a dummy data.jld2 for testing purposes...")
        try
            num_original_states = 10
            num_original_inputs_raw = 8
            min_len_intermediate = START_IDX + SIM_LEN + (NUM_EMBEDDING_DELAYS * EMBEDDING_LAG_ORIG_STEPS) + 100
            num_total_raw_samples = min_len_intermediate + 1

            _X_data_raw = rand(num_original_states, num_total_raw_samples)
            _U_data_raw = rand(num_original_inputs_raw, num_total_raw_samples)
            jldsave(DATA_PATH; X_data=_X_data_raw, U_data=_U_data_raw)
            println("Dummy data.jld2 created. X_data: $(size(_X_data_raw)), U_data: $(size(_U_data_raw)).")
            if maximum(RELEVANT_STATES) > num_original_states
                println("Warning: RELEVANT_STATES ($RELEVANT_STATES) out of bounds for dummy data ($num_original_states states).")
            end
        catch e_dummy
            println("Failed to create dummy data file: $e_dummy")
        end
    end
    main()
end
