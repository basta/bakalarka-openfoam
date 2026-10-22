#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays.
It constructs an augmented state-space model (A_full, B_full)
and compares its simulation results against the original DMDc formulation
and the actual data.

It also includes functions to calculate and plot:
- Autocorrelation Function (ACF)
- Average Mutual Information (AMI)
for a selected state.

Key Functions:
- `load_and_preprocess_data`: Loads data from JLD2, preprocesses control input U.
- `calculate_and_plot_acf`: Calculates and plots the ACF.
- `calculate_ami`: Computes AMI between two series using histograms.
- `calculate_and_plot_ami`: Calculates and plots AMI vs. lag.
- `delay_embed`: Creates delay-embedded matrices.
- `calculate_dmdc_matrices`: Computes DMDc system matrices (As, Bs, M) from data.
- `construct_augmented_matrices`: Builds the augmented A_full, B_full matrices.
- `create_abc_system_from_data`: Orchestrates data loading, analysis, DMDc, and ABC matrix construction.
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
const DELAY = 500      # Number of delay steps for states and inputs, also max lag for ACF/AMI
const START_IDX = 2000   # Starting index for simulation comparison
const SIM_LEN = 700      # Default length of the simulation
const PLOT_STATE_IDX = 1 # Which state component to plot (and for ACF/AMI)
const RELEVANT_STATES = [6] # Relevant states for analysis
const NORMALIZE = false
const NUM_BINS_AMI = 32 # Number of bins for AMI histogram estimation


# --- Core Data Processing and Model Creation Functions ---

"""
Loads X and U data, preprocesses U.
"""
function load_and_preprocess_data(data_path::String)
    println("Loading data from $data_path...")
    try
        X_data_raw = jldopen(data_path, "r")["X_data"]
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

        N_X = size(X_data, 1)
        N_U = size(U, 1)
        TOTAL_SAMPLES = size(X_data, 2)
        println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")

        return X_data, U, N_X, N_U, TOTAL_SAMPLES
    catch e
        println("Error loading data from $data_path: $e")
        println("Please ensure the JLD2 file exists, the path is correct, and contains 'X_data' and 'U_data'.")
        rethrow(e)
    end
end

"""
Calculates and plots the Autocorrelation Function (ACF) for a given time series.
"""
function calculate_and_plot_acf(time_series::AbstractVector, max_lag::Int, state_label::String, output_dir::String=".")
    println("Calculating ACF for state '$state_label' up to lag $max_lag...")
    n = length(time_series)
    if n <= max_lag
        println("Warning: Time series length ($n) is not sufficiently greater than max_lag ($max_lag) for ACF calculation. Skipping ACF plot.")
        return
    end

    mean_ts = mean(time_series)
    var_ts = var(time_series, corrected=false)

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
        acf_values[lag+1] = (autocov_sum / n) / var_ts
    end
    acf_values[1] = 1.0

    println("Plotting ACF...")
    lags = 0:max_lag
    p = plot(lags, acf_values, seriestype=:bar,
        title="ACF for State $state_label (Max Lag = $max_lag)",
        xlabel="Lag (\$\\tau\$)",
        ylabel="Autocorrelation",
        legend=false,
        xticks=0:round(Int, max_lag / 10):max_lag,
        ylims=(-1, 1.1))

    conf_bound = 1.96 / sqrt(n)
    hline!(p, [conf_bound, -conf_bound], linestyle=:dash, color=:red, label="95% Conf. Bounds")

    figname = joinpath(output_dir, "acf_plot_state_$(state_label)_delay$(max_lag).png")
    try
        savefig(p, figname)
        println("ACF plot saved to $figname")
    catch e
        println("Error saving ACF plot: $e")
    end
end

"""
Calculates Average Mutual Information (AMI) between two series using histograms.
AMI(X,Y) = H(X) + H(Y) - H(X,Y)
"""
function calculate_ami(series_x::AbstractVector, series_y::AbstractVector, num_bins::Int)
    n_samples = length(series_x)
    if n_samples != length(series_y)
        error("Series_x and series_y must have the same length for AMI calculation.")
    end
    if n_samples == 0
        return 0.0
    end

    min_x, max_x = minimum(series_x), maximum(series_x)
    min_y, max_y = minimum(series_y), maximum(series_y)

    # Add a small epsilon to max to ensure the max value falls into the last bin
    epsilon = 1e-9
    bin_width_x = (max_x - min_x + epsilon) / num_bins
    bin_width_y = (max_y - min_y + epsilon) / num_bins

    if bin_width_x == 0 && bin_width_y == 0 # Both series are constant
        return 0.0 # No uncertainty, no information gain
    elseif bin_width_x == 0 # series_x is constant
        # H(X) = 0, H(X,Y) = H(Y). AMI = 0 + H(Y) - H(Y) = 0
        return 0.0
    elseif bin_width_y == 0 # series_y is constant
        # H(Y) = 0, H(X,Y) = H(X). AMI = H(X) + 0 - H(X) = 0
        return 0.0
    end


    joint_counts = zeros(Int, num_bins, num_bins)
    marginal_counts_x = zeros(Int, num_bins)
    marginal_counts_y = zeros(Int, num_bins)

    for i in 1:n_samples
        # For bin_width_x == 0 case, (series_x[i] - min_x) is 0. Avoid division by zero.
        # This is handled by early exit if bin_width is 0.
        bin_idx_x = min(num_bins, floor(Int, (series_x[i] - min_x) / bin_width_x) + 1)
        bin_idx_y = min(num_bins, floor(Int, (series_y[i] - min_y) / bin_width_y) + 1)

        # Ensure indices are at least 1
        bin_idx_x = max(1, bin_idx_x)
        bin_idx_y = max(1, bin_idx_y)

        joint_counts[bin_idx_x, bin_idx_y] += 1
        marginal_counts_x[bin_idx_x] += 1
        marginal_counts_y[bin_idx_y] += 1
    end

    # Calculate probabilities
    joint_probs = joint_counts ./ n_samples
    marginal_probs_x = marginal_counts_x ./ n_samples
    marginal_probs_y = marginal_counts_y ./ n_samples

    # Calculate entropies
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
    return max(0.0, ami) # AMI should be non-negative; floating point issues might make it slightly negative
end

"""
Calculates and plots Average Mutual Information (AMI) for a given time series.
"""
function calculate_and_plot_ami(time_series::AbstractVector, max_lag::Int, state_label::String, num_bins::Int, output_dir::String=".")
    println("Calculating AMI for state '$state_label' up to lag $max_lag using $num_bins bins...")
    n = length(time_series)
    if n <= max_lag
        println("Warning: Time series length ($n) is not sufficiently greater than max_lag ($max_lag) for AMI. Skipping AMI plot.")
        return
    end

    ami_values = zeros(Float64, max_lag) # AMI for lag 0 is not typically calculated or is H(X)
    lags = 1:max_lag

    for lag in lags
        series_t = time_series[lag+1:end]
        series_t_minus_lag = time_series[1:end-lag]
        ami_values[lag] = calculate_ami(series_t, series_t_minus_lag, num_bins)
    end

    println("Plotting AMI...")
    p_ami = plot(lags, ami_values, seriestype=:line, marker=:circle,
        title="AMI for State $state_label (Max Lag = $max_lag, Bins = $num_bins)",
        xlabel="Lag (\$\\tau\$)",
        ylabel="Average Mutual Information (bits)",
        legend=false,
        xticks=0:round(Int, max_lag / 10):max_lag)

    # Find and mark the first local minimum
    first_min_lag = -1
    first_min_val = Inf
    if length(ami_values) >= 3
        for i in 2:(length(ami_values)-1)
            if ami_values[i] < ami_values[i-1] && ami_values[i] < ami_values[i+1]
                # Check if it's significantly lower or just noise; for now, any local minimum
                if ami_values[i] < first_min_val # Store the first one encountered
                    first_min_lag = lags[i]
                    first_min_val = ami_values[i]
                    # Mark it on the plot
                    scatter!(p_ami, [first_min_lag], [first_min_val], marker=:star, color=:red, markersize=8, label="First Local Minimum")
                    annotate!(p_ami, first_min_lag, first_min_val, text("  \$\\tau\$=$first_min_lag", :red, :left, 8))
                    break # Found the first one
                end
            end
        end
    end
    if first_min_lag != -1
        println("First local minimum of AMI found at lag = $first_min_lag with value = $first_min_val")
    else
        println("No distinct first local minimum found for AMI in the plotted range.")
    end


    figname = joinpath(output_dir, "ami_plot_state_$(state_label)_delay$(max_lag)_bins$(num_bins).png")
    try
        savefig(p_ami, figname)
        println("AMI plot saved to $figname")
    catch e
        println("Error saving AMI plot: $e")
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
    embedded_data = zeros(eltype(data), embedded_dim, n_embedded_samples)

    for t in 1:n_embedded_samples
        for d in 0:delay
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features
            source_col = t + delay - d
            embedded_data[row_start:row_end, t] = data[:, source_col]
        end
    end
    return embedded_data
end

function input_lift(U::AbstractMatrix)
    U_lifted = U
    return U_lifted
end

"""
Calculates DMDc matrices (As, Bs, M) from data.
"""
function calculate_dmdc_matrices(X_data::AbstractMatrix, U::AbstractMatrix, delay::Int, N_X::Int, N_U::Int, TOTAL_SAMPLES::Int)
    println("Creating DMDc model with delay $delay...")
    max_delay = delay

    _X_hist = X_data[:, 1:(TOTAL_SAMPLES-1)]
    _U_hist = U[:, 1:(TOTAL_SAMPLES-1)]

    if size(_X_hist, 2) <= max_delay || size(_U_hist, 2) <= max_delay
        error("Not enough samples for history embedding. X_hist size: $(size(_X_hist,2)), U_hist size: $(size(_U_hist,2)), max_delay: $max_delay")
    end

    X_embed_dmd = delay_embed(_X_hist, max_delay)
    U_embed_dmd = delay_embed(_U_hist, max_delay)

    num_embedded_cols = size(X_embed_dmd, 2)
    X_next_dmd = X_data[:, (max_delay+1)+1:(max_delay+1)+num_embedded_cols]

    common_len = size(X_next_dmd, 2)
    if common_len != size(X_embed_dmd, 2) || common_len != size(U_embed_dmd, 2)
        error("Dimension mismatch after slicing for DMDc: X_next $(size(X_next_dmd,2)), X_embed $(size(X_embed_dmd,2)), U_embed $(size(U_embed_dmd,2))")
    end
    if common_len <= 0
        error("Resulting common length for DMDc matrices is non-positive ($common_len).")
    end

    U_embed_dmd = abs.(U_embed_dmd)

    local mean_X_val
    if NORMALIZE
        mean_X_val = mean(X_data, dims=2)
    else
        mean_X_val = zeros(N_X, 1)
    end
    X_next_centered = X_next_dmd .- mean_X_val

    X_embed_centered = similar(X_embed_dmd)
    for i in 0:max_delay
        row_slice = (i*N_X+1):((i+1)*N_X)
        X_embed_centered[row_slice, :] = X_embed_dmd[row_slice, :] .- mean_X_val
    end

    Omega = [X_embed_centered; U_embed_dmd]

    println("Performing DMDc regression...")
    M = X_next_centered * pinv(Omega)
    println("DMDc matrix M size: $(size(M))")

    A_matrices_flat = M[:, 1:N_X*(max_delay+1)]
    B_matrices_flat = M[:, N_X*(max_delay+1)+1:end]

    As = [A_matrices_flat[:, i*N_X+1:(i+1)*N_X] for i in 0:max_delay]
    Bs = [B_matrices_flat[:, i*N_U+1:(i+1)*N_U] for i in 0:max_delay]

    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, vec(mean_X_val)
end

"""
Constructs augmented state-space matrices A_full, B_full.
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, delay::Int)
    println("Constructing augmented state-space model (A_full, B_full)...")
    D = delay

    N_STATES_AUG = N_X * (D + 1) + N_U * D
    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    B_full = zeros(N_STATES_AUG, N_U)

    for i in 0:D
        A_full[1:N_X, (i*N_X+1):((i+1)*N_X)] = As[i+1]
    end
    if D > 0
        for i in 1:D
            u_col_offset_in_A = N_X * (D + 1) + (i - 1) * N_U
            A_full[1:N_X, (u_col_offset_in_A+1):(u_col_offset_in_A+N_U)] = Bs[i+1]
        end
    end
    B_full[1:N_X, :] = Bs[1]

    if D > 0
        for i in 0:(D-1)
            dest_row_start = (i + 1) * N_X + 1
            dest_row_end = (i + 2) * N_X
            src_col_start = i * N_X + 1
            src_col_end = (i + 1) * N_X
            A_full[dest_row_start:dest_row_end, src_col_start:src_col_end] = I(N_X)
        end
    end

    u_target_first_delay_slot_row_start = N_X * (D + 1) + 1
    if D > 0
        B_full[u_target_first_delay_slot_row_start:u_target_first_delay_slot_row_start+N_U-1, :] = I(N_U)
    end

    if D > 1
        for i in 1:(D-1)
            src_col_start_u = N_X * (D + 1) + (i - 1) * N_U + 1
            src_col_end_u = src_col_start_u + N_U - 1
            dest_row_start_u = N_X * (D + 1) + i * N_U + 1
            dest_row_end_u = dest_row_start_u + N_U - 1
            A_full[dest_row_start_u:dest_row_end_u, src_col_start_u:src_col_end_u] = I(N_U)
        end
    end

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end

"""
Main function to create the ABC system from data, including ACF and AMI analysis.
"""
function create_abc_system_from_data(data_path::String, delay::Int)
    X_data, U, N_X, N_U, TOTAL_SAMPLES = load_and_preprocess_data(data_path)

    original_state_label = "State Unknown"
    if 1 <= PLOT_STATE_IDX <= N_X && PLOT_STATE_IDX <= length(RELEVANT_STATES)
        original_state_label = string(RELEVANT_STATES[PLOT_STATE_IDX])
    else
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is out of bounds for processed X_data with $N_X states or RELEVANT_STATES. Using generic label.")
        original_state_label = "State (Row $PLOT_STATE_IDX)"
    end

    if N_X > 0 && 1 <= PLOT_STATE_IDX <= N_X
        # Call ACF plotting function
        calculate_and_plot_acf(X_data[PLOT_STATE_IDX, :], DELAY, original_state_label)

        # Call AMI plotting function
        calculate_and_plot_ami(X_data[PLOT_STATE_IDX, :], DELAY, original_state_label, NUM_BINS_AMI)
    else
        println("Warning: Cannot plot ACF/AMI as PLOT_STATE_IDX ($PLOT_STATE_IDX) is invalid for processed X_data with $N_X states.")
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
"""
function simulate_orig_model(initial_state_orig_model_centered::AbstractVector,
    inputs::AbstractMatrix,
    M_dmdc::AbstractMatrix,
    sim_length::Int,
    delay::Int,
    n_x::Int,
    n_u::Int)

    states_out_centered = zeros(n_x, sim_length)

    current_x_embed_centered = initial_state_orig_model_centered[1:n_x*(delay+1)]
    current_u_delayed_embed = delay > 0 ? initial_state_orig_model_centered[n_x*(delay+1)+1:end] : zeros(0)

    for k in 1:sim_length
        current_input_uk = inputs[:, k]

        u_embed_k = zeros(n_u * (delay + 1))
        u_embed_k[1:n_u] = current_input_uk
        if delay > 0
            u_embed_k[n_u+1:end] = current_u_delayed_embed
        end

        u_embed_k_lifted = abs.(u_embed_k)

        omega_k = [current_x_embed_centered; u_embed_k_lifted]
        x_next_centered = M_dmdc * omega_k
        states_out_centered[:, k] = x_next_centered

        if delay > 0
            current_x_embed_centered = [x_next_centered; current_x_embed_centered[1:n_x*delay]]
        else
            current_x_embed_centered = x_next_centered
        end

        if delay > 0
            if delay == 1
                current_u_delayed_embed = current_input_uk
            else
                current_u_delayed_embed = [current_input_uk; current_u_delayed_embed[1:n_u*(delay-1)]]
            end
        end
    end
    return states_out_centered
end

"""
Creates the initial state vector for the ABC model.
"""
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time_idx::Int, delay::Int, mean_x_vec::AbstractVector, N_X::Int, N_U::Int)
    D = delay
    n_states_aug = N_X * (D + 1) + N_U * D
    init_state_aug = zeros(n_states_aug)

    min_required_X_idx = start_time_idx - D
    min_required_U_idx = start_time_idx - D

    if min_required_X_idx < 1
        error("Cannot create initial state for ABC: start_time_idx $start_time_idx with delay $D requires X data back to index $min_required_X_idx, which is < 1.")
    end
    if D > 0 && min_required_U_idx < 1
        error("Cannot create initial state for ABC: start_time_idx $start_time_idx with delay $D requires U data back to index $min_required_U_idx, which is < 1.")
    end

    for i in 0:D
        current_x_original_time_idx = start_time_idx - i
        dest_idx_start = i * N_X + 1
        dest_idx_end = (i + 1) * N_X
        init_state_aug[dest_idx_start:dest_idx_end] = X[:, current_x_original_time_idx] .- mean_x_vec
    end

    if D > 0
        u_part_start_in_aug_state = N_X * (D + 1) + 1
        for i in 1:D
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
"""
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time_idx::Int, delay::Int, mean_x_vec::AbstractVector, N_X::Int, N_U::Int)
    D = delay
    init_state_len = N_X * (D + 1) + N_U * D
    init_state_orig = zeros(init_state_len)

    min_required_X_idx = start_time_idx - D
    if min_required_X_idx < 1
        error("Cannot create initial state for Orig: start_time_idx $start_time_idx with delay $D requires X data back to index $min_required_X_idx, which is < 1.")
    end

    for i in 0:D
        current_x_original_time_idx = start_time_idx - i
        dest_idx_start = i * N_X + 1
        dest_idx_end = (i + 1) * N_X
        init_state_orig[dest_idx_start:dest_idx_end] = X[:, current_x_original_time_idx] .- mean_x_vec
    end

    if D > 0
        min_required_U_idx = start_time_idx - D
        if min_required_U_idx < 1
            error("Cannot create initial state for Orig: start_time_idx $start_time_idx with delay $D requires U data back to index $min_required_U_idx, which is < 1.")
        end
        u_part_start_in_state = N_X * (D + 1) + 1
        for i in 1:D
            current_u_original_time_idx = start_time_idx - i
            dest_idx_start = u_part_start_in_state + (i - 1) * N_U
            dest_idx_end = idx_start + N_U - 1
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
    time_axis = 1:sim_len
    original_state_label = "X?"
    if 1 <= plot_state_row_idx <= length(RELEVANT_STATES)
        original_state_label = "X$(RELEVANT_STATES[plot_state_row_idx])"
    else
        original_state_label = "X (Row $plot_state_row_idx)"
    end

    p_comp = plot(time_axis, actual_data_segment[plot_state_row_idx, :],
        label="Actual Data ($original_state_label)",
        linewidth=2, color=:black)
    plot!(p_comp, time_axis, states_abc[plot_state_row_idx, :],
        label="ABC Model ($original_state_label)",
        linewidth=1.5, linestyle=:dash, color=:red)
    plot!(p_comp, time_axis, states_orig[plot_state_row_idx, :],
        label="Original DMDc Model ($original_state_label)",
        linewidth=1.5, linestyle=:dot, color=:blue)

    title!(p_comp, "Model Comparison (Delay = $delay_val)")
    xlabel!(p_comp, "Time Step (relative to original data index $sim_start_abs_time)")
    ylabel!(p_comp, "State Value")
    plot!(p_comp, legend=:outertopright)

    figname = joinpath(output_dir, "dmdc_comparison_delay$(delay_val)_start$(sim_start_abs_time).png")
    try
        savefig(p_comp, figname)
        println("Plot saved to $figname")
    catch e
        println("Error saving comparison plot: $e")
    end
end

# --- Main Execution ---
function main()
    println("Starting main execution...")
    A_full, B_full, As, Bs, M_dmdc, mean_X_vec, N_X_processed, N_U_processed, X_data_processed, U_processed, TOTAL_SAMPLES_processed =
        create_abc_system_from_data(DATA_PATH, DELAY)

    sim_start_time = START_IDX
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

    println("Creating initial states for simulation starting at processed data index $sim_start_time...")
    init_state_abc = create_abc_init_state(X_data_processed, U_processed, sim_start_time, DELAY, mean_X_vec, N_X_processed, N_U_processed)
    init_state_orig = create_orig_init_state(X_data_processed, U_processed, sim_start_time, DELAY, mean_X_vec, N_X_processed, N_U_processed)

    println("Running simulations for $actual_sim_len steps...")
    states_abc_full_augmented = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, actual_sim_len)
    states_abc_simulated = states_abc_full_augmented[1:N_X_processed, :] .+ mean_X_vec

    states_orig_centered_simulated = simulate_orig_model(init_state_orig, sim_inputs, M_dmdc, actual_sim_len, DELAY, N_X_processed, N_U_processed)
    states_orig_simulated = states_orig_centered_simulated .+ mean_X_vec

    actual_data_for_comparison = X_data_processed[:, sim_start_time:sim_start_time+actual_sim_len-1]
    println("Simulations complete.")

    if 1 <= PLOT_STATE_IDX <= N_X_processed
        plot_comparison(actual_data_for_comparison, states_abc_simulated, states_orig_simulated,
            PLOT_STATE_IDX, sim_start_time, actual_sim_len, DELAY)
    else
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is invalid for plotting $N_X_processed processed states. Skipping comparison plot.")
    end

    println("Saving model data to ABC.jld2...")
    try
        jldsave("ABC.jld2";
            A=A_full,
            B=B_full,
            delay=DELAY,
            mean_vec=mean_X_vec,
            M_dmdc=M_dmdc,
            N_X=N_X_processed,
            N_U=N_U_processed,
            TOTAL_SAMPLES_PROCESSED=TOTAL_SAMPLES_processed,
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
        println("Error: Data file not found at $DATA_PATH")
        println("Creating a dummy data.jld2 for testing purposes...")
        try
            num_original_states = 10
            num_original_inputs_raw = 8
            num_total_raw_samples = START_IDX + SIM_LEN + DELAY + 100

            _X_data_raw = rand(num_original_states, num_total_raw_samples)
            _U_data_raw = rand(num_original_inputs_raw, num_total_raw_samples)
            jldsave(DATA_PATH; X_data=_X_data_raw, U_data=_U_data_raw)
            println("Dummy data.jld2 created with X_data: $(size(_X_data_raw)), U_data: $(size(_U_data_raw)).")
            if maximum(RELEVANT_STATES) > num_original_states
                println("Warning: RELEVANT_STATES ($RELEVANT_STATES) is out of bounds for dummy data ($num_original_states states). Adjust RELEVANT_STATES or dummy data.")
            end
        catch e_dummy
            println("Failed to create dummy data file: $e_dummy")
        end
    end
    main()
end
