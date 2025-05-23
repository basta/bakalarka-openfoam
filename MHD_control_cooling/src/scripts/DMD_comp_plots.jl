#=
This script loads pre-computed DMDc model matrices and original time-series data.
- Model data (A, B, M_dmdc, mean, delay, N_X, N_U, RELEVANT_STATES) is loaded from MODEL_DATA_PATH (e.g., "ABC.jld2").
- Original full time-series data (X, U) is loaded from ORIGINAL_TIMESERIES_DATA_PATH (e.g., "data.jld2").

It selects the relevant states from the full X_data based on RELEVANT_STATES,
then simulates both the augmented (ABC) and original DMDc models starting
from multiple specified time points.

It generates two figures:
1. A publication-quality figure containing comparison plots for long-term simulations
   (defined by START_INDICES and SIM_LEN) arranged in a grid, with per-subplot titles, labels, and legends.
2. A publication-quality figure containing comparison plots for many short-term simulations
   (defined by SHORT_TERM_START_INDICES and SHORT_TERM_SIM_LEN) arranged in a grid.
   This figure features a global title, a global y-axis label, a single global legend,
   and each subplot's x-axis will represent the global time index. Subplot titles are removed for this figure.

Both figures compare predictions against the selected actual data using Plots.jl.

Required JLD2 Files:
- "ABC.jld2" (or path set by MODEL_DATA_PATH): Contains :A, :B, :M_dmdc, :delay, :mean_vec, :N_X, :N_U, :RELEVANT_STATES
- "data.jld2" (or path set by ORIGINAL_TIMESERIES_DATA_PATH): Contains the full :X_data, :U_data, and :TOTAL_SAMPLES
=#

using JLD2
using LinearAlgebra
using Statistics
using Plots

# --- Configuration ---

# Input file paths
const MODEL_DATA_PATH = "./ABC.jld2" # File with A, B, M_dmdc, delay, mean_vec, N_X, N_U, RELEVANT_STATES
const ORIGINAL_TIMESERIES_DATA_PATH = "./data.jld2" # File with FULL X_data, U_data

# Simulation and Plotting Parameters for Long-Term Predictions
const START_INDICES = [500, 1000, 1500, 2000] # List of starting time indices (1-based) for simulations
const SIM_LEN = 700                     # Default length of each long-term simulation run

# Simulation and Plotting Parameters for Short-Term Predictions
const SHORT_TERM_START_INDICES = 1000:100:3000 # Range of starting time indices for short-term predictions
const SHORT_TERM_SIM_LEN = 30                # Length of each short-term simulation run

# Common Plotting Parameters
# Index (1-based) of the state variable to plot *within the selected RELEVANT_STATES*.
# E.g., if RELEVANT_STATES=[6], PLOT_STATE_IDX=1 plots the 6th state of the original data.
# If RELEVANT_STATES=[3, 6], PLOT_STATE_IDX=1 plots the 3rd state, PLOT_STATE_IDX=2 plots the 6th state.
const PLOT_STATE_IDX = 1
const OUTPUT_DIR = "./figures/comparison_plots_julia" # Directory to save plots

# Plotting Backend and Style (using GR backend for good defaults and saving options)
gr()
default(
    fontfamily="Computer Modern", # Or another preferred font
    linewidth=1.5,
    markersize=3,
    legendfontsize=8,
    tickfontsize=8,
    guidefontsize=10,
    titlefontsize=10,
    dpi=300, # High resolution for saving
    grid=true,
    framestyle=:box
)

# --- Data Loading ---

"""
Loads model components and original time-series data from specified JLD2 files.
Selects relevant states from X_data based on indices loaded from model file.
"""
function load_all_data(model_path::String, timeseries_path::String)
    println("Loading model data from $model_path...")
    local A_full, B_full, M_dmdc, delay, mean_X_vec, N_X, N_U, RELEVANT_STATES
    try
        jldopen(model_path, "r") do file
            A_full = read(file, "A")
            B_full = read(file, "B")
            delay = read(file, "delay")
            mean_X_vec = read(file, "mean_vec")
            M_dmdc = read(file, "M_dmdc")
            N_X = read(file, "N_X")
            N_U = read(file, "N_U")
            RELEVANT_STATES = read(file, "RELEVANT_STATES")

            if ndims(mean_X_vec) == 1
                mean_X_vec = reshape(mean_X_vec, :, 1)
            end
        end
        println("Model data loaded successfully.")
    catch e
        println("Error loading model data from $model_path: $e")
        println("Please ensure the JLD2 file exists and contains the keys :A, :B, :M_dmdc, :delay, :mean_vec, :N_X, :N_U, :RELEVANT_STATES.")
        rethrow(e)
    end

    delay = Int(delay)
    N_X = Int(N_X)
    N_U = Int(N_U)
    if !(typeof(RELEVANT_STATES) <: AbstractVector || typeof(RELEVANT_STATES) <: AbstractRange)
        try
            RELEVANT_STATES = vec(RELEVANT_STATES)
        catch
            error("RELEVANT_STATES loaded from $model_path is not a vector or range.")
        end
    end


    println("Loading original full time-series data from $timeseries_path...")
    local X_data_full, U_data, TOTAL_SAMPLES
    try
        jldopen(timeseries_path, "r") do file
            X_data_full = read(file, "X_data")
            U_data_raw = read(file, "U_data")

            if N_U == 16 && size(U_data_raw, 1) == 8
                println("Transforming U_data from 8 to 16 dimensions (outer product).")
                U_data_transformed = zeros(16, size(U_data_raw, 2))
                for i in 1:size(U_data_transformed, 2)
                    if size(U_data_raw, 1) >= 8
                        U_data_transformed[:, i] = vec(U_data_raw[1:4, i] * U_data_raw[5:8, i]')
                    else
                        @error "U_data_raw has fewer than 8 rows, inconsistent with expected transformation for N_U=16. Halting."
                        error("Cannot perform 8->16 transformation for U_data due to insufficient rows in raw data.")
                    end
                end
                U_data = U_data_transformed
            else
                if size(U_data_raw, 1) != N_U
                    @warn "Dimension mismatch: N_U from model file is $N_U, but U_data loaded from $timeseries_path has $(size(U_data_raw,1)) rows. Proceeding with loaded U_data dimensions, but this may cause issues if N_U is not accurately reflecting U_data's true dimension for the model."
                end
                U_data = U_data_raw
            end


            if haskey(file, "TOTAL_SAMPLES")
                TOTAL_SAMPLES = read(file, "TOTAL_SAMPLES")
            else
                TOTAL_SAMPLES = size(X_data_full, 2)
                if size(U_data, 2) != TOTAL_SAMPLES
                    error("Inferred TOTAL_SAMPLES from X_data ($(size(X_data_full, 2))) does not match U_data columns ($(size(U_data, 2))) in $timeseries_path")
                end
                @warn "TOTAL_SAMPLES key not found in $timeseries_path. Inferred as $TOTAL_SAMPLES from data dimensions."
            end
        end
        println("Original time-series data loaded successfully.")
    catch e
        println("Error loading original time-series data from $timeseries_path: $e")
        println("Please ensure the JLD2 file exists and contains the keys :X_data, :U_data (and optionally :TOTAL_SAMPLES).")
        rethrow(e)
    end

    TOTAL_SAMPLES = Int(TOTAL_SAMPLES)

    println("Selecting relevant states $(RELEVANT_STATES) from X_data...")
    local X_data
    try
        X_data = X_data_full[RELEVANT_STATES, :]
    catch e
        println("Error selecting RELEVANT_STATES $(RELEVANT_STATES) from X_data_full (size $(size(X_data_full))): $e")
        rethrow(e)
    end
    println("Selected X_data size: $(size(X_data))")


    if size(X_data, 1) != N_X
        error("Inconsistent N_X: Model data ($N_X) vs selected X_data ($(size(X_data, 1))) using RELEVANT_STATES=$RELEVANT_STATES")
    end
    if size(mean_X_vec, 1) != N_X
        error("Inconsistent N_X: mean_X_vec ($(size(mean_X_vec, 1))) vs N_X ($N_X)")
    end

    expected_omega_len = N_X * (delay + 1) + N_U * (delay + 1)
    if size(M_dmdc, 2) != expected_omega_len
        @warn """Dimension mismatch for M_dmdc columns vs expected Omega vector length:
               M_dmdc columns: $(size(M_dmdc, 2))
               Expected Omega length (N_X*(delay+1) + N_U*(delay+1)): $expected_omega_len
               N_X=$N_X, N_U=$N_U, delay=$delay.
               Proceeding with the loaded M_dmdc dimensions."""
    end

    if size(U_data, 1) != N_U
        @warn "Final U_data rows ($(size(U_data, 1))) do not match N_U ($(N_U)) from $model_path. Ensure U_data processing and N_U are consistent."
    end

    n_states_aug_a = N_X * (1 + delay) + N_U * delay
    if size(A_full) != (n_states_aug_a, n_states_aug_a)
        error("A_full dimensions ($(size(A_full))) are inconsistent with N_X=$N_X, N_U=$N_U, delay=$delay (expected ($n_states_aug_a, $n_states_aug_a))")
    end
    if size(B_full) != (n_states_aug_a, N_U)
        error("B_full dimensions ($(size(B_full))) are inconsistent with N_X=$N_X, N_U=$N_U, delay=$delay (expected ($n_states_aug_a, $N_U))")
    end
    if size(X_data, 2) != TOTAL_SAMPLES || size(U_data, 2) != TOTAL_SAMPLES
        error("Selected X_data ($(size(X_data,2))) or U_data ($(size(U_data,2))) length doesn't match TOTAL_SAMPLES ($TOTAL_SAMPLES)")
    end

    return A_full, B_full, M_dmdc, mean_X_vec, delay, N_X, N_U, X_data, U_data, TOTAL_SAMPLES, RELEVANT_STATES
end


# --- Initial State Creation ---
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVecOrMat, N_X::Int, N_U::Int)
    n_states_aug = N_X * (1 + delay) + N_U * delay
    init_state = zeros(n_states_aug)
    min_required_x_time = start_time - delay
    if min_required_x_time < 1 || start_time > size(X, 2)
        error("Cannot create initial ABC state at time $start_time with delay $delay: requires X data in range [$min_required_x_time, $start_time], available X=[1, $(size(X,2))]")
    end
    if delay > 0
        u_min_hist_idx = start_time - delay
        u_max_hist_idx = start_time - 1
        if u_min_hist_idx < 1 || u_max_hist_idx > size(U, 2)
            error("Cannot create initial ABC state at time $start_time with delay $delay: requires U data in range [$(u_min_hist_idx), $(u_max_hist_idx)], available U=[1, $(size(U,2))]")
        end
    end

    init_state[1:N_X] = X[:, start_time] .- mean_x
    for i in 1:delay
        state_idx_start = N_X * i + 1
        state_idx_end = N_X * (i + 1)
        init_state[state_idx_start:state_idx_end] = X[:, start_time-i] .- mean_x
    end
    if delay > 0
        u_del_start_index_in_state = N_X * (1 + delay) + 1
        for i in 1:delay
            input_idx_start = u_del_start_index_in_state + N_U * (i - 1)
            input_idx_end = input_idx_start + N_U - 1
            init_state[input_idx_start:input_idx_end] = U[:, start_time-i]
        end
    end
    return init_state
end

function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVecOrMat, N_X::Int, N_U::Int)
    state_len = N_X * (delay + 1) + N_U * (delay + 1)
    init_omega = zeros(state_len)

    min_required_x_time = start_time - delay
    if min_required_x_time < 1 || start_time > size(X, 2)
        error("Cannot create initial state for original model at time $start_time with delay $delay: requires X data indices [$min_required_x_time, $start_time], available X=[1, $(size(X,2))]")
    end

    if delay > 0
        max_required_u_time = start_time - 1
        min_required_u_time = start_time - delay
        if min_required_u_time < 1 || max_required_u_time > size(U, 2)
            error("Cannot create initial state for original model at time $start_time with delay $delay: requires U data indices [$min_required_u_time, $max_required_u_time] for history, available U=[1, $(size(U,2))]")
        end
    end

    for i in 0:delay
        idx_start = i * N_X + 1
        idx_end = (i + 1) * N_X
        init_omega[idx_start:idx_end] = X[:, start_time-i] .- mean_x
    end

    u_embed_start_idx_in_omega = N_X * (delay + 1) + 1
    if delay > 0
        for i in 1:delay
            idx_in_u_aug_offset = i * N_U
            actual_idx_start = u_embed_start_idx_in_omega + idx_in_u_aug_offset
            actual_idx_end = actual_idx_start + N_U - 1
            init_omega[actual_idx_start:actual_idx_end] = U[:, start_time-i]
        end
    end
    return init_omega
end


# --- Simulation Functions ---
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state::AbstractVector, inputs::AbstractMatrix, sim_length::Int)
    n_state_vars_aug = size(A_full, 1)
    n_input_steps = size(inputs, 2)
    if n_input_steps < sim_length
        error("Not enough input data provided for ABC simulation ($sim_length steps needed, $n_input_steps provided).")
    end
    states_aug = zeros(n_state_vars_aug, sim_length)
    current_state = copy(initial_state)
    for k in 1:sim_length
        current_state = A_full * current_state + B_full * inputs[:, k]
        states_aug[:, k] = current_state
    end
    return states_aug
end

function simulate_orig_model(initial_omega::AbstractVector,
    inputs::AbstractMatrix,
    M_dmdc::AbstractMatrix,
    sim_length::Int,
    delay::Int,
    n_x::Int,
    n_u::Int)

    states_out_centered = zeros(n_x, sim_length)
    current_omega = copy(initial_omega)

    x_aug_len = n_x * (delay + 1)
    u_aug_len = n_u * (delay + 1)

    u_current_start_idx_in_omega = x_aug_len + 1
    u_current_end_idx_in_omega = u_current_start_idx_in_omega + n_u - 1

    if size(M_dmdc, 2) != length(current_omega)
        error("Dimension mismatch: M_dmdc columns ($(size(M_dmdc, 2))) != initial Omega length ($(length(current_omega)))")
    end
    if size(inputs, 2) < sim_length
        error("Not enough input data for original DMDc simulation ($sim_length steps needed, $(size(inputs,2)) provided).")
    end

    for k_sim_step in 1:sim_length
        current_input_for_omega = inputs[:, k_sim_step]
        current_omega[u_current_start_idx_in_omega:u_current_end_idx_in_omega] = current_input_for_omega

        x_next_centered = M_dmdc * current_omega
        states_out_centered[:, k_sim_step] = x_next_centered

        if x_aug_len > n_x
            current_omega[n_x+1:x_aug_len] = current_omega[1:x_aug_len-n_x]
        end
        current_omega[1:n_x] = x_next_centered

        if u_aug_len > n_u
            u_history_target_start = u_current_start_idx_in_omega + n_u
            u_history_source_start = u_current_start_idx_in_omega
            len_to_shift = u_aug_len - n_u
            current_omega[u_history_target_start:u_history_target_start+len_to_shift-1] = current_omega[u_history_source_start:u_history_source_start+len_to_shift-1]
        end
    end
    return states_out_centered
end


# --- Plotting Function for Subplot ---
"""
Plots a single comparison. For short-term plots, x-axis uses global time indices
and subplot titles are omitted.
"""
function plot_single_comparison_subplot!(p::Plots.Plot,
    subplot_idx::Int,
    actual_data_segment::AbstractMatrix,
    states_abc::AbstractMatrix,
    states_orig::AbstractMatrix,
    plot_state_idx::Int,
    sim_start_time::Int,
    sim_len::Int,
    relevant_states::AbstractVector,
    add_legend_labels::Bool;
    is_short_term_plot::Bool=false)

    # For short-term plots, x-axis is global time; otherwise, it's steps from sim_start_time
    time_steps = is_short_term_plot ? (sim_start_time:(sim_start_time+sim_len-1)) : (1:sim_len)

    num_selected_states = size(actual_data_segment, 1)
    if plot_state_idx < 1 || plot_state_idx > num_selected_states
        @error "Invalid PLOT_STATE_IDX ($plot_state_idx) for subplot $subplot_idx. Must be between 1 and $num_selected_states."
        annotate!(p[subplot_idx], (0.5, 0.5), text("Invalid PLOT_STATE_IDX", :red, :center, 8))
        return
    end

    original_state_label = relevant_states[plot_state_idx]

    plot!(p[subplot_idx], time_steps, actual_data_segment[plot_state_idx, :],
        label=add_legend_labels ? "Actual" : "",
        color=:black, linewidth=2.0)

    plot!(p[subplot_idx], time_steps, states_abc[plot_state_idx, :],
        label=add_legend_labels ? "ABC" : "",
        color=:red, linestyle=:dash, linewidth=1.2)

    plot!(p[subplot_idx], time_steps, states_orig[plot_state_idx, :],
        label=add_legend_labels ? "DMDc" : "",
        color=:blue, linestyle=:dot, linewidth=1.2)

    if !is_short_term_plot
        title!(p[subplot_idx], "Start = $sim_start_time")
        xlabel!(p[subplot_idx], "Time Step from Start") # Specific x-label for long-term
        ncols_in_layout = size(p.layout.grid, 2)
        if subplot_idx % ncols_in_layout == 1 || ncols_in_layout == 1
            ylabel!(p[subplot_idx], "State $(original_state_label)")
        end
        plot!(p[subplot_idx], legend=:best)
    else
        # For short-term plots, no individual titles or x/y labels here.
        # Ticks on x-axis will reflect global time due to `time_steps` definition.
        # Ensure grid lines are on if not default for subplots
        plot!(p[subplot_idx], grid=true)
    end
end


# --- Main Execution ---
function main()
    local A_full, B_full, M_dmdc, mean_X_vec, delay, N_X, N_U, X_data, U_data, TOTAL_SAMPLES, RELEVANT_STATES
    try
        A_full, B_full, M_dmdc, mean_X_vec, delay, N_X, N_U, X_data, U_data, TOTAL_SAMPLES, RELEVANT_STATES = load_all_data(MODEL_DATA_PATH, ORIGINAL_TIMESERIES_DATA_PATH)
    catch e
        println("Failed to load and process data: $e")
        return
    end

    num_selected_states = size(X_data, 1)
    if PLOT_STATE_IDX < 1 || PLOT_STATE_IDX > num_selected_states
        @error "Invalid PLOT_STATE_IDX ($PLOT_STATE_IDX). Must be between 1 and $num_selected_states (number of selected states: $(RELEVANT_STATES)). Aborting."
        return
    end
    original_state_label_str = string(RELEVANT_STATES[PLOT_STATE_IDX])

    if !isdir(OUTPUT_DIR)
        println("Creating output directory: $OUTPUT_DIR")
        mkpath(OUTPUT_DIR)
    end

    # --- Generate Long-Term Prediction Figure ---
    println("\n--- Generating Long-Term Prediction Figure (SimLen: $SIM_LEN) ---")
    num_long_term_plots = length(START_INDICES)
    if num_long_term_plots > 0
        grid_cols_long = ceil(Int, sqrt(num_long_term_plots))
        grid_rows_long = ceil(Int, num_long_term_plots / grid_cols_long)
        plot_layout_long = (grid_rows_long, grid_cols_long)
        fig_width_long = 400 * grid_cols_long
        fig_height_long = 300 * grid_rows_long
        long_term_combined_plot = plot(layout=plot_layout_long, size=(fig_width_long, fig_height_long), legend=false)
        long_term_subplot_counter = 0

        for sim_start_time in START_INDICES
            println("\n--- Processing Long-Term, Start Time: $sim_start_time ---")
            long_term_subplot_counter += 1
            effective_sim_len = SIM_LEN
            min_hist_time = sim_start_time - delay
            max_sim_time = sim_start_time + effective_sim_len - 1

            if min_hist_time < 1
                annotate!(long_term_combined_plot[long_term_subplot_counter], (0.5, 0.5), text("Skipped: Start too early", :orange, :center, 8))
                title!(long_term_combined_plot[long_term_subplot_counter], "Start = $sim_start_time")
                continue
            end
            if max_sim_time > TOTAL_SAMPLES
                effective_sim_len = TOTAL_SAMPLES - sim_start_time + 1
                if effective_sim_len <= 0
                    annotate!(long_term_combined_plot[long_term_subplot_counter], (0.5, 0.5), text("Skipped: Start beyond data", :orange, :center, 8))
                    title!(long_term_combined_plot[long_term_subplot_counter], "Start = $sim_start_time")
                    continue
                end
                println("Adjusting sim length to $effective_sim_len for start $sim_start_time.")
            end
            if effective_sim_len <= 0
                annotate!(long_term_combined_plot[long_term_subplot_counter], (0.5, 0.5), text("Skipped: Zero sim length", :orange, :center, 8))
                title!(long_term_combined_plot[long_term_subplot_counter], "Start = $sim_start_time")
                continue
            end

            sim_indices = sim_start_time:(sim_start_time+effective_sim_len-1)
            if maximum(sim_indices) > size(U_data, 2)
                annotate!(long_term_combined_plot[long_term_subplot_counter], (0.5, 0.5), text("Error: U_data bounds", :red, :center, 8))
                title!(long_term_combined_plot[long_term_subplot_counter], "Start = $sim_start_time")
                continue
            end
            sim_inputs = U_data[:, sim_indices]
            if maximum(sim_indices) > size(X_data, 2)
                annotate!(long_term_combined_plot[long_term_subplot_counter], (0.5, 0.5), text("Error: X_data bounds", :red, :center, 8))
                title!(long_term_combined_plot[long_term_subplot_counter], "Start = $sim_start_time")
                continue
            end
            actual_data_segment = X_data[:, sim_indices]

            local init_state_abc, init_state_orig_omega
            try
                init_state_abc = create_abc_init_state(X_data, U_data, sim_start_time, delay, mean_X_vec, N_X, N_U)
                init_state_orig_omega = create_orig_init_state(X_data, U_data, sim_start_time, delay, mean_X_vec, N_X, N_U)
            catch e
                annotate!(long_term_combined_plot[long_term_subplot_counter], (0.5, 0.5), text("Error: Init state failed", :red, :center, 8))
                title!(long_term_combined_plot[long_term_subplot_counter], "Start = $sim_start_time")
                println("Error creating initial state for LT start $sim_start_time: $e")
                continue
            end

            println("Running LT sims for start $sim_start_time (len $effective_sim_len)...")
            states_abc_aug = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, effective_sim_len)
            states_abc = states_abc_aug[1:N_X, :] .+ mean_X_vec
            states_orig_centered = simulate_orig_model(init_state_orig_omega, sim_inputs, M_dmdc, effective_sim_len, delay, N_X, N_U)
            states_orig = states_orig_centered .+ mean_X_vec

            plot_single_comparison_subplot!(long_term_combined_plot, long_term_subplot_counter,
                actual_data_segment, states_abc, states_orig,
                PLOT_STATE_IDX, sim_start_time, effective_sim_len, RELEVANT_STATES, true, is_short_term_plot=false)
        end

        plot!(long_term_combined_plot, plot_title="DMDc vs ABC Model Comparison (Delay=$delay, State=$(original_state_label_str))",
            plot_titlefontsize=14, top_margin=10Plots.mm)
        figname_long = joinpath(OUTPUT_DIR, "dmdc_comparison_grid_delay$(delay)_state$(original_state_label_str).png")
        try
            savefig(long_term_combined_plot, figname_long)
            println("\nLong-term comparison plot saved to $figname_long")
        catch e
            println("\nError saving long-term comparison plot $figname_long: $e")
        end
    else
        println("No start indices defined for long-term predictions. Skipping long-term plot.")
    end


    # --- Generate Short-Term Prediction Figure ---
    println("\n\n--- Generating Short-Term Prediction Figure (SimLen: $SHORT_TERM_SIM_LEN) ---")
    num_short_term_plots = length(SHORT_TERM_START_INDICES)
    if num_short_term_plots > 0
        grid_cols_short = ceil(Int, sqrt(num_short_term_plots))
        grid_rows_short = ceil(Int, num_short_term_plots / grid_cols_short)
        plot_layout_short = (grid_rows_short, grid_cols_short)
        fig_width_short = max(1200, 250 * grid_cols_short)
        fig_height_short = max(900, 200 * grid_rows_short)
        short_term_combined_plot = plot(layout=plot_layout_short, size=(fig_width_short, fig_height_short), legend=:topright)
        short_term_subplot_counter = 0

        for sim_start_time_short in SHORT_TERM_START_INDICES
            println("\n--- Processing Short-Term, Start Time: $sim_start_time_short ---")
            short_term_subplot_counter += 1
            effective_sim_len_short = SHORT_TERM_SIM_LEN
            min_hist_time_short = sim_start_time_short - delay
            max_sim_time_short = sim_start_time_short + effective_sim_len_short - 1

            if min_hist_time_short < 1
                annotate!(short_term_combined_plot[short_term_subplot_counter], (0.5, 0.5), text("Skipped: Start too early", :orange, :center, 6))
                # title!(short_term_combined_plot[short_term_subplot_counter], "Start = $sim_start_time_short") # Title removed for ST plots
                continue
            end
            if max_sim_time_short > TOTAL_SAMPLES
                effective_sim_len_short = TOTAL_SAMPLES - sim_start_time_short + 1
                if effective_sim_len_short <= 0
                    annotate!(short_term_combined_plot[short_term_subplot_counter], (0.5, 0.5), text("Skipped: Start beyond data", :orange, :center, 6))
                    # title!(short_term_combined_plot[short_term_subplot_counter], "Start = $sim_start_time_short")
                    continue
                end
                println("Adjusting ST sim length to $effective_sim_len_short for start $sim_start_time_short.")
            end
            if effective_sim_len_short <= 0
                annotate!(short_term_combined_plot[short_term_subplot_counter], (0.5, 0.5), text("Skipped: Zero sim length", :orange, :center, 6))
                # title!(short_term_combined_plot[short_term_subplot_counter], "Start = $sim_start_time_short")
                continue
            end

            sim_indices_short = sim_start_time_short:(sim_start_time_short+effective_sim_len_short-1)
            if maximum(sim_indices_short) > size(U_data, 2)
                annotate!(short_term_combined_plot[short_term_subplot_counter], (0.5, 0.5), text("Error: U_data bounds", :red, :center, 6))
                # title!(short_term_combined_plot[short_term_subplot_counter], "Start = $sim_start_time_short")
                continue
            end
            sim_inputs_short = U_data[:, sim_indices_short]
            if maximum(sim_indices_short) > size(X_data, 2)
                annotate!(short_term_combined_plot[short_term_subplot_counter], (0.5, 0.5), text("Error: X_data bounds", :red, :center, 6))
                # title!(short_term_combined_plot[short_term_subplot_counter], "Start = $sim_start_time_short")
                continue
            end
            actual_data_segment_short = X_data[:, sim_indices_short]

            local init_state_abc_short, init_state_orig_omega_short
            try
                init_state_abc_short = create_abc_init_state(X_data, U_data, sim_start_time_short, delay, mean_X_vec, N_X, N_U)
                init_state_orig_omega_short = create_orig_init_state(X_data, U_data, sim_start_time_short, delay, mean_X_vec, N_X, N_U)
            catch e
                annotate!(short_term_combined_plot[short_term_subplot_counter], (0.5, 0.5), text("Error: Init state", :red, :center, 6))
                # title!(short_term_combined_plot[short_term_subplot_counter], "Start = $sim_start_time_short")
                println("Error creating initial state for ST start $sim_start_time_short: $e")
                continue
            end

            println("Running ST sims for start $sim_start_time_short (len $effective_sim_len_short)...")
            states_abc_aug_short = simulate_abc_model(A_full, B_full, init_state_abc_short, sim_inputs_short, effective_sim_len_short)
            states_abc_short = states_abc_aug_short[1:N_X, :] .+ mean_X_vec
            states_orig_centered_short = simulate_orig_model(init_state_orig_omega_short, sim_inputs_short, M_dmdc, effective_sim_len_short, delay, N_X, N_U)
            states_orig_short = states_orig_centered_short .+ mean_X_vec

            plot_single_comparison_subplot!(short_term_combined_plot, short_term_subplot_counter,
                actual_data_segment_short, states_abc_short, states_orig_short,
                PLOT_STATE_IDX, sim_start_time_short, effective_sim_len_short, RELEVANT_STATES,
                short_term_subplot_counter == 1, is_short_term_plot=true)
        end

        plot!(short_term_combined_plot, plot_title="Short-Term DMDc vs ABC (Delay=$delay, State=$(original_state_label_str), Pred Len=$(SHORT_TERM_SIM_LEN))",
            plot_titlefontsize=14, top_margin=10Plots.mm)
        ylabel!(short_term_combined_plot, "State $(original_state_label_str) (Global)")

        figname_short = joinpath(OUTPUT_DIR, "dmdc_short_term_grid_delay$(delay)_state$(original_state_label_str)_len$(SHORT_TERM_SIM_LEN).png")
        try
            savefig(short_term_combined_plot, figname_short)
            println("\nShort-term prediction plot saved to $figname_short")
        catch e
            println("\nError saving short-term prediction plot $figname_short: $e")
        end
    else
        println("No start indices defined for short-term predictions. Skipping short-term plot.")
    end

    println("\n--- Script finished ---")
end

# Execute the main function
main()
