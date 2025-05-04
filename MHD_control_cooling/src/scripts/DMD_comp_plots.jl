#=
This script loads pre-computed DMDc model matrices and original time-series data.
- Model data (A, B, M_dmdc, mean, delay, N_X, N_U, RELEVANT_STATES) is loaded from MODEL_DATA_PATH (e.g., "ABC.jld2").
- Original full time-series data (X, U) is loaded from ORIGINAL_TIMESERIES_DATA_PATH (e.g., "data.jld2").

It selects the relevant states from the full X_data based on RELEVANT_STATES,
then simulates both the augmented (ABC) and original DMDc models starting
from multiple specified time points. It generates a single publication-quality
figure containing comparison plots for each start time arranged in a grid,
comparing predictions against the selected actual data using Plots.jl.

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
const ORIGINAL_TIMESERIES_DATA_PATH = "./data.jld2" # File with FULL X_data, U_data, TOTAL_SAMPLES

# Simulation and Plotting Parameters
const START_INDICES = [500, 1000, 1500, 2000] # List of starting time indices (1-based) for simulations
const SIM_LEN = 700                  # Default length of each simulation run
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
    legendfontsize=8, # Adjusted for potentially smaller subplots
    tickfontsize=8,
    guidefontsize=10,
    titlefontsize=10, # Adjusted for subplot titles
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
            U_data_8 = read(file, "U_data")
            U_data_16 = zeros(16, size(U_data_8, 2))
            for i in 1:size(U_data_16, 2)
                U_data_16[:, i] = vec(U_data_8[1:4, i] * U_data_8[5:8, i]')
            end
            U_data = U_data_16
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

    # --- Select Relevant States ---
    println("Selecting relevant states $(RELEVANT_STATES) from X_data...")
    local X_data
    try
        X_data = X_data_full[RELEVANT_STATES, :]
    catch e
        println("Error selecting RELEVANT_STATES $(RELEVANT_STATES) from X_data_full (size $(size(X_data_full))): $e")
        rethrow(e)
    end
    println("Selected X_data size: $(size(X_data))")


    # --- Dimension Checks (After State Selection) ---
    if size(X_data, 1) != N_X
        error("Inconsistent N_X: Model data ($N_X) vs selected X_data ($(size(X_data, 1))) using RELEVANT_STATES=$RELEVANT_STATES")
    end
    if size(mean_X_vec, 1) != N_X
        error("Inconsistent N_X: mean_X_vec ($(size(mean_X_vec, 1))) vs N_X ($N_X)")
    end
    expected_u_dim_m = size(M_dmdc, 2) - N_X * (delay + 1)
    actual_u_embed_dim = N_U * (delay + 1)
    if actual_u_embed_dim != expected_u_dim_m
         @warn """Input dimension mismatch for M_dmdc:
                  N_U ($(N_U)) implies embedded U dim $(actual_u_embed_dim).
                  M_dmdc columns imply embedded U dim $(expected_u_dim_m).
                  Check if N_U in $model_path is correct or if U was lifted non-linearly before creating M_dmdc.
                  Proceeding with N_U = $(N_U)."""
         if size(U_data, 1) != N_U
            @warn "U_data rows ($(size(U_data, 1))) also don't match N_U ($(N_U)). Ensure U_data format is correct for simulation."
         end
    elseif size(U_data, 1) != N_U
         @warn "U_data rows ($(size(U_data, 1))) do not match N_U ($(N_U)) from $model_path, but M_dmdc input dimension is consistent with N_U. Ensure U_data format is correct."
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

"""
Creates the initial state vector for the ABC model.
State: [x(k); x(k-1); ...; x(k-D); u(k-1); ...; u(k-D)] (x centered)
Uses the already selected X_data.
"""
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVecOrMat, N_X::Int, N_U::Int)
    n_states_aug = N_X * (1 + delay) + N_U * delay
    init_state = zeros(n_states_aug)
    max_required_time = start_time
    min_required_time = start_time - delay
    if min_required_time < 1 || max_required_time > size(X, 2) || max_required_time > size(U, 2)
        error("Cannot create initial state at time $start_time with delay $delay: requires data indices [$min_required_time, $max_required_time], available X=[1, $(size(X,2))], U=[1, $(size(U,2))]")
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

"""
Creates the initial state vector for the original DMDc model simulation.
State: [x(k); ... x(k-D) centered; u(k-1); ... u(k-D)]
Uses the already selected X_data.
"""
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVecOrMat, N_X::Int, N_U::Int)
    state_len = N_X * (delay + 1) + N_U * (delay + 1)
    init_omega_minus_current_u = zeros(state_len)
    max_required_time = start_time
    min_required_time = start_time - delay
    if min_required_time < 1 || max_required_time > size(X, 2) || max_required_time > size(U, 2)
        error("Cannot create initial state for original model at time $start_time with delay $delay: requires data indices [$min_required_time, $max_required_time], available X=[1, $(size(X,2))], U=[1, $(size(U,2))]")
    end
    for i in 0:delay
        idx_start = i * N_X + 1
        idx_end = (i + 1) * N_X
        init_omega_minus_current_u[idx_start:idx_end] = X[:, start_time-i] .- mean_x
    end
    u_embed_start_idx = N_X * (delay + 1) + 1
    if delay > 0
        for i in 1:delay
            idx_start = u_embed_start_idx + i * N_U
            idx_end = idx_start + N_U - 1
            init_omega_minus_current_u[idx_start:idx_end] = U[:, start_time-i]
        end
    end
    return init_omega_minus_current_u
end


# --- Simulation Functions ---

"""
Simulates the ABC model (augmented state-space).
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state::AbstractVector, inputs::AbstractMatrix, sim_length::Int)
    n_state_vars_aug = size(A_full, 1)
    n_inputs, n_input_steps = size(inputs)
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

"""
Simulates the original DMDc model formulation using the M matrix.
"""
function simulate_orig_model(initial_omega_minus_current_u::AbstractVector,
                             inputs::AbstractMatrix,
                             M_dmdc::AbstractMatrix,
                             sim_length::Int,
                             delay::Int,
                             n_x::Int,
                             n_u::Int)
    states_out_centered = zeros(n_x, sim_length)
    current_omega = copy(initial_omega_minus_current_u)
    x_embed_len = n_x * (delay + 1)
    u_embed_len = n_u * (delay + 1)
    u_current_start_idx = x_embed_len + 1
    u_current_end_idx = u_current_start_idx + n_u - 1
    if size(M_dmdc, 2) != length(current_omega)
        error("Dimension mismatch: M_dmdc columns ($(size(M_dmdc, 2))) != initial Omega length ($(length(current_omega)))")
    end
    for k in 1:sim_length
        current_input = inputs[:, k]
        current_omega[u_current_start_idx:u_current_end_idx] = current_input
        x_next_centered = M_dmdc * current_omega
        states_out_centered[:, k] = x_next_centered
        new_x_embed_centered = circshift(current_omega[1:x_embed_len], n_x)
        new_x_embed_centered[1:n_x] = x_next_centered
        new_u_embed = circshift(current_omega[x_embed_len+1:end], n_u)
        current_omega[1:x_embed_len] = new_x_embed_centered
        current_omega[x_embed_len+1:end] = new_u_embed
    end
    return states_out_centered
end


# --- Plotting Function for Subplot ---

"""
Plots a single comparison (Actual vs ABC vs Original DMDc) onto a specified subplot
of a larger plot object `p`.
"""
function plot_single_comparison_subplot!(p::Plots.Plot, # The main plot object to modify
                                        subplot_idx::Int,
                                        actual_data_segment::AbstractMatrix,
                                        states_abc::AbstractMatrix, # Un-centered
                                        states_orig::AbstractMatrix, # Un-centered
                                        plot_state_idx::Int,
                                        sim_start_time::Int,
                                        sim_len::Int,
                                        relevant_states::AbstractVector)

    println("Plotting results for start time $sim_start_time onto subplot $subplot_idx...")
    time_steps = 1:sim_len

    num_selected_states = size(actual_data_segment, 1)
    if plot_state_idx < 1 || plot_state_idx > num_selected_states
         @error "Invalid PLOT_STATE_IDX ($plot_state_idx) for subplot $subplot_idx. Must be between 1 and $num_selected_states."
         # Optionally add an annotation to the subplot indicating the error
         annotate!(p[subplot_idx], (0.5, 0.5), text("Invalid PLOT_STATE_IDX", :red, :center, 8))
         return # Skip plotting this subplot
    end

    original_state_label = relevant_states[plot_state_idx]

    # Plot data onto the specified subplot using plot!(p[subplot_idx], ...)
    plot!(p[subplot_idx], time_steps, actual_data_segment[plot_state_idx, :],
          label="Actual", # Keep labels short for subplots
          color=:black, linewidth=2.0) # Slightly thinner line for subplots

    plot!(p[subplot_idx], time_steps, states_abc[plot_state_idx, :],
          label="ABC",
          color=:red, linestyle=:dash, linewidth=1.2)

    plot!(p[subplot_idx], time_steps, states_orig[plot_state_idx, :],
          label="DMDc",
          color=:blue, linestyle=:dot, linewidth=1.2)

    # Add title and labels specific to this subplot
    title!(p[subplot_idx], "Start = $sim_start_time")
    xlabel!(p[subplot_idx], "Time Step")
    # Only add y-label to the first column plots for cleaner look
    if subplot_idx % Int(sqrt(length(p.layout.grid))) == 1 || length(p.layout.grid) <= 2 # Adjust logic based on layout if needed
         ylabel!(p[subplot_idx], "State $(original_state_label)")
    end
    plot!(p[subplot_idx], legend = :best) # Let Plots.jl try to find the best spot
end


# --- Main Execution ---
function main()
    # 1. Load Data
    local A_full, B_full, M_dmdc, mean_X_vec, delay, N_X, N_U, X_data, U_data, TOTAL_SAMPLES, RELEVANT_STATES
    try
        A_full, B_full, M_dmdc, mean_X_vec, delay, N_X, N_U, X_data, U_data, TOTAL_SAMPLES, RELEVANT_STATES = load_all_data(MODEL_DATA_PATH, ORIGINAL_TIMESERIES_DATA_PATH)
    catch e
        println("Failed to load and process data: $e")
        return
    end

    # Determine the original state index for labeling/saving
    num_selected_states = size(X_data, 1)
     if PLOT_STATE_IDX < 1 || PLOT_STATE_IDX > num_selected_states
         @error "Invalid PLOT_STATE_IDX ($PLOT_STATE_IDX). Must be between 1 and $num_selected_states (number of selected states: $(RELEVANT_STATES)). Aborting."
         return
    end
    original_state_label = RELEVANT_STATES[PLOT_STATE_IDX]


    # 2. Setup Combined Plot
    num_plots = length(START_INDICES)
    grid_cols = ceil(Int, sqrt(num_plots))
    grid_rows = ceil(Int, num_plots / grid_cols)
    plot_layout = (grid_rows, grid_cols)
    # Adjust figure size based on grid size (heuristic)
    fig_width = 400 * grid_cols
    fig_height = 300 * grid_rows
    combined_plot = plot(layout = plot_layout, size = (fig_width, fig_height), legend=false) # Initialize empty plot with layout, disable global legend
    subplot_counter = 0 # To track which subplot to plot on

    # 3. Loop Through Start Indices and Simulate/Plot onto Subplots
    for sim_start_time in START_INDICES
        println("\n--- Processing Start Time: $sim_start_time ---")
        subplot_counter += 1 # Increment subplot index

        # Simulation Setup & Boundary Checks
        effective_sim_len = SIM_LEN
        min_hist_time = sim_start_time - delay
        max_sim_time = sim_start_time + effective_sim_len - 1

        if min_hist_time < 1
            println("Warning: Start time $sim_start_time is too early for delay $delay. Skipping.")
            annotate!(combined_plot[subplot_counter], (0.5, 0.5), text("Skipped: Start too early", :orange, :center, 8))
            title!(combined_plot[subplot_counter], "Start = $sim_start_time") # Add title even if skipped
            continue
        end
        if max_sim_time > TOTAL_SAMPLES
            println("Warning: Simulation length extends beyond available data (needs up to $max_sim_time, max is $TOTAL_SAMPLES).")
            effective_sim_len = TOTAL_SAMPLES - sim_start_time + 1
            if effective_sim_len <= 0
                 println("Warning: Start time $sim_start_time is beyond data range. Skipping.")
                 annotate!(combined_plot[subplot_counter], (0.5, 0.5), text("Skipped: Start beyond data", :orange, :center, 8))
                 title!(combined_plot[subplot_counter], "Start = $sim_start_time")
                 continue
            end
            println("Adjusting simulation length to $effective_sim_len.")
        end
        if effective_sim_len <= 0
             println("Warning: Calculated simulation length is non-positive ($effective_sim_len). Skipping start time $sim_start_time.")
             annotate!(combined_plot[subplot_counter], (0.5, 0.5), text("Skipped: Zero sim length", :orange, :center, 8))
             title!(combined_plot[subplot_counter], "Start = $sim_start_time")
             continue
        end

        sim_indices = sim_start_time:(sim_start_time + effective_sim_len - 1)
        if maximum(sim_indices) > size(U_data, 2)
            println("Error: Simulation indices $sim_indices exceed U_data columns $(size(U_data, 2)). Skipping start time $sim_start_time.")
            annotate!(combined_plot[subplot_counter], (0.5, 0.5), text("Error: U_data bounds", :red, :center, 8))
            title!(combined_plot[subplot_counter], "Start = $sim_start_time")
            continue
        end
        sim_inputs = U_data[:, sim_indices]

        # Create Initial States
        local init_state_abc, init_state_orig_omega
        try
             init_state_abc = create_abc_init_state(X_data, U_data, sim_start_time, delay, mean_X_vec, N_X, N_U)
             init_state_orig_omega = create_orig_init_state(X_data, U_data, sim_start_time, delay, mean_X_vec, N_X, N_U)
        catch e
             println("Error creating initial state for start time $sim_start_time: $e")
             annotate!(combined_plot[subplot_counter], (0.5, 0.5), text("Error: Init state failed", :red, :center, 8))
             title!(combined_plot[subplot_counter], "Start = $sim_start_time")
             continue
        end

        # Run Simulations
        println("Running simulations for start time $sim_start_time...")
        states_abc_aug = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, effective_sim_len)
        states_abc = states_abc_aug[1:N_X, :] .+ mean_X_vec

        states_orig_centered = simulate_orig_model(init_state_orig_omega, sim_inputs, M_dmdc, effective_sim_len, delay, N_X, N_U)
        states_orig = states_orig_centered .+ mean_X_vec
        println("Simulations complete.")

        # Get Actual Data for Comparison
         if maximum(sim_indices) > size(X_data, 2)
            println("Error: Simulation indices $sim_indices exceed selected X_data columns $(size(X_data, 2)). Skipping start time $sim_start_time.")
             annotate!(combined_plot[subplot_counter], (0.5, 0.5), text("Error: X_data bounds", :red, :center, 8))
             title!(combined_plot[subplot_counter], "Start = $sim_start_time")
            continue
        end
        actual_data_segment = X_data[:, sim_indices]

        # Plot results onto the current subplot
        plot_single_comparison_subplot!(combined_plot, subplot_counter,
                                        actual_data_segment, states_abc, states_orig,
                                        PLOT_STATE_IDX, sim_start_time, effective_sim_len, RELEVANT_STATES)

    end

    # 4. Finalize and Save Combined Plot
    # Add an overall title (adjust top margin if needed)
    plot!(combined_plot, plot_title = "DMDc vs ABC Model Comparison (Delay=$delay, State=$(original_state_label))",
          plot_titlefontsize=14, top_margin=10Plots.mm) # Add some margin for the main title

    # Ensure output directory exists
    if !isdir(OUTPUT_DIR)
        println("Creating output directory: $OUTPUT_DIR")
        mkpath(OUTPUT_DIR)
    end

    # Save the combined plot
    figname = joinpath(OUTPUT_DIR, "dmdc_comparison_grid_delay$(delay)_state$(original_state_label).png")
    try
        savefig(combined_plot, figname)
        println("\nCombined plot saved to $figname")
    catch e
        println("\nError saving combined plot $figname: $e")
    end

    println("\n--- Script finished ---")
end

# Execute the main function
main()
