using JLD2
using LinearAlgebra
using Statistics
using CairoMakie

# --- Configuration ---
const MODEL_DATA_PATH = "./ABC.jld2"
const ORIGINAL_TIMESERIES_DATA_PATH = "./data.jld2"

const START_INDICES = [1500, 4000, 1700, 1800, 3500, 2000, 3600, 2200, 2300]
const SIM_LEN = 50
const PLOT_STATE_IDX = 1
const OUTPUT_DIR = "./figures/comparison_plots_julia_makie"

# --- Data Loading --- (This function remains the same)
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

# --- Initial State Creation --- (These functions remain the same)
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
            idx_start_in_state = u_embed_start_idx + (i * N_U)
            idx_end_in_state = idx_start_in_state + N_U - 1
            init_omega_minus_current_u[idx_start_in_state:idx_end_in_state] = U[:, start_time-i]
        end
    end
    return init_omega_minus_current_u
end

# --- Simulation Functions --- (These functions remain the same)
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

        current_omega[1:x_embed_len] = new_x_embed_centered
        current_omega[x_embed_len+1+n_u : end] = current_omega[x_embed_len+1 : end-n_u]
    end
    return states_out_centered
end

# --- CairoMakie Plotting Function for Subplot ---
legend_shown = false;
function plot_single_comparison_makie!(ax::Axis,
                                     actual_data_segment::AbstractMatrix,
                                     states_abc::AbstractMatrix,
                                     states_orig::AbstractMatrix,
                                     plot_state_idx::Int,
                                     sim_start_time::Int,
                                     sim_len::Int,
                                     effective_sim_len::Int) # Added effective_sim_len for title

    absolute_time_steps = sim_start_time:(sim_start_time + effective_sim_len - 1) # Use effective_sim_len

    num_selected_states = size(actual_data_segment, 1)
    if plot_state_idx < 1 || plot_state_idx > num_selected_states
        @error "Invalid PLOT_STATE_IDX ($plot_state_idx). Must be between 1 and $num_selected_states."
        # Check if absolute_time_steps is empty before calling mean
        title_pos_x = isempty(absolute_time_steps) ? 0.5 : mean(absolute_time_steps)
        text!(ax, "Invalid PLOT_STATE_IDX", position = (title_pos_x, 0.5), align = (:center, :center), color=:red)
        return
    end

    lines!(ax, absolute_time_steps, actual_data_segment[plot_state_idx, :], label="Simulation", color=:black, linewidth=2.0)
    # lines!(ax, absolute_time_steps, states_abc[plot_state_idx, :], label="ABC", color=:red, linestyle=:dash, linewidth=1.2)
    lines!(ax, absolute_time_steps, states_orig[plot_state_idx, :], label="DMDc", color=:blue, linestyle=:dot, linewidth=1.2)

    # MODIFICATION: Subplot title shows the time range
    ax.title = "$(first(absolute_time_steps)) - $(last(absolute_time_steps))"
    ax.titlesize = 14
    ax.titlefont = "Computer New Modern"

    global legend_shown
    if !legend_shown
        axislegend(ax, labelsize=8, position=:lt, labelfont="Computer New Modern")
    end
    legend_shown = true

end

# --- Main Execution with CairoMakie ---
function main_makie()
    # 1. Load Data
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
    original_state_label = RELEVANT_STATES[PLOT_STATE_IDX]

    # 2. Setup Makie Figure and Grid
    num_plots = length(START_INDICES)
    grid_cols = ceil(Int, sqrt(num_plots))
    grid_rows = ceil(Int, num_plots / grid_cols)

    fig_width = 250 * grid_cols
    fig_height = 200 * grid_rows + 120 # Adjusted for potentially more visible tick labels + titles

    fig = Figure(size = (fig_width, fig_height), figure_padding = (5, 25, 20, 25)) # Increased bottom/left padding

    Label(fig[0, 1:grid_cols], "DMDc vs Simulation Comparison (Delay=$(Int(delay)), heatflux)",
          fontsize=16, font="Computer New Modern", tellwidth=false, halign=:center)

    Label(fig[grid_rows + 1, 1:grid_cols], "Time (s)", fontsize=18, font="Computer New Modern", tellwidth=false, halign=:center, padding=(0,0,10,0)) # Increased bottom padding for label

    Label(fig[1:grid_rows, 0], "Heatflux (W/m²)", fontsize=18, font="Computer New Modern", rotation=pi/2, tellheight=false, valign=:center, padding=(0,10,0,0)) # Increased left padding for label

    axes_grid = Matrix{Any}(undef, grid_rows, grid_cols)

    subplot_counter = 0
    for r in 1:grid_rows
        for c in 1:grid_cols
            subplot_counter += 1
            if subplot_counter > num_plots
                ax = Axis(fig[r,c])
                hidedecorations!(ax)
                hidespines!(ax)
                axes_grid[r,c] = ax
                continue
            end

            sim_start_time = START_INDICES[subplot_counter]
            println("\n--- Processing Start Time: $sim_start_time (Subplot $subplot_counter at [$r, $c]) ---")

            ax = Axis(fig[r, c],
                      xticklabelsize=12, yticklabelsize=12,
                      xticklabelfont="Computer New Modern", yticklabelfont="Computer New Modern"
            )
            axes_grid[r,c] = ax

            # MODIFICATION: Remove hidexdecorations! to show all x-axis ticks/labels
            # if r < grid_rows
            #     hidexdecorations!(ax, grid=false, minorgrid=false)
            # end
            if c > 1
                hideydecorations!(ax, grid=false, minorgrid=false) # Keep for y-axis
            end

            effective_sim_len = SIM_LEN
            min_hist_time = sim_start_time - delay
            max_sim_time = sim_start_time + effective_sim_len - 1

            current_abs_times_for_title = sim_start_time:(sim_start_time + SIM_LEN -1) # For title in case of skip

            if min_hist_time < 1
                println("Warning: Start time $sim_start_time is too early for delay $delay. Skipping.")
                text!(ax, "Skipped:\nStart too early", position = (0.5, 0.5), align = (:center, :center), justification=:center, color=:orange, textsize=10)
                ax.title = "$(first(current_abs_times_for_title)) - $(last(current_abs_times_for_title)) (Skipped)"
                ax.titlesize = 12
                ax.titlefont = "Computer New Modern"
                continue
            end

            if max_sim_time > TOTAL_SAMPLES
                println("Warning: Simulation length extends beyond available data (needs up to $max_sim_time, max is $TOTAL_SAMPLES).")
                effective_sim_len = TOTAL_SAMPLES - sim_start_time + 1
                current_abs_times_for_title = sim_start_time:(sim_start_time + effective_sim_len -1) # Update for title
                if effective_sim_len <= 0
                    println("Warning: Start time $sim_start_time is beyond data range. Skipping.")
                    text!(ax, "Skipped:\nStart beyond data", position = (0.5,0.5), align = (:center, :center), justification=:center, color=:orange, textsize=10)
                    ax.title = "$sim_start_time (Skipped)"
                    ax.titlesize = 10
                    ax.titlefont = "Computer New Modern"
                    continue
                end
                println("Adjusting simulation length to $effective_sim_len.")
            else # If not adjusted, use original SIM_LEN for title consistency if plot proceeds
                 current_abs_times_for_title = sim_start_time:(sim_start_time + SIM_LEN -1)
            end

            if effective_sim_len <= 0
                println("Warning: Calculated simulation length is non-positive ($effective_sim_len). Skipping start time $sim_start_time.")
                text!(ax, "Skipped:\nZero sim length", position = (0.5,0.5), align = (:center, :center), justification=:center, color=:orange, textsize=10)
                # Attempt to make a title even if skipped
                title_range_str = isempty(current_abs_times_for_title) ? "$sim_start_time" : "$(first(current_abs_times_for_title)) - $(last(current_abs_times_for_title))"
                ax.title = "$title_range_str (Skipped)"
                ax.titlesize = 10
                ax.titlefont = "Computer New Modern"
                continue
            end

            sim_indices = sim_start_time:(sim_start_time + effective_sim_len - 1)
            if maximum(sim_indices) > size(U_data, 2)
                println("Error: Simulation indices $sim_indices exceed U_data columns $(size(U_data, 2)). Skipping start time $sim_start_time.")
                text!(ax, "Error:\nU_data bounds", position = (0.5,0.5), align = (:center, :center), justification=:center, color=:red, textsize=10)
                title_range_str = isempty(current_abs_times_for_title) ? "$sim_start_time" : "$(first(current_abs_times_for_title)) - $(last(current_abs_times_for_title))"
                ax.title = "$title_range_str (Error)"
                ax.titlesize = 10
                ax.titlefont = "Computer New Modern"
                continue
            end
            sim_inputs = U_data[:, sim_indices]

            local init_state_abc, init_state_orig_omega
            try
                init_state_abc = create_abc_init_state(X_data, U_data, sim_start_time, delay, mean_X_vec, N_X, N_U)
                init_state_orig_omega = create_orig_init_state(X_data, U_data, sim_start_time, delay, mean_X_vec, N_X, N_U)
            catch e
                println("Error creating initial state for start time $sim_start_time: $e")
                text!(ax, "Error:\nInit state failed", position = (0.5,0.5), align = (:center, :center), justification=:center, color=:red, textsize=10)
                title_range_str = isempty(current_abs_times_for_title) ? "$sim_start_time" : "$(first(current_abs_times_for_title)) - $(last(current_abs_times_for_title))"
                ax.title = "$title_range_str (Error)"
                ax.titlesize = 10
                ax.titlefont = "Computer New Modern"
                continue
            end

            println("Running simulations for start time $sim_start_time...")
            states_abc_aug = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, effective_sim_len)
            states_abc = states_abc_aug[1:N_X, :] .+ mean_X_vec

            states_orig_centered = simulate_orig_model(init_state_orig_omega, sim_inputs, M_dmdc, effective_sim_len, delay, N_X, N_U)
            states_orig = states_orig_centered .+ mean_X_vec
            println("Simulations complete.")

            if maximum(sim_indices) > size(X_data, 2)
                println("Error: Simulation indices $sim_indices exceed selected X_data columns $(size(X_data, 2)). Skipping start time $sim_start_time.")
                text!(ax, "Error:\nX_data bounds", position = (0.5,0.5), align = (:center, :center), justification=:center, color=:red, textsize=10)
                title_range_str = isempty(current_abs_times_for_title) ? "$sim_start_time" : "$(first(current_abs_times_for_title)) - $(last(current_abs_times_for_title))"
                ax.title = "$title_range_str (Error)"
                ax.titlesize = 10
                ax.titlefont = "Computer New Modern"
                continue
            end
            actual_data_segment = X_data[:, sim_indices]

            plot_single_comparison_makie!(ax, actual_data_segment, states_abc, states_orig,
                                          PLOT_STATE_IDX, sim_start_time, SIM_LEN, effective_sim_len) # Pass effective_sim_len
        end
    end

    if !isdir(OUTPUT_DIR)
        println("Creating output directory: $OUTPUT_DIR")
        mkpath(OUTPUT_DIR)
    end

    figname = joinpath(OUTPUT_DIR, "dmdc_comparison_grid_delay$(delay)_state$(original_state_label)_makie_allxticks.png")
    try
        save(figname, fig, px_per_unit=2)
        println("\nCombined Makie plot saved to $figname")
    catch e
        println("\nError saving Makie plot $figname: $e")
    end

    println("\n--- CairoMakie script finished ---")
end

main_makie()
