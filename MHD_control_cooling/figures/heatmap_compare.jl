#=
Script to compare DMDc model simulation with actual system measurements
by generating state heatmaps for multiple prediction horizons.
Outputs are saved as separate SVG files for actual and predicted states
into a specified directory, using a consistent color scale across all plots.
Heatmaps are generated without labels or colorbars for direct use in documents.
Additionally, a single combined figure with all heatmaps (and titles) is produced.
Aligned with the specific data loading of the provided training script.
=#

using JLD2
using LinearAlgebra
using Plots
# Ensure you have a plotting backend, e.g., Plots.gr() or Plots.pyplot()
# Plots.gr() # Example to set a backend

# --- Configuration ---
const ABC_MODEL_FILE_PATH = "./ABC.jld2"
# This is the specific file and keys used by the TRAINING script for X and U
const ACTUAL_MEASUREMENT_SOURCE_FILE = "./openfoam_surface_data_processed_N20xM20.jld2"
const ACTUAL_X_DATA_KEY = "sampled_temperature_matrix_on_grid"
const ACTUAL_U_DATA_KEY = "experimental_input_matrix"

# Output Directory for Figures
const OUTPUT_FIGURE_DIR = "./figures/tpredict/"

# Simulation and Plotting Parameters
const SIM_START_TIME_ABS = 10 # Absolute time index in the measurement data where the
                               # initial condition for the simulation is constructed.
                               # This is time 'k' for x(k) in the initial history.
                               # Must be > loaded_model_delay to allow for history.

const PREDICTION_HORIZONS_STEPS = [5, 10, 20, 40] # Prediction horizons in time steps
const SIM_LEN = maximum(PREDICTION_HORIZONS_STEPS) # Simulation runs for the longest horizon

# Expected state dimension for heatmap
const HEATMAP_ROWS = 20
const HEATMAP_COLS = 20
const EXPECTED_N_X = HEATMAP_ROWS * HEATMAP_COLS # Should match N_X_effective from ABC.jld2

# --- Helper Functions ---

"""
Loads the specific X and U data. Uses the provided source_file, x_key, and u_key arguments.
"""
function load_actual_data_for_comparison(source_file::String, x_key::String, u_key::String)
    println("Loading actual X and U data from $source_file (keys: '$x_key', '$u_key')...")
    local X_data_from_file::AbstractArray
    local U_data_from_file::AbstractArray

    try
        jldopen(source_file, "r") do file
            if !haskey(file, x_key)
                error("Key '$x_key' not found in $source_file")
            end
            loaded_X = file[x_key]
            if !(loaded_X isa AbstractArray)
                error("Data for '$x_key' in $source_file is not an array. Found type: $(typeof(loaded_X))")
            end
            X_data_from_file = loaded_X

            if !haskey(file, u_key)
                error("Key '$u_key' not found in $source_file")
            end
            loaded_U = file[u_key]
            if !(loaded_U isa AbstractArray)
                error("Data for '$u_key' in $source_file is not an array. Found type: $(typeof(loaded_U))")
            end
            U_data_from_file = loaded_U
        end # file is automatically closed here
    catch e
        println("Error opening or reading from JLD2 file: $source_file")
        rethrow(e)
    end

    println("Actual data loaded: X_data_from_file $(size(X_data_from_file)), U_data_from_file $(size(U_data_from_file))")
    return X_data_from_file, U_data_from_file
end

"""
Preprocesses raw U data into the format expected by the model.
"""
function preprocess_U_data_for_model(U_raw::AbstractMatrix)
    println("Preprocessing U data for model...")
    cols = []
    for u_col in eachcol(U_raw)
        if size(u_col, 1) == 8
             push!(cols, vec(u_col[1:4] * u_col[5:8]'))
        else
            println("Warning: Unexpected size for column in U_raw during preprocessing. Expected 8 rows, got $(size(u_col,1)). Skipping column.")
        end
    end
    if isempty(cols)
        error("No U data columns were processed successfully. Check U_raw format or preprocessing logic.")
    end
    return hcat(cols...)
end

"""
Creates the initial augmented state vector for the ABC model.
'start_time_k' is the time index k for the "current" values in the augmented state history.
"""
function create_abc_init_state(X_model_states::AbstractMatrix, U_model_inputs::AbstractMatrix,
                               start_time_k::Int, delay_D::Int, mean_x_vec::AbstractVector,
                               N_X::Int, N_U::Int)

    n_states_aug = N_X * (1 + delay_D) + N_U * delay_D
    init_state_z0 = zeros(n_states_aug)

    earliest_X_idx_needed = start_time_k - delay_D
    earliest_U_idx_needed = start_time_k - delay_D

    if earliest_X_idx_needed < 1
        error("Not enough historical X data for initial state. Need index $earliest_X_idx_needed, but data starts at 1.")
    end
    if start_time_k > size(X_model_states, 2)
        error("start_time_k ($start_time_k) for X initial state is beyond available X samples ($(size(X_model_states, 2))).")
    end

    if delay_D > 0
        if earliest_U_idx_needed < 1
            error("Not enough historical U data for initial state. Need index $earliest_U_idx_needed, but data starts at 1.")
        end
        if (start_time_k - 1) < 1
            error("start_time_k-1 results in an invalid index < 1 for U data history in initial state.")
        end
        if (start_time_k - 1) > size(U_model_inputs, 2)
             error("start_time_k-1 ($start_time_k-1) for U initial state is beyond available U samples ($(size(U_model_inputs, 2))).")
        end
    end

    init_state_z0[1:N_X] = X_model_states[:, start_time_k] .- mean_x_vec

    for i in 1:delay_D
        state_idx_start_in_z = N_X * i + 1
        state_idx_end_in_z = N_X * (i + 1)
        time_idx_in_X = start_time_k - i
        init_state_z0[state_idx_start_in_z:state_idx_end_in_z] = X_model_states[:, time_idx_in_X] .- mean_x_vec
    end

    if delay_D > 0
        u_del_block_start_in_z = N_X * (1 + delay_D) + 1
        for i in 1:delay_D
            input_idx_start_in_z = u_del_block_start_in_z + N_U * (i - 1)
            input_idx_end_in_z = input_idx_start_in_z + N_U - 1
            time_idx_in_U = start_time_k - i
            init_state_z0[input_idx_start_in_z:input_idx_end_in_z] = U_model_inputs[:, time_idx_in_U]
        end
    end
    return init_state_z0
end

"""
Simulates the LTI ABC model.
"""
function simulate_abc_model(A_full_matrix::AbstractMatrix, B_full_matrix::AbstractMatrix,
                            initial_aug_state_z0::AbstractVector, U_input_sequence::AbstractMatrix,
                            simulation_length::Int)

    num_aug_states = size(A_full_matrix, 1)
    num_sim_inputs_avail = size(U_input_sequence, 2)

    if num_sim_inputs_avail < simulation_length
        error("Not enough input U data provided for simulation ($simulation_length steps needed, $num_sim_inputs_avail available).")
    end

    simulated_aug_states_history = zeros(num_aug_states, simulation_length)
    current_aug_state_zk = copy(initial_aug_state_z0)

    for k_sim_step in 1:simulation_length
        current_model_input_uk = U_input_sequence[:, k_sim_step]
        next_aug_state_zkplus1 = A_full_matrix * current_aug_state_zk + B_full_matrix * current_model_input_uk
        simulated_aug_states_history[:, k_sim_step] = next_aug_state_zkplus1
        current_aug_state_zk = next_aug_state_zkplus1
    end
    return simulated_aug_states_history
end

# --- Main Script Logic ---
function main_comparison_script()
    println("Starting heatmap comparison script...")

    if !isdir(OUTPUT_FIGURE_DIR)
        println("Output directory '$OUTPUT_FIGURE_DIR' does not exist. Creating it...")
        mkpath(OUTPUT_FIGURE_DIR)
    end

    println("Loading LTI model parameters from: $ABC_MODEL_FILE_PATH...")
    local A_full, B_full, model_delay, mean_vec, N_X_model, N_U_model

    jldopen(ABC_MODEL_FILE_PATH, "r") do file
        for key in ["A", "B", "delay", "mean_vec", "N_X_effective", "N_U_effective"]
            if !haskey(file, key) error("Key '$key' not found in $ABC_MODEL_FILE_PATH") end
        end
        A_full = file["A"]
        B_full = file["B"]
        model_delay = file["delay"]
        mean_vec = file["mean_vec"]
        N_X_model = file["N_X_effective"]
        N_U_model = file["N_U_effective"]
    end
    println("Model loaded: N_X_model=$N_X_model, N_U_model=$N_U_model, Delay=$model_delay")

    if N_X_model != EXPECTED_N_X
        error("Loaded N_X_model ($N_X_model) from $ABC_MODEL_FILE_PATH does not match EXPECTED_N_X ($EXPECTED_N_X) for $HEATMAP_ROWS x $HEATMAP_COLS heatmap.")
    end

    X_actual_model_states, U_raw_actual = load_actual_data_for_comparison(
        ACTUAL_MEASUREMENT_SOURCE_FILE, ACTUAL_X_DATA_KEY, ACTUAL_U_DATA_KEY
    )
    U_actual_model_inputs = preprocess_U_data_for_model(U_raw_actual)

    if size(X_actual_model_states,1) != N_X_model
        error("Dimension mismatch: Actual X data rows ($(size(X_actual_model_states,1))) vs N_X_model ($N_X_model).")
    end
    if size(U_actual_model_inputs,1) != N_U_model
        error("Dimension mismatch: Processed actual U data rows ($(size(U_actual_model_inputs,1))) vs N_U_model ($N_U_model).")
    end

    total_measurement_samples_X = size(X_actual_model_states, 2)
    total_measurement_samples_U = size(U_actual_model_inputs, 2)

    if SIM_START_TIME_ABS <= model_delay
        error("SIM_START_TIME_ABS ($SIM_START_TIME_ABS) must be > model_delay ($model_delay) to construct initial history for augmented state.")
    end
    if SIM_START_TIME_ABS + SIM_LEN - 1 > total_measurement_samples_U
        error("Simulation period requires U inputs up to index $(SIM_START_TIME_ABS + SIM_LEN - 1), but only $total_measurement_samples_U U-samples available.")
    end
    if SIM_START_TIME_ABS + maximum(PREDICTION_HORIZONS_STEPS) > total_measurement_samples_X
         error("Comparison requires actual X states up to index $(SIM_START_TIME_ABS + maximum(PREDICTION_HORIZONS_STEPS)), but only $total_measurement_samples_X X-samples available.")
    end

    println("Creating initial augmented state for simulation, based on data up to absolute time index $SIM_START_TIME_ABS...")
    initial_aug_state_z0 = create_abc_init_state(
        X_actual_model_states, U_actual_model_inputs,
        SIM_START_TIME_ABS, model_delay, mean_vec, N_X_model, N_U_model
    )

    simulation_input_U_sequence = U_actual_model_inputs[:, SIM_START_TIME_ABS : SIM_START_TIME_ABS + SIM_LEN - 1]

    println("Running simulation for $SIM_LEN steps...")
    simulated_aug_states_history_z = simulate_abc_model(
        A_full, B_full, initial_aug_state_z0, simulation_input_U_sequence, SIM_LEN
    )
    simulated_x_model_states_centered = simulated_aug_states_history_z[1:N_X_model, :]
    simulated_x_model_states = simulated_x_model_states_centered .+ mean_vec

    # --- Determine global min/max for consistent color scaling ---
    println("\nDetermining global min/max for color scaling across all plots...")
    global_min_val = Inf
    global_max_val = -Inf

    for horizon_step in PREDICTION_HORIZONS_STEPS
        actual_data_target_abs_idx = SIM_START_TIME_ABS + horizon_step
        if actual_data_target_abs_idx > total_measurement_samples_X
            println("Warning: Skipping horizon $horizon_step for global min/max calc due to insufficient actual X data.")
            continue
        end

        actual_X_snapshot = X_actual_model_states[:, actual_data_target_abs_idx]
        simulated_X_snapshot = simulated_x_model_states[:, horizon_step]

        global_min_val = min(global_min_val, minimum(actual_X_snapshot), minimum(simulated_X_snapshot))
        global_max_val = max(global_max_val, maximum(actual_X_snapshot), maximum(simulated_X_snapshot))
    end
    println("Global color scale determined: [$global_min_val, $global_max_val]")

    # --- Loop through horizons for plotting ---
    all_individual_plots = [] # To store plot objects for the combined figure

    for horizon_step in PREDICTION_HORIZONS_STEPS
        println("\nProcessing horizon: $horizon_step steps...")

        actual_data_target_abs_idx = SIM_START_TIME_ABS + horizon_step

        if actual_data_target_abs_idx > total_measurement_samples_X
            println("Warning: Actual X data not available for horizon $horizon_step (target abs_idx $actual_data_target_abs_idx). Max X samples: $total_measurement_samples_X. Skipping this horizon's plots.")
            continue
        end

        actual_X_snapshot = X_actual_model_states[:, actual_data_target_abs_idx]
        simulated_X_snapshot = simulated_x_model_states[:, horizon_step]

        actual_heatmap_matrix = reshape(actual_X_snapshot, HEATMAP_ROWS, HEATMAP_COLS)
        simulated_heatmap_matrix = reshape(simulated_X_snapshot, HEATMAP_ROWS, HEATMAP_COLS)

        println("Generating separate heatmaps for target absolute time index $actual_data_target_abs_idx (Horizon: $horizon_step steps)...")

        # --- Plotting for individual label-less files ---
        heatmap_actual_labelless = heatmap(actual_heatmap_matrix,
                                 aspect_ratio=:equal, clims=(global_min_val, global_max_val),
                                 yflip=false, colorbar=false, axis=nothing,
                                 border=:none, framestyle=:none,
                                 size=(HEATMAP_COLS*10,HEATMAP_ROWS*10)
                                 )
        actual_output_filename = joinpath(OUTPUT_FIGURE_DIR, "actual_state_start_t$(SIM_START_TIME_ABS)_horizon$(horizon_step)steps_delay$(model_delay).svg")
        try
            savefig(heatmap_actual_labelless, actual_output_filename)
            println("Actual state (label-less) heatmap saved to: $actual_output_filename")
        catch e
            println("Error saving actual state (label-less) heatmap '$actual_output_filename': $e.")
        end

        heatmap_simulated_labelless = heatmap(simulated_heatmap_matrix,
                                    aspect_ratio=:equal, clims=(global_min_val, global_max_val),
                                    yflip=false, colorbar=false, axis=nothing,
                                    border=:none, framestyle=:none,
                                    size=(HEATMAP_COLS*10,HEATMAP_ROWS*10)
                                    )
        predicted_output_filename = joinpath(OUTPUT_FIGURE_DIR, "predicted_state_start_t$(SIM_START_TIME_ABS)_horizon$(horizon_step)steps_delay$(model_delay).svg")
        try
            savefig(heatmap_simulated_labelless, predicted_output_filename)
            println("Predicted state (label-less) heatmap saved to: $predicted_output_filename")
        catch e
            println("Error saving predicted state (label-less) heatmap '$predicted_output_filename': $e.")
        end

        # --- Plotting for combined figure (with titles) ---
        actual_title_combined = "Actual (t=$actual_data_target_abs_idx)"
        heatmap_actual_for_combined = heatmap(actual_heatmap_matrix,
                                 title=actual_title_combined, titlefontsize=8,
                                 aspect_ratio=:equal, clims=(global_min_val, global_max_val),
                                 yflip=false, colorbar=true, # Keep colorbar for combined context if desired, or set to false
                                 size=(HEATMAP_COLS*15,HEATMAP_ROWS*15) # Slightly larger for title
                                 )
        push!(all_individual_plots, heatmap_actual_for_combined)

        predicted_title_combined = "Pred. (t_init + $horizon_step, target t=$actual_data_target_abs_idx)"
        heatmap_simulated_for_combined = heatmap(simulated_heatmap_matrix,
                                    title=predicted_title_combined, titlefontsize=8,
                                    aspect_ratio=:equal, clims=(global_min_val, global_max_val),
                                    yflip=false, colorbar=true, # Keep colorbar for combined context
                                    size=(HEATMAP_COLS*15,HEATMAP_ROWS*15)
                                    )
        push!(all_individual_plots, heatmap_simulated_for_combined)

    end # end loop over PREDICTION_HORIZONS_STEPS

    # --- Generate and save the combined figure ---
    if !isempty(all_individual_plots)
        println("\nGenerating combined heatmap figure...")
        num_horizons = length(PREDICTION_HORIZONS_STEPS)
        # Adjust layout size dynamically; e.g., 2 columns, num_horizons rows
        # Each subplot was sized individually, so the overall plot size needs to accommodate them.
        # Example: If each subplot is ~300x300, a 4x2 grid might be 700x1300 (width x height)
        combined_plot_width = 2 * (HEATMAP_COLS * 15 + 100) # +100 for spacing/colorbar
        combined_plot_height = num_horizons * (HEATMAP_ROWS*15 + 100) # +100 for spacing/titles

        combined_figure = plot(all_individual_plots...,
                               layout=(num_horizons, 2), # Rows = num_horizons, Columns = 2 (Actual, Predicted)
                               size=(combined_plot_width, combined_plot_height), # Adjust overall size
                               plot_title="Overall State Comparison (t_init = $SIM_START_TIME_ABS, Delay = $model_delay)",
                               plot_titlefontsize=12
                               )
        combined_output_filename = joinpath(OUTPUT_FIGURE_DIR, "combined_heatmap_comparison_start_t$(SIM_START_TIME_ABS)_delay$(model_delay).svg")
        try
            savefig(combined_figure, combined_output_filename)
            println("Combined heatmap figure saved to: $combined_output_filename")
        catch e
            println("Error saving combined heatmap figure '$combined_output_filename': $e.")
        end
    else
        println("\nNo plots were generated for the combined figure (possibly due to data limitations for all horizons).")
    end


    println("\nScript finished.")
end

# --- Script Execution ---
main_comparison_script()
