#=
Script to generate an animation sequence comparing DMDc model simulation
with actual system measurements. For each time step, it produces a single
figure containing side-by-side heatmaps of the actual and predicted
temperature fields.

Outputs are saved as sequentially numbered image files (e.g., frame_0001.svg)
into a specified directory, ready for compilation into a video. A consistent
color scale is used across all frames for accurate visual comparison.
=#

using JLD2
using LinearAlgebra
using Plots
using Printf
using ProgressMeter

# Ensure you have a plotting backend. GR is a good default.
# To install, run in Julia REPL: using Pkg; Pkg.add("GR")
gr()

# --- Configuration ---
const ABC_MODEL_FILE_PATH = "./ABC.jld2"
# This is the specific file and keys used by the TRAINING script for X and U
const ACTUAL_MEASUREMENT_SOURCE_FILE = "./openfoam_surface_data_processed_N20xM20.jld2"
const ACTUAL_X_DATA_KEY = "sampled_temperature_matrix_on_grid"
const ACTUAL_U_DATA_KEY = "experimental_input_matrix"

# Output Directory for Animation Frames
const OUTPUT_ANIMATION_DIR = "./figures/tpredict_animation/"

# Simulation and Animation Parameters
const SIM_START_TIME_ABS = 10 # Absolute time index in the measurement data for the initial condition.
# This is time 'k' for x(k) in the augmented state's history.
# Must be > loaded_model_delay to allow for history construction.

const ANIMATION_LENGTH = 100 # Number of frames/time steps to simulate and generate.

# Expected state dimension for heatmap
const HEATMAP_ROWS = 20
const HEATMAP_COLS = 20
const EXPECTED_N_X = HEATMAP_ROWS * HEATMAP_COLS # Should match N_X_effective from ABC.jld2

# --- Helper Functions (largely unchanged from original script) ---

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
                error("Data for '$x_key' is not an array.")
            end
            X_data_from_file = loaded_X

            if !haskey(file, u_key)
                error("Key '$u_key' not found in $source_file")
            end
            loaded_U = file[u_key]
            if !(loaded_U isa AbstractArray)
                error("Data for '$u_key' is not an array.")
            end
            U_data_from_file = loaded_U
        end
    catch e
        println("Error opening or reading from JLD2 file: $source_file")
        rethrow(e)
    end

    println("Actual data loaded: X_data $(size(X_data_from_file)), U_data $(size(U_data_from_file))")
    return X_data_from_file, U_data_from_file
end

"""
Preprocesses raw U data into the format expected by the model.
"""
function preprocess_U_data_for_model(U_raw::AbstractMatrix)
    println("Preprocessing U data for model...")
    cols = []
    for u_col in eachcol(U_raw)
        # This function is specific to the user's data structure
        if size(u_col, 1) == 8
            push!(cols, vec(u_col[1:4] * u_col[5:8]'))
        else
            println("Warning: Unexpected size for column in U_raw during preprocessing. Skipping.")
        end
    end
    if isempty(cols)
        error("No U data columns were processed successfully.")
    end
    return hcat(cols...)
end

"""
Creates the initial augmented state vector for the ABC model.
"""
function create_abc_init_state(X_model_states::AbstractMatrix, U_model_inputs::AbstractMatrix,
    start_time_k::Int, delay_D::Int, mean_x_vec::AbstractVector,
    N_X::Int, N_U::Int)

    n_states_aug = N_X * (1 + delay_D) + N_U * delay_D
    init_state_z0 = zeros(n_states_aug)

    # Check if there is enough historical data available
    if (start_time_k - delay_D) < 1
        error("Not enough historical data for initial state. Need data back to index $(start_time_k - delay_D).")
    end

    init_state_z0[1:N_X] = X_model_states[:, start_time_k] .- mean_x_vec

    for i in 1:delay_D
        state_idx_start = N_X * i + 1
        state_idx_end = N_X * (i + 1)
        init_state_z0[state_idx_start:state_idx_end] = X_model_states[:, start_time_k-i] .- mean_x_vec
    end

    if delay_D > 0
        u_del_block_start = N_X * (1 + delay_D) + 1
        for i in 1:delay_D
            input_idx_start = u_del_block_start + N_U * (i - 1)
            input_idx_end = input_idx_start + N_U - 1
            init_state_z0[input_idx_start:input_idx_end] = U_model_inputs[:, start_time_k-i]
        end
    end
    return init_state_z0
end

"""
Simulates the LTI ABC model for a given number of steps.
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix,
    initial_aug_state_z0::AbstractVector, U_input_sequence::AbstractMatrix,
    simulation_length::Int)

    num_aug_states = size(A_full, 1)
    if size(U_input_sequence, 2) < simulation_length
        error("Not enough input U data for simulation ($simulation_length steps needed, $(size(U_input_sequence, 2)) available).")
    end

    simulated_aug_states_history = zeros(num_aug_states, simulation_length)
    current_aug_state_zk = copy(initial_aug_state_z0)

    for k_sim in 1:simulation_length
        current_input_uk = U_input_sequence[:, k_sim]
        next_aug_state_zkplus1 = A_full * current_aug_state_zk + B_full * current_input_uk
        simulated_aug_states_history[:, k_sim] = next_aug_state_zkplus1
        current_aug_state_zk = next_aug_state_zkplus1
    end
    return simulated_aug_states_history
end

# --- Main Script Logic ---
function main_animation_script()
    println("Starting heatmap animation generation script...")

    if isdir(OUTPUT_ANIMATION_DIR)
        println("Warning: Output directory '$OUTPUT_ANIMATION_DIR' already exists. Files may be overwritten.")
    else
        println("Creating output directory: '$OUTPUT_ANIMATION_DIR'...")
        mkpath(OUTPUT_ANIMATION_DIR)
    end

    println("Loading LTI model from: $ABC_MODEL_FILE_PATH...")
    local A_full, B_full, model_delay, mean_vec, N_X_model, N_U_model
    jldopen(ABC_MODEL_FILE_PATH, "r") do file
        A_full = file["A"]
        B_full = file["B"]
        model_delay = file["delay"]
        mean_vec = file["mean_vec"]
        N_X_model = file["N_X_effective"]
        N_U_model = file["N_U_effective"]
    end
    println("Model loaded: N_X=$N_X_model, N_U=$N_U_model, Delay=$model_delay")
    if N_X_model != EXPECTED_N_X
        error("Model N_X ($N_X_model) does not match expected ($EXPECTED_N_X).")
    end

    X_actual, U_raw_actual = load_actual_data_for_comparison(
        ACTUAL_MEASUREMENT_SOURCE_FILE, ACTUAL_X_DATA_KEY, ACTUAL_U_DATA_KEY
    )
    U_actual = preprocess_U_data_for_model(U_raw_actual)

    # Validate data dimensions and simulation length
    if size(X_actual, 1) != N_X_model
        error("Dimension mismatch: Actual X data rows vs N_X_model.")
    end
    if size(U_actual, 1) != N_U_model
        error("Dimension mismatch: Processed actual U data rows vs N_U_model.")
    end
    if SIM_START_TIME_ABS <= model_delay
        error("SIM_START_TIME_ABS must be > model_delay.")
    end
    if SIM_START_TIME_ABS + ANIMATION_LENGTH > size(X_actual, 2)
        error("Animation period exceeds available actual X data. Needed up to index $(SIM_START_TIME_ABS + ANIMATION_LENGTH).")
    end
    if SIM_START_TIME_ABS + ANIMATION_LENGTH - 1 > size(U_actual, 2)
        error("Animation period exceeds available U inputs. Needed up to index $(SIM_START_TIME_ABS + ANIMATION_LENGTH - 1).")
    end


    println("Creating initial state for simulation at absolute time t=$SIM_START_TIME_ABS...")
    initial_aug_state_z0 = create_abc_init_state(
        X_actual, U_actual, SIM_START_TIME_ABS, model_delay, mean_vec, N_X_model, N_U_model
    )

    simulation_input_U = U_actual[:, SIM_START_TIME_ABS:(SIM_START_TIME_ABS+ANIMATION_LENGTH-1)]

    println("Running simulation for $ANIMATION_LENGTH steps...")
    simulated_aug_states_z = simulate_abc_model(
        A_full, B_full, initial_aug_state_z0, simulation_input_U, ANIMATION_LENGTH
    )
    # Extract the state part (x) and add the mean back
    simulated_x_states = simulated_aug_states_z[1:N_X_model, :] .+ mean_vec

    println("\nDetermining global min/max for consistent color scaling...")
    actual_x_to_compare = X_actual[:, (SIM_START_TIME_ABS+1):(SIM_START_TIME_ABS+ANIMATION_LENGTH)]
    global_min_val = min(minimum(actual_x_to_compare), minimum(simulated_x_states))
    global_max_val = max(maximum(actual_x_to_compare), maximum(simulated_x_states))
    println("Global color scale determined: [$global_min_val, $global_max_val]")

    println("\nGenerating $ANIMATION_LENGTH animation frames...")
    for k_frame in 1:ANIMATION_LENGTH
        actual_time_idx = SIM_START_TIME_ABS + k_frame

        actual_X_snapshot = X_actual[:, actual_time_idx]
        simulated_X_snapshot = simulated_x_states[:, k_frame]

        actual_heatmap_matrix = reshape(actual_X_snapshot, HEATMAP_ROWS, HEATMAP_COLS)
        simulated_heatmap_matrix = reshape(simulated_X_snapshot, HEATMAP_ROWS, HEATMAP_COLS)

        # Plot for the 'Actual' state
        p_actual = heatmap(actual_heatmap_matrix,
            title="Actual State (t=$actual_time_idx)",
            aspect_ratio=:equal,
            clims=(global_min_val, global_max_val),
            colorbar=true,
            yflip=false,
            framestyle=:box)

        # Plot for the 'Predicted' state
        p_predicted = heatmap(simulated_heatmap_matrix,
            title="Predicted State (t_init+$k_frame)",
            aspect_ratio=:equal,
            clims=(global_min_val, global_max_val),
            colorbar=true,
            yflip=false,
            framestyle=:box)

        # Combine into a single plot for the frame
        frame_plot = plot(p_actual, p_predicted,
            layout=(1, 2), # Arrange side-by-side
            plot_title="State Comparison at Frame $k_frame",
            plot_titlefontsize=12,
            size=(800, 400), # Adjust size as needed
            left_margin=5Plots.mm,
            bottom_margin=5Plots.mm)

        # Save the combined plot as a single frame
        # Using @sprintf to format the frame number with leading zeros
        output_filename = joinpath(OUTPUT_ANIMATION_DIR, @sprintf("frame_%04d.svg", k_frame))
        # For video generation, PNG is often faster. To switch, change the extension to ".png"
        # output_filename = joinpath(OUTPUT_ANIMATION_DIR, @sprintf("frame_%04d.png", k_frame))
        savefig(frame_plot, output_filename)
    end

    println("\nScript finished. $ANIMATION_LENGTH frames saved in '$OUTPUT_ANIMATION_DIR'.")
    println("You can now compile the frames into a video using a tool like FFmpeg:")
    println("Example FFmpeg command: ffmpeg -r 10 -i frame_%04d.svg -c:v libx264 -pix_fmt yuv420p -y animation.mp4")

end

# --- Script Execution ---
main_animation_script()
