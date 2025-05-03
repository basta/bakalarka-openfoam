#=
This script processes experimental data from a single specified run directory.
It reads control input data and system state data from CSV files within that directory,
aligns them based on timestamps, and extracts the relevant data.
Finally, it saves the feature data (X_data) and control input data (U_data)
into a JLD2 file named 'data.jld2' in the current working directory.

Usage:
julia process_single_dir.jl /path/to/your/experiment_runX

Arguments:
- The first command-line argument should be the full or relative path to the
  experiment directory containing the *_control.csv and *_control.csv files.
=#

using CSV
using DataFrames
using JLD2
using Glob
using LinearAlgebra # For transpose

println("Starting data processing for a single directory...")

# --- Check for Command-Line Argument ---
if isempty(ARGS)
    @error "Usage: julia process_single_dir.jl <path_to_experiment_directory>"
    exit(1)
end

experiment_dir = ARGS[1]

if !isdir(experiment_dir)
    @error "Provided path is not a valid directory: $experiment_dir"
    exit(1)
end

println("Processing directory: ", experiment_dir)

# --- Find CSV Files ---
control_files = glob("*_control.csv", experiment_dir)
data_files = glob("*_state.csv", experiment_dir)

# Basic error checking
if isempty(control_files)
    @error "No *_control.csv file found in $experiment_dir. Aborting."
    exit(1)
elseif length(control_files) > 1
    @warn "Multiple *_control.csv files found in $experiment_dir. Using first one: $(basename(control_files[1]))."
end
if isempty(data_files)
    @error "No *_control.csv file found in $experiment_dir. Aborting."
    exit(1)
elseif length(data_files) > 1
    @warn "Multiple *_control.csv files found in $experiment_dir. Using first one: $(basename(data_files[1]))."
end

control_file = control_files[1]
data_file = data_files[1]
println("  Control file: ", basename(control_file))
println("  Data file:    ", basename(data_file))

# --- Process Data ---
local X_data, U_data # Ensure variables are accessible outside try block

try
    df_control = CSV.read(control_file, DataFrame)
    df_data = CSV.read(data_file, DataFrame)

    # --- Identify Columns ---
    # Identify control columns (starting with 'u' followed by a digit)
    control_cols = filter(name -> startswith(string(name), 'u') && isdigit(string(name)[2]), names(df_control))
    if isempty(control_cols)
        @error "Could not automatically detect control columns (u1, u2, ...) in $(basename(control_file)). Please specify manually if needed. Aborting."
        exit(1)
    end
    println("    Detected Control Columns: ", control_cols)

    # Identify feature columns (all columns in data except 'Time')
    feature_cols = filter(name -> name != :Time, names(df_data))
    if isempty(feature_cols)
        @error "Could not automatically detect feature columns in $(basename(data_file)). Please specify manually if needed. Aborting."
        exit(1)
    end
    println("    Detected Feature Columns: ", feature_cols)

    # --- Data Alignment ---
    df_data_clean = dropmissing(df_data, feature_cols)
    if isempty(df_data_clean)
         @error "No valid (non-NaN) data rows found in $(basename(data_file)) after cleaning. Aborting."
         exit(1)
    end

    aligned_U_rows = []
    aligned_X_rows = []
    last_data_idx = 0
    min_data_time = df_data_clean.Time[1]

    for i in 1:nrow(df_control)
        t_control = df_control.Time[i]
        if t_control < min_data_time
            continue
        end

        current_data_idx = 0
        for j in (last_data_idx + 1):nrow(df_data_clean)
             if df_data_clean.Time[j] <= t_control
                 current_data_idx = j
             else
                 break
             end
             if j == nrow(df_data_clean) && df_data_clean.Time[j] <= t_control
                current_data_idx = j
             end
        end

         if current_data_idx > 0
             push!(aligned_U_rows, Vector(df_control[i, control_cols]))
             push!(aligned_X_rows, Vector(df_data_clean[current_data_idx, feature_cols]))
             last_data_idx = current_data_idx
        end
    end

    if isempty(aligned_U_rows)
        @error "No data points could be aligned for directory $(basename(experiment_dir)). Aborting."
        exit(1)
    end

    # Convert aligned rows to DataFrames
    run_U_df = DataFrame(hcat(aligned_U_rows...)', control_cols)
    run_X_df = DataFrame(hcat(aligned_X_rows...)', feature_cols)

    # Convert DataFrames to Matrices (features x snapshots) and (inputs x snapshots)
    global U_data = Matrix(run_U_df)'
    global X_data = Matrix(run_X_df)'

    println("  Aligned and processed snapshots: ", size(U_data, 2))
    println("  U_data dimensions (inputs x snapshots): ", size(U_data))
    println("  X_data dimensions (features x snapshots): ", size(X_data))

catch e
    @error "Error processing files in directory $experiment_dir: $e"
    exit(1)
end

# --- Save Data ---
output_filename = "data.jld2"
println("\nSaving data to $output_filename in the current directory...")
# Save the data to a JLD2 file in the current working directory
jldsave(output_filename; X_data, U_data)

println("Processing complete. Data saved to $output_filename.")
