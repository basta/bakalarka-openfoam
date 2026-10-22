# Julia script to sample surface temperature data onto a uniform grid
# using NearestNeighbors.jl and incorporate experimental inputs from a CSV file.
# Input:
#   1. JLD2 file from 'openfoam_julia_extractor.jl' (temperatures, original coords, sim times)
#   2. CSV file with experimental inputs (Time, u1, u2, ..., u8)
# Output: JLD2 file with sampled temperatures, grid coordinates, simulation times,
#         and an aligned experimental input matrix.

using Printf
using JLD2
using FileIO
using NearestNeighbors
using LinearAlgebra # For LinRange, norm, etc. if needed, though LinRange is Base.
using CSV
using DataFrames

"""
    create_uniform_grid_points(surface_coords::Vector{Vector{Float64}}, N::Int, M::Int)

Creates a list of 3D query points for a uniform N x M grid based on the
extent of the provided surface_coords.
All query points will share the same z-coordinate as the input surface_coords.
"""
function create_uniform_grid_points(surface_coords::Vector{Vector{Float64}}, N::Int, M::Int)
    if isempty(surface_coords)
        error("Surface coordinates are empty. Cannot determine grid extent.")
    end

    x_coords = [pt[1] for pt in surface_coords]
    y_coords = [pt[2] for pt in surface_coords]
    surface_z = surface_coords[1][3]

    min_x, max_x = minimum(x_coords), maximum(x_coords)
    min_y, max_y = minimum(y_coords), maximum(y_coords)

    println("Grid extent: X=[$min_x, $max_x], Y=[$min_y, $max_y], Z=$surface_z")
    println("Creating $N x $M uniform grid.")

    x_grid_range = (N == 1 && min_x == max_x) ? [min_x] : LinRange(min_x, max_x, N)
    y_grid_range = (M == 1 && min_y == max_y) ? [min_y] : LinRange(min_y, max_y, M)

    num_grid_points = N * M
    grid_query_points = Vector{Vector{Float64}}(undef, num_grid_points)

    idx = 1
    for y_val in y_grid_range
        for x_val in x_grid_range
            grid_query_points[idx] = [x_val, y_val, surface_z]
            idx += 1
        end
    end

    return grid_query_points
end

"""
    process_experimental_inputs(csv_filepath::String, simulation_time_values::Vector{Float64})

Reads experimental inputs from a CSV file and aligns them to simulation snapshot times.
Returns a matrix of inputs (num_inputs x num_snapshots).
"""
function process_experimental_inputs(csv_filepath::String, simulation_time_values::Vector{Float64})
    println("Processing experimental inputs from CSV: $csv_filepath")
    if !isfile(csv_filepath)
        error("CSV input file not found: $csv_filepath")
    end

    df = CSV.read(csv_filepath, DataFrame)
    if isempty(df)
        error("CSV file is empty: $csv_filepath")
    end

    # Assuming the first column is 'Time' and subsequent are inputs u1, u2, ...
    if !("Time" in names(df))
        error("CSV file must contain a 'Time' column.")
    end

    csv_times = df.Time
    # Identify input columns (e.g., u1, u2, ..., u8)
    # We assume they are all columns other than 'Time' for simplicity, or you can specify them.
    input_column_names = [name for name in names(df) if name != "Time" && occursin(r"^u\d+$", name)]
    if isempty(input_column_names)
        error("No input columns (e.g., 'u1', 'u2', ...) found in CSV file: $csv_filepath")
    end
    println("Found input columns in CSV: $input_column_names")
    num_inputs = length(input_column_names)

    num_simulation_snapshots = length(simulation_time_values)
    aligned_input_matrix = Matrix{Float64}(undef, num_inputs, num_simulation_snapshots)

    if isempty(simulation_time_values)
        @warn "No simulation time values provided. Returning empty input matrix."
        return Matrix{Float64}(undef, num_inputs, 0)
    end
    if isempty(csv_times)
        @warn "No time values in CSV. Cannot align. Returning matrix of NaNs."
        aligned_input_matrix .= NaN
        return aligned_input_matrix
    end

    for (snap_idx, sim_time) in enumerate(simulation_time_values)
        # Find the index in csv_times that is closest to sim_time
        time_diffs = abs.(csv_times .- sim_time)
        closest_csv_idx = argmin(time_diffs)

        # Extract input values for that closest time
        for (input_idx, col_name) in enumerate(input_column_names)
            aligned_input_matrix[input_idx, snap_idx] = df[closest_csv_idx, col_name]
        end
    end

    println("Successfully processed and aligned $num_inputs experimental inputs for $num_simulation_snapshots snapshots.")
    return aligned_input_matrix
end


"""
    main_sampler(jld2_input_path::String, csv_input_path::String, output_jld2_path::String, N::Int, M::Int)

Main function to load data, create grid, perform nearest neighbor search,
sample temperatures, process CSV inputs, and save results.
"""
function main_sampler(jld2_input_path::String, csv_input_path::String, output_jld2_path::String, N::Int, M::Int)
    println("Loading simulation data from: $jld2_input_path")
    if !isfile(jld2_input_path)
        error("Input JLD2 file not found: $jld2_input_path")
    end

    loaded_data = FileIO.load(jld2_input_path)
    surface_cell_coords = loaded_data["surface_cell_coordinates"]
    surface_temp_matrix = loaded_data["surface_temperature_matrix"]
    time_values = loaded_data["time_values"] # These are simulation snapshot times

    if isempty(surface_cell_coords)
        error("Loaded 'surface_cell_coordinates' is empty. Cannot proceed.")
    end
    if isempty(surface_temp_matrix)
        @warn "Loaded 'surface_temperature_matrix' is empty."
    end
    num_original_surface_cells = size(surface_temp_matrix, 1)
    num_snapshots = size(surface_temp_matrix, 2)

    println("Loaded $(length(surface_cell_coords)) original surface cell coordinates.")
    println("Loaded temperature matrix of size: $num_original_surface_cells x $num_snapshots.")
    println("Loaded $(length(time_values)) simulation time values.")

    # Process experimental inputs from CSV
    experimental_input_matrix = process_experimental_inputs(csv_input_path, time_values)
    # experimental_input_matrix will be num_inputs x num_snapshots

    # Create uniform grid points for querying
    uniform_grid_query_points_vec = create_uniform_grid_points(surface_cell_coords, N, M)
    num_uniform_grid_points = length(uniform_grid_query_points_vec)

    if isempty(surface_cell_coords)
        error("Cannot build KDTree: surface_cell_coords is empty.")
    end
    original_points_matrix = hcat(surface_cell_coords...)

    println("Building KDTree for original $(size(original_points_matrix, 2)) surface points...")
    kdtree = KDTree(original_points_matrix)

    if isempty(uniform_grid_query_points_vec)
        error("Cannot perform KNN search: uniform_grid_query_points_vec is empty.")
    end
    query_points_matrix = hcat(uniform_grid_query_points_vec...)

    println("Performing kNN search (k=1) for $num_uniform_grid_points grid points...")
    idxs_knn, dists_knn = knn(kdtree, query_points_matrix, 1, true)
    nearest_neighbor_indices = [idx_vec[1] for idx_vec in idxs_knn]
    println("kNN search complete.")

    # Sample temperatures onto the uniform grid
    sampled_temps_on_grid = Matrix{Float64}(undef, num_uniform_grid_points, num_snapshots)
    if num_original_surface_cells > 0 && num_snapshots > 0
        println("Sampling temperatures onto the $N x $M grid...")
        for i in 1:num_uniform_grid_points
            original_cell_idx = nearest_neighbor_indices[i]
            if 1 <= original_cell_idx <= num_original_surface_cells
                sampled_temps_on_grid[i, :] = surface_temp_matrix[original_cell_idx, :]
            else
                @warn "Invalid nearest neighbor index $original_cell_idx for grid point $i. Max index is $num_original_surface_cells. Filling with NaN."
                sampled_temps_on_grid[i, :] .= NaN
            end
        end
        println("Temperature sampling complete.")
    elseif num_snapshots == 0
         println("No snapshots in original data. Sampled temperature matrix will be empty along time dimension.")
    else
        println("No original surface cells with temperature data. Sampled temperature matrix will be empty.")
    end

    # Save the results
    println("\nSaving all processed data to JLD2 file: $output_jld2_path")
    FileIO.save(output_jld2_path,
        "uniform_grid_coordinates", uniform_grid_query_points_vec,
        "sampled_temperature_matrix_on_grid", sampled_temps_on_grid,
        "simulation_time_values", time_values,
        "experimental_input_matrix", experimental_input_matrix, # num_inputs x num_snapshots
        "grid_dimensions_NxM", (N,M)
    )
    println("Data successfully saved to $output_jld2_path")

    println("\nSummary:")
    println("  Input JLD2 (Sim Data): $jld2_input_path")
    println("  Input CSV (Exp Data): $csv_input_path")
    println("  Output JLD2: $output_jld2_path")
    println("  Uniform Grid: $N x $M points")
    println("  Sampled Temperature Matrix Dimensions: $(size(sampled_temps_on_grid))")
    println("  Experimental Input Matrix Dimensions: $(size(experimental_input_matrix))")
    println("  Uniform Grid Coordinates Stored: $(length(uniform_grid_query_points_vec)) points")
end

# --- Main Execution Logic ---
function main()
    N_grid_default = 10
    M_grid_default = 10

    if length(ARGS) < 2 || length(ARGS) > 4
        println("Usage: julia script_name.jl <input_jld2_file> <input_csv_file> [N_grid M_grid]")
        println("Example: julia $(@__FILE__) openfoam_surface_data.jld2 experiment_inputs.csv 20 20")
        println("If N_grid and M_grid are not provided, defaults (N=$N_grid_default, M=$M_grid_default) will be used.")
        return
    end

    jld2_input_file = ARGS[1]
    csv_input_file = ARGS[2]

    N_grid = N_grid_default
    M_grid = M_grid_default

    if length(ARGS) >= 3
        try
            N_grid = parse(Int, ARGS[3])
        catch e
            println("Error parsing N_grid value. Using default $N_grid_default. Error: $e")
        end
    end
    if length(ARGS) == 4
        try
            M_grid = parse(Int, ARGS[4])
        catch e
            println("Error parsing M_grid value. Using default $M_grid_default. Error: $e")
        end
    end

    base_jld2_name = first(splitext(basename(jld2_input_file)))
    output_jld2_file = "$(base_jld2_name)_processed_N$(N_grid)xM$(M_grid).jld2"

    println("--- Uniform Grid Sampler & CSV Input Processor ---")
    println("Input JLD2 (Sim): $jld2_input_file")
    println("Input CSV (Exp): $csv_input_file")
    println("Output JLD2: $output_jld2_file")
    println("Grid Dimensions (N x M): $N_grid x $M_grid")

    try
        main_sampler(jld2_input_file, csv_input_file, output_jld2_file, N_grid, M_grid)
    catch e
        println(stderr, "\nAn error occurred during script execution: ")
        showerror(stderr, e)
        println(stderr)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
