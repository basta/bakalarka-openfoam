# Julia script to extract surface temperature data from OpenFOAM case files
# Output: A JLD2 file containing surface cell coordinates, temperature matrix, and time values.

using Printf
using DelimitedFiles
using JLD2 # For saving data in JLD2 format
using FileIO # Works with JLD2 for saving

"""
    parse_openfoam_field(filepath::String)

Parses an OpenFOAM field file (scalar or vector).
It automatically detects if the data is scalar or vector based on the first data line.
Returns a tuple: (data_array, is_vector_field, num_points)
data_array: Array of scalars or Array of 3-element vectors.
is_vector_field: Boolean, true if data is vector.
num_points: Integer, number of data points.
"""
function parse_openfoam_field(filepath::String)
    lines = readlines(filepath)
    data_start_index = -1
    num_points = 0

    # Find the line that just contains the number of points
    for (i, line) in enumerate(lines)
        line_stripped = strip(line)
        if occursin(r"^\d+$", line_stripped) # Check if the line is purely a number
            # The next line should be "("
            if i + 1 <= length(lines) && strip(lines[i+1]) == "("
                num_points = parse(Int, line_stripped)
                data_start_index = i + 2 # Data starts after the number and "("
                break
            end
        end
    end

    if data_start_index == -1 || num_points == 0
        error("Could not find data start or number of points in file: $filepath")
    end

    data_lines = lines[data_start_index : data_start_index + num_points - 1]

    if isempty(data_lines)
        error("No data lines found after parsing header in file: $filepath")
    end

    # Determine if it's vector or scalar from the first data line
    first_data_line_stripped = strip(data_lines[1])
    is_vector_field = false
    if occursin(r"^\(.*\)$", first_data_line_stripped) # Looks like (x y z)
        is_vector_field = true
        # Remove parentheses and split
        parsed_data = [parse.(Float64, split(replace(strip(line), r"[\(\)]" => ""))) for line in data_lines]
        # Ensure all are 3-element vectors
        if any(length(v) != 3 for v in parsed_data)
            error("Inconsistent vector dimensions in file: $filepath")
        end
        data_array = [collect(p) for p in parsed_data] # Convert to Vector{Vector{Float64}}
    else # Scalar field
        is_vector_field = false
        parsed_data = [parse(Float64, strip(line)) for line in data_lines]
        data_array = collect(parsed_data) # Convert to Vector{Float64}
    end

    # Final check for closing parenthesis
    if data_start_index + num_points > length(lines) || strip(lines[data_start_index + num_points]) != ")"
        @warn "Expected closing parenthesis ')' not found or at unexpected position in file: $filepath. Data parsing might be incomplete."
    end

    return data_array, is_vector_field, num_points
end

"""
    extract_surface_temperatures(case_dir::String)

Main function to extract surface temperatures.
case_dir: Path to the OpenFOAM case directory.
Returns a tuple:
  (surface_temp_matrix, surface_cell_coords, time_directory_names)
  surface_temp_matrix: Matrix of surface temperatures (rows: surface cells, cols: time snapshots)
  surface_cell_coords: Vector of surface cell coordinates.
  time_directory_names: Sorted list of time directory names (as strings).
"""
function extract_surface_temperatures(case_dir::String)
    println("Processing OpenFOAM case directory: $case_dir")

    # 1. Parse cell centers from 0/C
    c_file_path = joinpath(case_dir, "0", "C")
    if !isfile(c_file_path)
        error("Cell centers file not found: $c_file_path")
    end

    println("Parsing cell centers from: $c_file_path")
    cell_centers_data, is_vector, num_cells = parse_openfoam_field(c_file_path)

    if !is_vector
        error("Cell centers file $c_file_path does not seem to contain vector data.")
    end
    cell_centers = cell_centers_data

    println("Successfully parsed $num_cells cell centers.")

    # 2. Determine maximum z-coordinate and identify surface cells
    if isempty(cell_centers)
        error("No cell centers found.")
    end

    z_coordinates = [coords[3] for coords in cell_centers]
    if isempty(z_coordinates)
        error("Could not extract z-coordinates from cell centers.")
    end

    max_z = maximum(z_coordinates)
    println("Maximum z-coordinate (surface level): $max_z")

    epsilon = 1e-9
    surface_cell_indices = findall(z -> abs(z - max_z) < epsilon, z_coordinates)

    if isempty(surface_cell_indices)
        @warn "No surface cells found at z = $max_z. Check `epsilon` or data. Using all cells with z == max_z."
        surface_cell_indices = findall(z -> z == max_z, z_coordinates)
        if isempty(surface_cell_indices)
            error("Strict equality also found no surface cells. Problem with z-coordinates or max_z calculation.")
        end
    end
    num_surface_cells = length(surface_cell_indices)
    println("Found $num_surface_cells surface cells.")

    surface_cell_coords = cell_centers[surface_cell_indices]

    # 3. Get list of time-step directories (numeric, excluding "0")
    time_dirs_str = []
    for item in readdir(case_dir)
        item_path = joinpath(case_dir, item)
        if isdir(item_path)
            try
                time_val = parse(Float64, item) # Check if item is a number
                if item != "0"
                    push!(time_dirs_str, item)
                end
            catch e
                # Not a numeric directory, ignore
            end
        end
    end

    sort!(time_dirs_str, by=x->parse(Float64,x))

    if isempty(time_dirs_str)
        @warn "No time-step directories found (excluding '0')."
        # Return empty matrix and the time_dirs_str (which is empty)
        return Matrix{Float64}(undef, num_surface_cells, 0), surface_cell_coords, time_dirs_str
    end
    println("Found time-step directories: $time_dirs_str")

    # 4. Initialize matrix for surface temperatures
    num_snapshots = length(time_dirs_str)
    surface_temp_matrix = Matrix{Float64}(undef, num_surface_cells, num_snapshots)

    # 5. Loop through each time-step directory and extract temperatures
    for (snap_idx, time_dir_name) in enumerate(time_dirs_str)
        t_file_path = joinpath(case_dir, time_dir_name, "T")
        println("Processing snapshot $(snap_idx)/$(num_snapshots): $t_file_path")

        if !isfile(t_file_path)
            @warn "Temperature file not found: $t_file_path. Skipping this snapshot."
            surface_temp_matrix[:, snap_idx] .= NaN
            continue
        end

        temp_data, is_temp_vector, num_temp_points = parse_openfoam_field(t_file_path)

        if is_temp_vector
            @warn "Temperature file $t_file_path appears to be a vector field. Expected scalar."
            if eltype(temp_data) <: AbstractVector && !isempty(temp_data) && eltype(first(temp_data)) <: Number
                 @warn "Taking the first component of the vector temperature field."
                 temp_values = [v[1] for v in temp_data]
            else
                error("Temperature field $t_file_path is vector but not in expected format of list of vectors of numbers.")
            end
        else
            temp_values = temp_data
        end

        if num_temp_points != num_cells
            @warn "Mismatch in number of points between C ($num_cells) and T ($num_temp_points) in $t_file_path."
            if maximum(surface_cell_indices) > num_temp_points
                error("Surface cell indices are out of bounds for temperature data in $t_file_path.")
            end
        end

        surface_temps_for_snapshot = temp_values[surface_cell_indices]
        surface_temp_matrix[:, snap_idx] = surface_temps_for_snapshot
    end

    println("Finished processing all snapshots.")
    return surface_temp_matrix, surface_cell_coords, time_dirs_str
end

# --- Main Execution Logic ---
function main()
    actual_case_dir = ""

    if length(ARGS) == 1 && isdir(ARGS[1])
        actual_case_dir = ARGS[1]
        println("Using case directory from command line argument: $actual_case_dir")
    elseif length(ARGS) == 1 && (ARGS[1] == "--help" || ARGS[1] == "-h")
        println("Usage: julia script_name.jl <path_to_openfoam_case_directory>")
        println("Output: Saves 'openfoam_surface_data.jld2' in the current directory.")
        return
    elseif length(ARGS) > 0
        println("Error: Invalid argument or too many arguments.")
        println("Usage: julia script_name.jl <path_to_openfoam_case_directory>")
        return
    else
        println("Error: No OpenFOAM case directory provided.")
        println("Please provide the case directory as a command line argument:")
        println("  julia your_script_name.jl /path/to/your/openfoam_case")
        return
    end

    try
        surface_temp_matrix, surface_cell_coords, time_dirs_string_array = extract_surface_temperatures(actual_case_dir)

        if !isempty(surface_temp_matrix)
            println("\nSurface Temperature Matrix (Cells x Snapshots):")
            if size(surface_temp_matrix, 1) > 10 || size(surface_temp_matrix, 2) > 5
                println("Matrix dimensions: ", size(surface_temp_matrix))
                println("Showing first 10 rows and 5 columns if available:")
                max_rows_to_show = min(10, size(surface_temp_matrix, 1))
                max_cols_to_show = min(5, size(surface_temp_matrix, 2))
                display(surface_temp_matrix[1:max_rows_to_show, 1:max_cols_to_show])
            else
                display(surface_temp_matrix)
            end

            println("\nSurface Cell Coordinates (corresponding to rows of the matrix):")
            if length(surface_cell_coords) > 10
                 println("Showing first 10 surface cell coordinates:")
                 for i in 1:min(10, length(surface_cell_coords))
                    @printf "Cell %d: (%.6f, %.6f, %.6f)\n" i surface_cell_coords[i][1] surface_cell_coords[i][2] surface_cell_coords[i][3]
                 end
            else
                for (i, coords) in enumerate(surface_cell_coords)
                    @printf "Cell %d: (%.6f, %.6f, %.6f)\n" i coords[1] coords[2] coords[3]
                end
            end

            # Prepare time_values array
            time_values = Float64[]
            if !isempty(time_dirs_string_array)
                time_values = [parse(Float64, t_str) for t_str in time_dirs_string_array]
                println("\nTime values for snapshots: ", time_values)
            else
                println("\nNo time directories processed, time_values array will be empty.")
            end

            # Save data to JLD2 file
            jld2_output_filename = "openfoam_surface_data.jld2"
            println("\nSaving data to JLD2 file: $jld2_output_filename")
            FileIO.save(jld2_output_filename,
                "surface_temperature_matrix", surface_temp_matrix,
                "surface_cell_coordinates", surface_cell_coords,
                "time_values", time_values)
            println("Data successfully saved to $jld2_output_filename")

        else
            println("No surface temperature data was extracted (matrix is empty). JLD2 file not created.")
        end

    catch e
        println(stderr, "An error occurred during script execution: ")
        showerror(stderr, e)
        println(stderr)
    end
end

# Run the main function if the script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
