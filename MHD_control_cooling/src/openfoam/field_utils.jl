using OteraEngine
using LinearAlgebra
using Logging


function lineToVec(line)
	coords = split(strip(line, ['(', ')']))
	return [parse(Float64, x) for x in coords]
end

function normalize_cells(cells)
    # Get min and max for each dimension separately
    min_coords = minimum.(zip(cells...))  # Get minimum of each dimension
    max_coords = maximum.(zip(cells...))  # Get maximum of each dimension
    
    # Calculate range for each dimension
    ranges = max_coords .- min_coords
    
    # Avoid division by zero by replacing zero ranges with 1
    # Normalize each cell
    normalized_cells = [
        (cell .- min_coords) ./ ranges
        for cell in cells
    ]
    for (i, cell) in enumerate(normalized_cells)
        for dim in 1:3
            if isnan(cell[dim])
                normalized_cells[i][dim] = 0.5
            end 
        end
    end
    
    return normalized_cells
end

function create_force_field_string(F::Function, centers_field_path::String, field_template_path::String, normalize::Bool=false)
    cellsText = read(centers_field_path, String)
    cellStart = 0
    nCells = 0
    lines = split(cellsText, "\n")
    for i in eachindex(lines)
        if (lines[i] == "(")
            cellStart = i+1
            nCells = parse(Int32, lines[i-1])
            break
        end
    end

    cells = [lineToVec(line) for line in lines[cellStart:cellStart+nCells-1]]
    min_coords = minimum(eachrow(cells))
    max_coords = maximum(eachrow(cells))
    min_x, min_y, min_z = min_coords[1]
    max_x, max_y, max_z = max_coords[1]

    if normalize
        normalized_cells = normalize_cells(cells)

        A = F.(normalized_cells) 

    else 
        A = F.(cells) 
    end

    template = read(field_template_path, String)
	data = Dict(:nCells=>nCells, :A=>A)
    tmp = Template(template, path=false)
    render = Base.invokelatest(tmp, init=data)
    return render
end

function read_field_vector(field_path::String)::Matrix
    cellsText = read(field_path, String)
    cellStart = 0
    nCells = 0
    lines = split(cellsText, "\n")
    for i in eachindex(lines)
        if (lines[i] == "(")
            cellStart = i+1
            nCells = parse(Int32, lines[i-1])
            break
        end
    end
    cells = [lineToVec(line) for line in lines[cellStart:cellStart+nCells-1]]

    cells_matrix = zeros(3,size(cells,1))
    for (i, cell) in enumerate(cells)
        cells_matrix[:, i] = cell
    end
    return cells_matrix
end

function read_field_scalar(field_path::String)::Vector{Float64}
    cellsText = read(field_path, String)
    cellStart = 0
    nCells = 0
    lines = split(cellsText, "\n")
    for i in eachindex(lines)
        if (lines[i] == "(")
            cellStart = i+1
            nCells = parse(Int32, lines[i-1])
            break
        end
    end
    cells = [parse(Float32, line) for line in lines[cellStart:cellStart+nCells-1]]
    cells_matrix = zeros(size(cells,1))
    for (i, cell) in enumerate(cells)
        cells_matrix[i] = cell
    end
    return cells_matrix
end

function get_slice_xmin(cells::Matrix)
    min_coords = minimum.(eachrow(cells))
    tol = 0.0025
    is_xmin = map(cell -> cell[1] < min_coords[1]+tol, eachcol(cells))
    return findall(is_xmin)
end


function set_field_at_time(case_path::String, time::String, field_file_string::String, field_name::String)
    if !isdir(case_path)
        error("Case path $case_path does not exist.")
    end

    time_path = joinpath(case_path, time)
    if !isdir(time_path)
        error("Time path $time_path does not exist.")
    end
    field_path = joinpath(case_path, time, field_name)
    open(field_path, "w+") do f
        write(f, field_file_string)
    end
end

function run_case(case_path::String, solver_name::String, verbose = false)
    if !isdir(case_path)
        error("Case path $case_path does not exist.")
    end

    if verbose
        run(setenv(`pwd`, dir=case_path))
        try
            run(setenv(`openfoam2406-run -c $solver_name`, dir=case_path))
        catch e
            if isa(e, Base.IOError)
                @info "openfoam2406-run failed, trying the other one (For my desktop pc)"
                run(setenv(`openfoam2406 -c $solver_name`, dir=case_path))
            else
                rethrow(e)
            end
        end
    else
        open("output.txt", "w") do file
            redirect_stdout(file) do
                run(setenv(`pwd`, dir=case_path))
                try
                    run(setenv(`openfoam2406-run -c $solver_name`, dir=case_path))
                catch e
                    if isa(e, Base.IOError)
                        @info "openfoam2406-run failed, trying the other one (For my desktop pc)"
                        run(setenv(`openfoam2406 -c $solver_name`, dir=case_path))
                    else
                        rethrow(e)
                    end
                end
            end
        end
    end    
end

