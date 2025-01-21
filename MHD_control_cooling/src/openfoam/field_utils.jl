module FieldUtils

using OteraEngine
using LinearAlgebra


export create_force_field_string, set_field_at_time, run_case

function lineToVec(line)
	coords = split(strip(line, ['(', ')']))
	return [parse(Float64, x) for x in coords]
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
    min_coords = minimum.(eachrow(cells))
    max_coords = maximum.(eachrow(cells))
    min_x, min_y, min_z = min_coords
    max_x, max_y, max_z = max_coords

    if normalize
        normalized_cells = [(cell .- [min_x, min_y, min_z]) ./ ([max_x - min_x, max_y - min_y, max_z - min_z]) for cell in cells]
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

function run_case(case_path::String, solver_name::String)
    if !isdir(case_path)
        error("Case path $case_path does not exist.")
    end

    run(setenv(`pwd`, dir=case_path))
    run(setenv(`openfoam2406-run -c $solver_name`, dir=case_path))
    
end


end