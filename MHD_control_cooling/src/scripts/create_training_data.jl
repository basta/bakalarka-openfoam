using Revise, JSON, JLD2, ProgressLogging, DataFrames, CSV
#include("../MHD_control_cooling.jl")
include("../openfoam/field_utils.jl")
include("../openfoam/real_field.jl")
includet("../MHD_control_cooling.jl")

using .MHD_control_cooling

EXAMPLE_2D_SAMPLE_POINTS = stack([
    [x;y;0.005] for x in 0.025:0.025:0.1 for y in 0.025:0.025:0.1
])




function T_out_sampler(tree, cells, scalar_field, sample_points, case_path)
    critreria_indices = get_criteria_cell_indices(joinpath(case_path, "centers.csv"))
    critreria_indices .+= 1
    criterium = sum(scalar_field[critreria_indices])/size(scalar_field[critreria_indices], 1)
    return [interpolate_tree(tree, scalar_field, sample_points) criterium]
end

function find_largest_smaller(arr::Vector{Int}, num::Int)
    idx = searchsortedlast(arr, num)
    idx > 0 && arr[idx] >= num && (idx -= 1)
    return idx
end

function extract_inputs_from_file(file_path::String)
    # Read the file content
    content = read(file_path, String)
    
    # Find the INPUTS section
    inputs_start = findfirst("// INPUTS:", content)
    if inputs_start === nothing
        return Float64[] # Return empty array if INPUTS section not found
    end
    
    # Find the end of the INPUTS section (either "// Test" or "dimensions")
    inputs_end = findfirst("// Test", content)
    if inputs_end === nothing
        inputs_end = findfirst("dimensions", content)
    end
    
    if inputs_end === nothing
        return Float64[] # Return empty array if end marker not found
    end
    
    # Extract the text between INPUTS: and the end marker
    inputs_section = content[inputs_start[end]+1:inputs_end[1]-1]
    
    # Extract numbers using regex
    numbers = []
    for line in split(inputs_section, '\n')
        # Remove comments and trim whitespace
        line = strip(replace(line, r"//\s*" => ""))
        # Skip empty lines
        if isempty(line)
            continue
        end
        # Try to parse the number
        try
            push!(numbers, parse(Float64, line))
        catch
            # Skip lines that don't contain valid numbers
        end
    end
    
    return numbers
end


function create_u_Y_matrix_for_case(case_path, sample_points; inputs_file=nothing, csv_file=nothing)

    if !isnothing(inputs_file)
        datas = jldopen(inputs_file)
    end

    if !isnothing(csv_file)
        df = CSV.read(csv_file, DataFrame)
        df[:, :Time] = Int32.(round.(df[:, :Time]))
    end



    time_dirs = filter(
        x -> isdir(joinpath(case_path, x)) && isnumeric(x) && parse(Float64, x) > 0,
    readdir(case_path))

    cells = read_field_vector(joinpath(case_path, "0/C"))
    cell_tree = KDTree(cells)

    sort!(time_dirs, by = x -> parse(Float64, x))
    Y = zeros(size(sample_points, 2)+1, length(time_dirs))
    times = []
    inputs = zeros(8, length(time_dirs))
    last_input = zeros(8)
    @progress for (i, time_dir) in enumerate(time_dirs)
        t = parse(Float64, time_dir)
        push!(times, t)
        time_path = joinpath(case_path, time_dir)
        T = read_field_scalar(joinpath(time_path,
        "T"))
        @info time_path
        if !isnothing(csv_file)
            idx = findlast(df[:, :Time] .<= t)
            if idx === nothing
                idx = 1
            end 

            input_t = Vector(df[idx, 2:end])
        else 
            input_t = extract_inputs_from_file(joinpath(time_path, "F"))
        end
        
        if !isempty(input_t)
            last_input = input_t
        else 
            input_t = last_input
        end
        inputs[:, i] = input_t

        Y[:, i] = T_out_sampler(cell_tree, cells, T, sample_points, case_path)
    end
    if !isnothing(inputs_file)
        inputs_in = datas["inputs"]
        F_times = datas["F_times"]
        inputs = zeros(size(inputs_in[1], 1), length(times))
        for t in times
            idx = find_largest_smaller(F_times, t-2) + 1 #TODO the jld was created wrong
            if (t < 100)
                @info "Time:$t idx:$idx"
            end
            inputs[:, t] = inputs_in[idx]
        end
    end
    return inputs, Y
end

function main()
    dataset_out = "./data/dataset-long.jld2"
    dataset = [create_u_Y_matrix_for_case(joinpath("./data/cases/2dexample-dynamic",
     case),
     EXAMPLE_2D_SAMPLE_POINTS) for case in readdir("./data/cases/2dexample-dynamic")]
    jldsave(dataset_out; dataset=dataset)
    println("Dataset saved: $dataset_out")
    return  dataset
end

function main_for_long()
    dataset_out = "./data/dataset-long.jld2"
    dataset = create_u_Y_matrix_for_case(("./data/cases/long2d-2"), 
    EXAMPLE_2D_SAMPLE_POINTS)
    jldsave(dataset_out; dataset=dataset)
    println("Dataset saved: $dataset_out")
    
    return dataset
end



function main_for_live(csv_path)
    dataset_out = "./data/dataset-live2.jld2"
    dataset = create_u_Y_matrix_for_case(("./data/cases/2d-example-live2"), 
    EXAMPLE_2D_SAMPLE_POINTS, csv_file=csv_path)
    jldsave(dataset_out; dataset=dataset)
    println("Dataset saved: $dataset_out")
    
    return dataset
end