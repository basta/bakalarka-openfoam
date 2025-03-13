using Revise, JSON, JLD2, ProgressLogging
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


function create_u_Y_matrix_for_case(case_path, sample_points; inputs_file=nothing)

    if !isnothing(inputs_file)
        datas = jldopen(inputs_file)
    else
        inputs = JSON.parse(open(joinpath(case_path, "sim_info.json")))["inputs"]
    end 



    time_dirs = filter(
        x -> isdir(joinpath(case_path, x)) && isnumeric(x) && parse(Int, x) > 0,
    readdir(case_path))

    cells = read_field_vector(joinpath(case_path, "0/C"))
    cell_tree = KDTree(cells)

    sort!(time_dirs, by = x -> parse(Int, x))
    Y = zeros(size(sample_points, 2)+1, length(time_dirs))
    times = []
    @progress for time_dir in time_dirs
        t = parse(Int, time_dir)
        push!(times, t)
        time_path = joinpath(case_path, time_dir)
        T = read_field_scalar(joinpath(time_path,
        "T"))
        Y[:, t] = T_out_sampler(cell_tree, cells, T, sample_points, case_path)
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

function main_for_long(inputs_file)
    dataset_out = "./data/dataset-long.jld2"
    dataset = create_u_Y_matrix_for_case(("./data/cases/long2d"), 
    EXAMPLE_2D_SAMPLE_POINTS; inputs_file=inputs_file)
    jldsave(dataset_out; dataset=dataset)
    println("Dataset saved: $dataset_out")
    
    return dataset
end
