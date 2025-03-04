using Revise, JSON, JLD2
#include("../MHD_control_cooling.jl")
include("../openfoam/field_utils.jl")
include("../openfoam/real_field.jl")
includet("../MHD_control_cooling.jl")

using .MHD_control_cooling

EXAMPLE_2D_SAMPLE_POINTS = stack([
    [x;y;0.005] for x in 0.025:0.025:0.1 for y in 0.025:0.025:0.1
])

function isnumeric(str)
    try
        parse(Float64, str)
        return true
    catch
        return false
    end
end


function T_out_sampler(tree, cells, scalar_field, sample_points, case_path)
    critreria_indices = get_criteria_cell_indices(joinpath(case_path, "centers.csv"))
    critreria_indices .+= 1
    criterium = sum(scalar_field[critreria_indices])/size(scalar_field[critreria_indices], 1)
    return [interpolate_tree(tree, scalar_field, sample_points) criterium]
end

function create_u_Y_matrix_for_case(case_path, sample_points) 

    inputs = JSON.parse(open(joinpath(case_path, "sim_info.json")))["inputs"]


    time_dirs = filter(
        x -> isdir(joinpath(case_path, x)) && isnumeric(x) && parse(Int, x) > 0,
    readdir(case_path))

    cells = read_field_vector(joinpath(case_path, "0/C"))
    cell_tree = KDTree(cells)

    sort!(time_dirs, by = x -> parse(Int, x))
    Y = zeros(size(sample_points, 2)+1, length(time_dirs))
    for time_dir in time_dirs
        t = parse(Int, time_dir)
        time_path = joinpath(case_path, time_dir)
        T = read_field_scalar(joinpath(time_path,
        "T"))
        Y[:, t] = T_out_sampler(cell_tree, cells, T, sample_points, case_path)
    end
    return inputs, Y
end

function main()
    dataset = [create_u_Y_matrix_for_case(joinpath("./data/cases/2dexample-1", case), EXAMPLE_2D_SAMPLE_POINTS) for case in readdir("./data/cases/2dexample-1")]
    jldsave("./data/dataset.jld2"; dataset=dataset)
    dataset
end

