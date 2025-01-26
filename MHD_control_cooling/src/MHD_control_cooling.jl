using Revise

includet("./openfoam/field_utils.jl") # goddamn julia development, stupid mistake probably

using LinearAlgebra, BlackBoxOptim, Dates

module MHD_control_cooling
using Revise, JSON, LinearAlgebra

include("./openfoam/field_utils.jl")

function evaluate_criterium(time_path::String, cell_indexes::Vector{Int64})::Float64
    # read all temperatures
    T = read_field_scalar(joinpath(time_path, "T"))
    # extract relevant fields
    T_slice = T[cell_indexes]
    # calculate criterium
    return sum(T_slice)/length(T_slice) 
end

function run_sim(F::Function)
    @info "Running simulation"
    centers_field_path = "../2d-example/0/C"
    field_template_path = "./data/fieldsTemplate.mustache"
    field_string = create_force_field_string(F, centers_field_path, field_template_path, true)
    set_field_at_time("../2d-example", "0", field_string, "F")
    run_case("../2d-example", "./icoHeatExternalForce")
end

function get_criterium_in_time(case_path::String)
    criteria = [];
    slice_idxs = get_slice_xmin(read_field_vector(joinpath(case_path, "0/C")));
    for i in 1:1000
        try
            time_name = string(i)
            time_path = joinpath(case_path, time_name)
            J = evaluate_criterium(time_path, slice_idxs)
            push!(criteria, evaluate_criterium(time_path, slice_idxs))
        catch LoadError
            break
        end
    end
    return criteria
end

function log_experiment_json(dir, name, criteria, desc)
    json = JSON.json(Dict("name"=>name, "criteria"=>criteria, "description"=>desc))
    open(joinpath(dir, "$name.json"), "w") do f
        write(f, json)
    end
end

function main()
    centers_field_path = "../../2d-example/0/C"
    field_template_path = "../data/fieldsTemplate.mustache"
    F = x -> [x[1], x[2], x[3]]
    field_string = create_force_field_string(F, centers_field_path, field_template_path)
    set_field_at_time("../../2d-example", "0", field_string, "F")
    run_case("../../2d-example", "./icoHeatExternalForce")
    println(field_string)
end

end # module MHD_control_cooling

function baseline_force(x)
    return [0.0, 0.0, 0.0]
end

function rotation_force(x) 
    println("forcefn: $(typeof(x)) $x\t");
    return cross(x, [0.5, 0.5, 0.5])*50
end

function central_force(x)
    return x .- [0.5, 0.5, 0.5]
end

function parametric_cells_generator(θ, xCells, yCells, dim=3)::Function
    @assert xCells * yCells * dim == length(θ)
    function F(X)::Vector{Float64}
        xPos = ceil(X[1]*xCells+0.01)
        yPos = ceil(X[2]*yCells+0.01)
        xPos = clamp(xPos, 1, xCells)
        yPos = clamp(yPos, 1, yCells)


        idx = Int(dim*((xPos-1)*yCells + (yPos-1)) + 1)
        return θ[idx:idx+(dim-1)]
    end
end

quartal_force = parametric_cells_generator([0.04045899007671505, 0.46618114797333665, -0.0009118347788213965, 0.9513413554885509, -0.5429034604752823, 0.8688778971762849, -0.9053939364779692, 0.8857732272567411, 0.31375231902429757, 0.8100247853948019, 0.1093953582080961, 0.5896836503093306, -0.7568366612862325, 0.28160218122185054, -0.9437382625074121, 0.5371512919099412, -0.31551113665600206, 0.8563497146424621, 0.9721859245051305, 0.754209782250431, 0.8677320742981145, 0.8451275884265578, 0.13708732702575843, 0.6131597814970461, -0.04357261653343791, -0.7211213086112536, -0.4732729848866935, 0.5568502379230311, -0.4908464947540193, -0.7370713733435196, 0.2402778390005072, -0.7353868860035349, 0.6945452160383623, -0.3596922091559357, 0.8924828776906942, -0.21973258257016784, -0.6819051740630332, 0.7090368781262668, 0.00818920092892711, -0.48746312068004777, -0.8707585436420283, -0.11206378977506903, 0.9739015828332115, -0.3600719510441425, 0.976073658012101, 0.46613050746984086, 0.8277529174039324, -0.4029413845152499, -0.5915814895174086, 0.9454767050706295, -0.2997066968307025, -0.09128217812603245, -0.28886350420505824, -0.741355776615034, -0.5947942639447024, 0.8196978787420647, 0.40109735507003275, 0.5974378758638021, 0.834443971236232, 0.7518661739829385, -0.4718585128664097, -0.9336703697748423, -0.5566139051954919, -0.9072687619813506, -0.844042694076521, -0.9130447378324676, -0.9738601577920518, -0.24471801742817204, 0.42779733698549965, -0.4545036095952283, 0.5689830501075652, -0.7382464768616737], 6, 6, 2)


function main()
    MHD_control_cooling.run_sim(x -> vcat(quartal_force(x), [0.0]).*50)
    criteria = MHD_control_cooling.get_criterium_in_time("../2d-example")
    MHD_control_cooling.log_experiment_json("./experiments", "quartal", criteria, "quartal")
end

function blackbox(θ)::Float64
    quartal_force = parametric_cells_generator(θ, 6, 6, 2)
    MHD_control_cooling.run_sim(x -> vcat(quartal_force(x), [0.0]).*50)
    criteria = MHD_control_cooling.get_criterium_in_time("../2d-example")
    timestamp = Dates.format(Dates.now(), "yyyy-mm-dd_HH-MM-SS")
    name = "blackbox_$(criteria[end])_$(timestamp)"
    desc = string(θ)
    MHD_control_cooling.log_experiment_json("./experiments/blackbox_6x6", name, criteria, desc)
    return criteria[end]
end

function optimize()
    res = bboptimize(blackbox, SearchRange = (-1.0, 1.0), NumDimensions = 6*6*2)
    println("Optimalization result is")
    println(res)
    return res
end

main()