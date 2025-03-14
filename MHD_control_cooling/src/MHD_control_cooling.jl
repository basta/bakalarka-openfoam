# using Revise, Glob, Logging, GLMakie
# if !endswith(pwd(), "MHD_control_cooling")
#     cd("/home/basta/Projects/bakalarka-openfoam/MHD_control_cooling")
# end


# includet("./openfoam/field_utils.jl") # goddamn julia development, stupid mistake probably
# includet("./openfoam/real_field.jl")

# using LinearAlgebra, BlackBoxOptim, Dates

module MHD_control_cooling

using Revise, JSON, LinearAlgebra
using Revise, Glob, Logging
using LinearAlgebra, BlackBoxOptim, Dates

logger = ConsoleLogger(stderr, Logging.Info)

global_logger(logger)




include("./openfoam/field_utils.jl")
include("./openfoam/real_field.jl")

elMagData = create_trees("./data/forceFields")



CASE_DIR = "../2d-example/"
TRANS_FUN = example2magman

export CASE_DIR, TRANS_FUN, run_sim, real_force_generator


function evaluate_criterium(time_path::String, cell_indexes::Vector{Int64})::Float64
    # read all temperatures
    T = read_field_scalar(joinpath(time_path, "T"))
    # extract relevant fields
    T_slice = T[cell_indexes]
    # calculate criterium
    return sum(T_slice)/length(T_slice)
end

function run_sim(inputs::AbstractArray{Float64, 2}, F_times::AbstractVector{Int}; ambient_T::Real=20)
    @info "Running simulation"
    centers_field_path = "$(CASE_DIR)0/C"
    field_template_path = "./data/fieldsTemplate.mustache"

    case_end_time = read_control_dict_entry(CASE_DIR, "endTime")
    @info "Running a case with endTime $case_end_time"

    @assert issorted(F_times)

    if F_times[end] > parse(Int, case_end_time)
        @warn "Last time value of F_times is higher than case_end_time"
        push!(F_times,F_times[end])
    else
        push!(F_times, parse(Int, case_end_time))
    end
    @info F_times
    clean_case(CASE_DIR)
    set_control_dict_entry(CASE_DIR, "startFrom", "latestTime")
    for i in 1:size(inputs, 2)
        force = real_force_generator(inputs[:, i])
        field_string = create_force_field_string(force, centers_field_path, field_template_path,inputs[:, i], false)
        @info "Simtime $(F_times[i+1])"
        set_field_at_time("$(CASE_DIR)", get_last_fields_time(CASE_DIR), field_string, "F")
        set_control_dict_entry(CASE_DIR, "endTime", F_times[i+1])
        run_case("$(CASE_DIR)", "./icoHeatExternalForce") # TODO use the other solver 
    end
end

function get_criterium_in_time(case_path::String)
    criteria = [];
    # slice_idxs = get_slice_xmin(read_field_vector(joinpath(case_path, "0/C")));
    slice_idxs = get_criteria_cell_indices(joinpath(case_path, "centers.csv"));
    for i in 5:5:1000
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

function real_force_generator(θ, per_partes_out=false)
    @assert length(θ) == 8
    function real_force(x)
        x_mg = MHD_control_cooling.TRANS_FUN(x)
        return get_force_at_point(elMagData, x_mg, θ, per_partes_out);
    end
end


function main()
    centers_field_path = "$(CASE_DIR)0/C"
    field_template_path = "../data/fieldsTemplate.mustache"
    F = x -> [x[1], x[2], x[3]]
    field_string = create_force_field_string(F, centers_field_path, field_template_path)
    set_field_at_time("../$(CASE_DIR)", "0", field_string, "F")
    run_sim("../$(CASE_DIR)", "./icoHeatExternalForce")
end

end # module MHD_control_cooling
