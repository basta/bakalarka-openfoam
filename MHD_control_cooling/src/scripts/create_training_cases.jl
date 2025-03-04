using Revise
include("../MHD_control_cooling.jl")
include("../openfoam/real_field.jl")
includet("../MHD_control_cooling.jl")

using .MHD_control_cooling

DATA_ID = "2dexample-2"

struct SimInfo
    inputs::Vector{Float64}
    data_id::String
    orig_case::String
end

using JSON, FileIO

function save_sim_info(sim_info::SimInfo, file_path::String)
    json_data = JSON.json(Dict(
        "inputs" => sim_info.inputs,
        "data_id" => sim_info.data_id,
        "orig_case" => sim_info.orig_case
        ))
    open(file_path, "w+") do file
        write(file, json_data)
    end
end


function save_case_data(case, data_path, sim_info)
    mkpath(data_path)
    dirs = filter(x -> isdir(joinpath(data_path, x)), readdir(data_path))
    numbers = [parse(Int, dir) for dir in dirs if tryparse(Int, dir) !== nothing]
    new_dir_number = isempty(numbers) ? 1 : maximum(numbers) + 1
    new_dir_path = joinpath(data_path, string(new_dir_number))
    mkpath(new_dir_path)
    cp(case, new_dir_path; force=true)
    save_sim_info(sim_info, joinpath(new_dir_path, "sim_info.json"))
end


function run_random_case(save_dir::String)
    inputs = (rand(8).-0.5)*10
    sim_info = SimInfo(inputs, DATA_ID, MHD_control_cooling.CASE_DIR)
    run_sim(real_force_generator(inputs));
    save_case_data(CASE_DIR, save_dir, sim_info)
end

function main()
    ctr = 0
    @info "Creating random cases for $DATA_ID"
    @info "Saving to $(pwd())/data/cases/$(DATA_ID)"
    while true
        run_random_case("./data/cases/$(DATA_ID)")
        ctr += 1
        if ctr % 10 == 0
            @info "Generated $ctr cases"
        end
    end
end
