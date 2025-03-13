using Revise
include("../MHD_control_cooling.jl")
include("../openfoam/real_field.jl")
includet("../MHD_control_cooling.jl")

using .MHD_control_cooling

DATA_ID = "2dexample-dynamic"

struct SimInfo
    inputs::AbstractArray{Float64}
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
    N_changes = Int(round.(rand().*9))+1
    inputs_arr = (rand(8,N_changes).-0.5)*10
    sim_info = SimInfo(inputs_arr, DATA_ID, MHD_control_cooling.CASE_DIR)
    Fs = [MHD_control_cooling.real_force_generator(inputs) for inputs in eachcol(inputs_arr)]
    case_end_time = parse(Int, read_control_dict_entry(MHD_control_cooling.CASE_DIR, "endTime"));
    F_times = collect(1:(case_end_time/N_changes):case_end_time)
    F_times = round.(F_times)
    @debug F_times
    F_times = Int.(F_times)
    run_sim(Fs, F_times);
    save_case_data(MHD_control_cooling.CASE_DIR, save_dir, sim_info)
end

function run_infinite_case()
    MAX_TIME = 100000
    time = 0
    Fs = []
    F_times::Vector{Int} = []
    inputs = []
    while time < MAX_TIME
        seg_len = Int(round(rand()*50))
        time += seg_len
        push!(F_times, time)
        input = (rand(8).-0.5)*10
        push!(Fs, MHD_control_cooling.real_force_generator(input))
        push!(inputs, input)
    end
    jldsave("infinite_case.jld2"; inputs=inputs, F_times=F_times)
    MHD_control_cooling.run_sim(Fs, F_times)
end

function main(ininite=false)
    ctr = 0
    @info "Creating random cases for $DATA_ID"
    @info "Saving to $(pwd())/data/cases/$(DATA_ID)"
    while true
        if ininite
            run_infinite_case()
        else
            run_random_case("./data/cases/$(DATA_ID)")
        end
        ctr += 1
        if ctr % 10 == 0
            @info "Generated $ctr cases"
        end
    end
end
