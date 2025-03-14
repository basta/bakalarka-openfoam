include("../MHD_control_cooling.jl")
using .MHD_control_cooling

function main()
    centers_field_path = "$(MHD_control_cooling.CASE_DIR)0/C"
    field_template_path = "../data/fieldsTemplate.mustache"
    Fs = [real_force_generator(rand(8)), real_force_generator(rand(8))]
    run_sim(Fs, [0, 20])
end