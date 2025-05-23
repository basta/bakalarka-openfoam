include("../MHD_control_cooling.jl")
using .MHD_control_cooling

function main()
    centers_field_path = "$(MHD_control_cooling.CASE_DIR)0/C"
    field_template_path = "../data/fieldsTemplate.mustache"
    force = [-10; 10;-10;10;10;10;-10;10] # min first
    force = [-10; 10;-10;10;-10;-10;10;-10] # max first
    force = Float64.(force)
    Fs = [force force]
    run_sim(Fs, [0,50])
end