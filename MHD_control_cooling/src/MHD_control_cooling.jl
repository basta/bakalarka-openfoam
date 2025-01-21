includet("./openfoam/field_utils.jl") # goddamn julia development, stupid mistake probably
using .FieldUtils

module MHD_control_cooling
using Revise

include("./openfoam/field_utils.jl")
using .FieldUtils

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
