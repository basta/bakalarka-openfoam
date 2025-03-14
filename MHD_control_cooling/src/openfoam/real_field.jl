using NearestNeighbors, LinearAlgebra, Logging, PrettyTables, CairoMakie, CSV, DataFrames
using JLD2
include("./field_utils.jl")

mutable struct ElMagneticInput
    poses::Array{Float64, 2}
    vecs::Array{Float64, 2}
    tree::KDTree
    isEl::Bool
end

mutable struct ElMagneticData
    inputs::Vector{ElMagneticInput}
end

function interpolate_tree(nn_tree, vecs, points)
    if ndims(vecs) == 1
        vecs = reshape(vecs, 1, :)
    end
	idxssE, distssE = knn(nn_tree, points[:, :], 5)
	vecsE = []
	for (i, idxs) in enumerate(idxssE)
			dists = distssE[i]
			pos = points[:, i]
			i_vecs = vecs[:, idxs]
			vec3 = zeros(size(i_vecs)[1])
			dist_sum = sum(dists)
			w_sum = 0
			for (dist, i_vec) in zip(dists, eachcol(i_vecs))
				w_i = 1/dist
				w_sum += w_i
				# i_vec = cross(collect(i_vec), [1;0;0])
				vec3 = vec3 .+ i_vec.*w_i
			end

			vec3 = vec3 ./ w_sum

			push!(vecsE,vec3)
		end
	return stack(vecsE)
end


function create_elmaginput(jld_path::String)
    data = load(jld_path)
    poses = data["coords"].parent'
    isEl = haskey(data, "Ez")
    if isEl
        vecs = stack(zip(data["Ex"], data["Ey"], data["Ez"]))
    else
        vecs = stack(zip(data["Bx"], data["By"], data["Bz"]))
    end
    @info "Creating $(isEl ? "el" : "mag" ) tree from $(jld_path)"
    return ElMagneticInput(poses, vecs, KDTree(poses), isEl)
end

function create_trees(jld_dir::String)
    @info "Calculating KDTrees from $(jld_dir)"

    maginputs = []

    for i in 1:100
        if ispath(joinpath(jld_dir, "E$i.jld"))
            push!(maginputs, create_elmaginput(joinpath(jld_dir, "E$i.jld")))
        else
            @debug "Found $(i-1) el inputs, notfound:$(joinpath(jld_dir, "E$i.jld"))"
            break
        end
    end

    for i in 1:100
        if ispath(joinpath(jld_dir, "B$i.jld"))
            push!(maginputs, create_elmaginput(joinpath(jld_dir, "B$i.jld")))
        else
            @debug "Found $(i-1) mag inputs, notfound:$(joinpath(jld_dir, "E$i.jld"))"
            break
        end
    end

    return ElMagneticData(maginputs)
end

function get_force_at_point(elMagData::ElMagneticData, point::Array{Float64, 1}, θ::AbstractVector{<:Number}, per_partes_out = false)
    F = zeros(3)
    Fs = zeros(3, 8*8)
    i = 0
    for (Bi, maginput) in enumerate(elMagData.inputs)
        if maginput.isEl
            continue
        end
        for (Ei, elinput) in enumerate(elMagData.inputs)
            if !elinput.isEl
                continue
            end
            E = interpolate_tree(elinput.tree, elinput.vecs, vec(point))[:, 1] .* θ[Ei]
            B = interpolate_tree(maginput.tree, maginput.vecs, vec(point))[:, 1] .* θ[Bi]
            F += cross(E, B)
            if θ[Ei] ≈ 0 || θ[Bi] ≈ 0
                continue
            end
            i+=1
            Fs[:, i] = F
        end
    end
    # pretty_table(Fs,
    # formatters = (v, i, j) -> round(v, digits=2),
    # tf = tf_borderless
    # )
    if !per_partes_out
        return F
    else
        return Fs
    end
end

function norm2magman(vec::AbstractArray{Float64, 1})
    vec = vec[:]
    vec[1] = vec[1] / 10 - 0.05
    vec[2] = vec[2] / 10 - 0.05
    vec[3] = 0.005
    return vec
end

function H2magman(vec::AbstractArray{Float64, 1})
    source_cube = (x=(0.05, 0.15), y=(0.0, 0.1), z=(0.0, 1.0))
    target_cube = (x=(-0.05, 0.05), y=(-0.05, 0.05), z=(0.0, 0.1))
    return transform_coordinate(source_cube, target_cube, vec)
end

function example2magman(vec::AbstractArray{Float64, 1})
    source_cube = (x=(0.00, 0.1), y=(0.0, 0.1), z=(0.0, 1.0))
    target_cube = (x=(-0.05, 0.05), y=(-0.05, 0.05), z=(0.0, 0.1))
    return transform_coordinate(source_cube, target_cube, vec)
end

function get_criteria_cell_indices(csv_path::String)::Vector{Int}
    df = CSV.read(csv_path, DataFrame)
    try 
        return df.vtkOriginalPointIds
    catch ArgumentError
        return df.vtkOriginalCellIds
    end
end

function sample_cells_scalar(cellsField::AbstractArray{Float64, 2}, cell_values::AbstractArray{Float64, 1}, sample_points::AbstractArray{Float64, 2})
    nn_tree = KDTree(cellsField)
    return interpolate_tree(nn_tree, cell_values, sample_points)
end

# function scalar_hotspots_F_generator(hotspots::AbstractArray{Real, 2}, sizes::AbstractVector{Real}, peaks::AbstractVector{Real})
#     return x -> sum(peaks .* exp.())
# end

function set_temp_field(F::Function, casepath::String, t::Int=0; normalize::Bool=false)
    centers_field_path = joinpath(casepath, string(t))
    field_string = create_force_field_string(F, centers_field_path, "./data/fieldsTemplateT", normalize)
end

function set_field_uniform(casepath::String, time::Int, field::String, value::String)
    input_file = joinpath(casepath, string(time), field)
    pattern = r"internalField\s+uniform\s+(\w+);"
    output_file = "/tmp/tmp_T"

    open(output_file, "w+") do out_file
        for line in eachline(input_file)
            if occursin("internalField", line) && (m = match(pattern, line)) !== nothing
                # Replace the matched value with the new value
                new_line = replace(line, "uniform $(m[1]);" => "uniform $value;")
                write(out_file, new_line * "\n")
            else
                write(out_file, line * "\n")
            end
        end
    end
    mv(output_file, input_file, force=true)
end

function set_control_dict_entry(casepath::String, key::String, value::Any)
    if typeof(value) != String
        value = string(value)
    end
    input_file = joinpath(casepath, "system/controlDict")
    pattern = Regex("$key\\s+(\\w+);")
    output_file = "/tmp/tmp_case"

    found = false

    open(output_file, "w+") do out_file
        for line in eachline(input_file)
            if (m = match(pattern, line)) !== nothing
                # Replace the matched value with the new value
                new_line = "$key $value; // written automatically"
                println(new_line)
                write(out_file, new_line * "\n")
                if found
                    @warn "Found multiple matches for $key setting"
                end
                found = true
            else
                write(out_file, line * "\n")
            end
        end
    end
    if !found
        @warn "Tried to change setting $key but found none"
    end
    mv(output_file, input_file, force=true)
end

function read_control_dict_entry(casepath::String, key::String)::String
    input_file = joinpath(casepath, "system/controlDict")
    pattern = Regex("$key\\s+(\\w+);")
    found = false
    val = ""
    for line in eachline(input_file)
        if (m = match(pattern, line)) !== nothing
            # Replace the matched value with the new value
            val = m[1]
            if found
                @warn "Found multiple matches for $key setting"
            end
            found = true
        end
    end
    
    if !found
        @error "Tried to read setting $key but found none"
    end
    return val
end
