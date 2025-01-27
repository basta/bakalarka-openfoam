using NearestNeighbors, JLD, LinearAlgebra

mutable struct ElMagneticInput
    poses::Array{Float64, 3}
    vecs::Array{Float64, 3}
    tree::KDTree
    isEl::Bool
end

mutable struct ElMagneticData
    inputs::Vector{MagneticInput}
end

function interpolate_tree(nn_tree, vecs, points)
	idxssE, distssE = knn(nn_tree, points[:, :], 5)
	vecsE = []
	for (i, idxs) in enumerate(idxssE)
			dists = distssE[i]
			pos = out_grid[:, i]
			i_vecs = vecs[idxs]
			vec2 = Vec3f(0)
			dist_sum = sum(dists)
			w_sum = 0
			for (dist, i_vec) in zip(dists, i_vecs)
				w_i = 1/dist
				w_sum += w_i
				# i_vec = cross(collect(i_vec), [1;0;0])
				vec2 += Vec3f(i_vec[1], i_vec[2], i_vec[3])*w_i
			end
			
			vec2 = vec2 / w_sum
			
			push!(vecsE,vec2)
		end
	return stack(vecsE)
end

function create_elmaginput(jld_path::String)
    data = load(jld_path)
    poses = data["coords"]
    isEl = haskey(data, "Ez")
    if isEl
        vecs = stack(zip(data["Ex"], data["Ey"], data["Ez"]))
    else
        vecs = stack(zip(data["Bx"], data["By"], data["Bz"]))
    end
    return MagneticInput(poses, vecs, KDTree(poses), isEl)
end

function create_trees(jld_paths::AbstractVector{String})
    maginputs = []
    for path in jld_paths
        push!(maginputs, create_maginput(path))
    end
    return MagneticData(maginputs)
end

function get_force_at_point(elMagData::ElMagneticData, point::Array{Float64, 1}, θ::Vector{Float64})
    F = Vec3f(0)
    for (Bi, maginput) in enumerate(elMagData.inputs)
        if maginput.isEl
            continue
        end
        for (Ei, elinput) in elMagData.inputs
            if !elinput.isEl
                continue
            end
            E = interpolate_tree(elinput.tree, elinput.vecs, [point])[1] * θ[Ei]
            B = interpolate_tree(maginput.tree, maginput.vecs, [point])[1] * θ[Bi]
            F += cross(E, B)
        end
    end
    return F
end