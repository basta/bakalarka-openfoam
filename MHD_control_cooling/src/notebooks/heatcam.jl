### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ cdcd3247-c1c8-4227-b104-20d4885edb17
using Revise

# ╔═╡ 330930ad-a341-4283-be87-e21a2ee45996
using CairoMakie

# ╔═╡ 3e638860-cf88-4e34-98b0-6bc05152cde9
using JLD

# ╔═╡ da91dd76-1138-4245-ac90-68fd4aeefc12
PKG_PATH = "/home/basta/Projects/bakalarka-openfoam/MHD_control_cooling"

# ╔═╡ 94bee77a-44ec-4b99-bed0-9261253d40ca
begin
    import Pkg
    Pkg.activate(temp=true)
    Pkg.develop(path=PKG_PATH)
	Pkg.add("WGLMakie")
	Pkg.add("CairoMakie")
	Pkg.add("JLD")
	
    using MHD_control_cooling 
end

# ╔═╡ 1141d878-98c0-484e-9b70-d0cf1a0243b2
case = "2d-example"

# ╔═╡ f6486da4-f4ef-4fb6-bbbd-96e13d14d14f
cells = MHD_control_cooling.read_field_vector(
	joinpath(PKG_PATH, "../$(case)/0/C")
)

# ╔═╡ 9bca1197-597a-4dee-b435-7bd25f0566e5
temps5 = MHD_control_cooling.read_field_scalar(
	joinpath(PKG_PATH, "../$(case)/5/T")
)

# ╔═╡ dee612f0-8aa0-4bc4-beb7-5f94d9eed4a2
temps150 = MHD_control_cooling.read_field_scalar(
	joinpath(PKG_PATH, "../$(case)/20/T")
)

# ╔═╡ 688d2f09-0ae3-4a4e-a572-d0f25dd89033
function poses_samples_to_mat(sample_poses, sample_vals, xmin, xstep, xmax, ymin, ystep, ymax)
	xsize = Int(ceil((xmax - xmin)/xstep))
	ysize = Int(ceil((ymax - ymin)/ystep))
	out_mat = zeros(xsize, ysize)
	print(size(out_mat))
	for (i, sample_pos) in enumerate(eachcol(sample_poses))
		out_idx = 
			(sample_pos - [xmin; ymin; 1]) ./
			([xmax - xmin; ymax - ymin; 1])
		out_idx .*= [xsize; ysize; 1]
		out_idx[1] = clamp(out_idx[1], 0, xsize-1)
		out_idx[2] = clamp(out_idx[2], 0, ysize-1)		
		out_mat[Int(round(out_idx[1])+1), Int(round(out_idx[2]))+1] = sample_vals[i]
	end
	return out_mat
end

# ╔═╡ 19ab629e-be36-445a-8751-e314d0f86cf4
function matrix_gradient(A::Matrix{T}) where T<:Number
    rows, cols = size(A)
    grad_y = zeros(T, rows, cols)
    grad_x = zeros(T, rows, cols)
    
    # Interior points - central difference
    for i in 2:rows-1, j in 2:cols-1
        grad_y[i,j] = (A[i+1,j] - A[i-1,j]) / 2
        grad_x[i,j] = (A[i,j+1] - A[i,j-1]) / 2
    end
    
    # Edges - forward/backward difference
    for j in 1:cols
        grad_y[1,j] = A[2,j] - A[1,j]
        grad_y[rows,j] = A[rows,j] - A[rows-1,j]
    end
    
    for i in 1:rows
        grad_x[i,1] = A[i,2] - A[i,1]
        grad_x[i,cols] = A[i,cols] - A[i,cols-1]
    end
    
    return grad_y, grad_x
end

# ╔═╡ 0298b1e2-9120-4592-a1a8-61743b8da9b8
begin
	samples_T = []
	z = 0.005
	for x in 0:(0.1/20):0.1
		for y in 0:(0.1/20):0.1
			push!(samples_T, [x;y;z])
		end
	end
	samples_mat = hcat(samples_T...)
end;

# ╔═╡ 61f19e09-5ee9-479f-b8ff-18bc401be62a
begin
	sampled_temps5 = MHD_control_cooling.sample_cells_scalar(cells, temps5, samples_mat)
	sampled_temps150 = MHD_control_cooling.sample_cells_scalar(cells, temps150, samples_mat)
end;

# ╔═╡ 62c120e3-1569-4369-aae3-a1d12a8d4fbb
heat_mat = poses_samples_to_mat(samples_mat, sampled_temps5, 0, 0.1/20 ,0.1, 0, 0.1/20, 0.1);

# ╔═╡ 37649575-2f96-40a9-8027-281ef3040dc0
let
	CairoMakie.activate!()
	case_noinput = "2d-example"
	case_input = "cases_bckup/2d-example-bb"
	case = case_input
	for t in 1:1:50
		temps = MHD_control_cooling.read_field_scalar(
			joinpath(PKG_PATH, "../$(case)/$(t)/T"))
		sampled_temps = MHD_control_cooling.sample_cells_scalar(cells, temps, 		samples_mat)
		heat_mat = poses_samples_to_mat(samples_mat, sampled_temps, 0, 0.1/20 ,0.1, 0, 0.1/20, 0.1);
		grad_x, grad_y = matrix_gradient(heat_mat)

		begin
			figh = Figure(resolution=(700, 700))	
			ax = Axis(figh[1, 1])
		    heatmap!(ax, heat_mat, cmap=:viridis)
			ax = Axis(figh[1, 2])
		    quiver!(ax, 0:1:32, 0:1:32, grad_x, grad_y, cmap=:viridis)
			ax = Axis(figh[2, 1])
		    heatmap!(ax, grad_x, cmap=:viridis)
			ax = Axis(figh[2, 2])
		    heatmap!(ax, grad_y, cmap=:viridis)
			save("/tmp/grad$(t).png", figh)
		end
	end
	JLD.save("/tmp/T_mat.jld", "T", heat_mat)	
end

# ╔═╡ 7b05973c-99af-46c0-8b5f-5a0dfaeadeda
begin
	x = [col[1] for col in eachcol(samples_mat)]
	y = [col[2] for col in eachcol(samples_mat)]
	fig = Figure(resolution = (400, 800))
	color_range = (20, 80)
	ax1 = Axis(fig[1, 1], title = "Scatter Plot Visualization")
	
	scatter!(ax1, x, y, color = vec(sampled_temps5), colormap = :viridis, markersize = 30, colorrange=color_range)
	
	ax2 = Axis(fig[2, 1], title = "Scatter Plot Visualization")
	scatter!(ax2, x, y, color = vec(sampled_temps150), colormap = :viridis, markersize = 30, colorrange=color_range)
	
	# Colorbar(fig[1, 2], ax, label = "Values")
	fig
end

# ╔═╡ 9fb36d4a-93ba-4d0e-a563-291706a938b4


# ╔═╡ Cell order:
# ╠═da91dd76-1138-4245-ac90-68fd4aeefc12
# ╠═cdcd3247-c1c8-4227-b104-20d4885edb17
# ╠═94bee77a-44ec-4b99-bed0-9261253d40ca
# ╠═330930ad-a341-4283-be87-e21a2ee45996
# ╠═1141d878-98c0-484e-9b70-d0cf1a0243b2
# ╠═f6486da4-f4ef-4fb6-bbbd-96e13d14d14f
# ╠═9bca1197-597a-4dee-b435-7bd25f0566e5
# ╠═dee612f0-8aa0-4bc4-beb7-5f94d9eed4a2
# ╟─688d2f09-0ae3-4a4e-a572-d0f25dd89033
# ╟─19ab629e-be36-445a-8751-e314d0f86cf4
# ╠═0298b1e2-9120-4592-a1a8-61743b8da9b8
# ╠═61f19e09-5ee9-479f-b8ff-18bc401be62a
# ╠═62c120e3-1569-4369-aae3-a1d12a8d4fbb
# ╠═3e638860-cf88-4e34-98b0-6bc05152cde9
# ╠═37649575-2f96-40a9-8027-281ef3040dc0
# ╠═7b05973c-99af-46c0-8b5f-5a0dfaeadeda
# ╠═9fb36d4a-93ba-4d0e-a563-291706a938b4
