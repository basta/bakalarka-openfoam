### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 33877cce-873b-4c97-b87a-242a66abbb66
using DifferenceEquations, WGLMakie, LinearAlgebra, Colors, JLD, PlutoUI, EquivariantOperators, ImageFiltering

# ╔═╡ 5955b343-5af1-43a2-90c0-0ae17ea9a41f


# ╔═╡ b8b450bb-0718-4868-964e-35b56d594c99
begin
	N = 20
	M = 20
	D = 0.1
	# T = rand(N,M)
	T = zeros(N,M)
	# T = [(j/10)^3 for i in 1:N, j in 1:M]
	# T = JLD.load("/tmp/T_mat.jld")["T"]
	T = JLD.load("../../data/nosample_T.jld")["T"]
	

	dx = .1/20
	dy = .1/20
	▽ = Del([dx 0; 0 dy])  # Cell size matrix
	Δ² = Lap((dx, dy); pad = :same, border = :smooth)
	cam_ΔT = 0.0000002
    
end;

# ╔═╡ c394313f-2e65-4287-a481-1b1d3f4216d9
function simulate_thermal_cam(mat, ΔT)
	return round.(mat ./ ΔT).*ΔT
end

# ╔═╡ d05cea9b-a4e0-49ee-80e5-407e19c7d8b5
let
	fig = Figure(resolution=(600, 300))
	ax = Axis(fig[1,1], title="T field")
	heatmap!(ax, simulate_thermal_cam(T, cam_ΔT))
	ax = Axis(fig[1,2], title="T field")
	heatmap!(ax, imfilter(simulate_thermal_cam(T, cam_ΔT), Kernel.gaussian(3)))
	fig
end

# ╔═╡ 89aa2199-75d6-44e4-83f8-3bd59c7a0e88
begin
	T_real = simulate_thermal_cam(T, cam_ΔT)
	imfilter(simulate_thermal_cam(T, cam_ΔT), Kernel.gaussian(1))
end

# ╔═╡ 59981483-90af-4586-b0f1-5a91e511bf45
## Calculate b vector
b, T_lap = let
	b = zeros(N*M*2+2*N+2*M)
	T_lap = real(Δ²(T_real))
	idx = 1 # skip first
	for i in 1:N
		for j in 1:M
			b[idx] = T_lap[i,j]
			idx += 1
		end
	end
	b, T_lap
end

# ╔═╡ 72128444-e1ac-4853-b902-efe8131738c1
ij_to_idx, idx_to_ij = let
function ij_to_idx(i, j, N, M, axis)::Int
	offset = 0
	if axis == :y
		offset = N*M
	end
	return (i-1)*M+j+offset
end
function idx_to_ij(idx, N, M, axis)::Tuple{Int, Int}
	if axis == :y
		idx -= N*M
	end
	i = div(idx - 1, M) + 1
	j = mod(idx - 1, M) + 1
	return i, j
end
	N = 32
	M = 32
	idx = 1
	for i in 1:N
		for j in 1:M
			if ij_to_idx(i,j,N,M, :x) != idx
				println("i:$(i), j:$(j), f:$(ij_to_idx(i,j,N,M, :x)), idx:$idx")
			end
			idx += 1
		end
	end
	for i in 1:N
		for j in 1:M
			if ij_to_idx(i,j,N,M, :y) != idx
				println("i:$(i), j:$(j), f:$(ij_to_idx(i,j,N,M, :y)), idx:$idx")
			end
			idx += 1
		end
	end
	for k in 1024
		i,j = idx_to_ij(k, 32, 32, :x)
		k_out = ij_to_idx(i, j, 32, 32, :x)
		if k_out != k
			println("$k, $i, $j != $k_out")
		end
	end
	ij_to_idx, idx_to_ij
end

# ╔═╡ fa395f8a-0938-4e7c-a19a-1a26a7612d8f
begin
function build_divergence_matrix(N, M)
    # Create matrix for divergence-free condition (∇·u = 0)
    ico_mat = zeros(N*M, 2*N*M)
    
	for i in 1:N
		for j in 1:M
			row = zeros(N*M*2)
			if i == 1
				row[ij_to_idx(i, j, N, M, :x)] = -1
				row[ij_to_idx(i+1, j, N, M, :x)] = 1
			elseif i == N
				row[ij_to_idx(i-1, j, N, M, :x)] = -1
				row[ij_to_idx(i, j, N, M, :x)] = 1
			else
				row[ij_to_idx(i-1, j, N, M, :x)] = -1
				row[ij_to_idx(i+1, j, N, M, :x)] = 1
			end
			
			if j == 1
				row[ij_to_idx(i, j, N, M, :y)] = -1
				row[ij_to_idx(i, j+1, N, M, :y)] = 1
			elseif j == M
				row[ij_to_idx(i, j-1, N, M, :y)] = -1
				row[ij_to_idx(i, j, N, M, :y)] = 1
			else
				row[ij_to_idx(i, j-1, N, M, :y)] = -1
				row[ij_to_idx(i, j+1, N, M, :y)] = 1
			end
       		ico_mat[ij_to_idx(i, j, N,M, :x), :] = row 
		end
	end
    return ico_mat
end

function build_boundary_mat_u(N,M)
	# Builds noslip conditions
	bound_mat = zeros(N*2+M*2, 2*N*M)
	idx = 1
	for i in 1:N
		for j in 1:M
			if i == 1 || i == N
				bound_mat[idx, ij_to_idx(i,j,N,M, :x)] = 1
				idx += 1
			end
			if j == 1 || j == M
				bound_mat[idx, ij_to_idx(i,j,N,M, :y)] = 1
				idx += 1
			end
		end
	end
	return bound_mat
end
end

# ╔═╡ e0e8aa4f-4f05-4afe-86b3-08990d1a1366
## Calculate A matrix
A, ico_mat, Dux, Duy, bound_mat = let
	diag_ux = zeros(N*M)
	diag_uy = zeros(N*M)
	idx = 1
	T_grad = real(▽(T_real))
	for i in 1:N
		for j in 1:M
			diag_ux[idx] = T_grad[i,j][1]
			diag_uy[idx] = T_grad[i,j][2]			
			idx += 1
       end
	end
	ico_mat = build_divergence_matrix(N,M)
	bound_mat = build_boundary_mat_u(N,M)
	
	Dux = Diagonal(diag_ux)
	Duy = Diagonal(diag_uy)
	# ico_mat = build_divergence_matrix(N,M)
	A = [
		Dux Duy;
		ico_mat;
		bound_mat;
	]
	A, ico_mat, Dux, Duy, bound_mat
end;

# ╔═╡ dac1b492-2831-4b53-a9d3-4db02b99ed48


# ╔═╡ 35e468e9-15f2-4860-ac72-66684c3005bf
begin
	fig = Figure()
	ax = Axis(fig[1,1], title="Dux")
	spy!(Dux')
	ax = Axis(fig[1,2], title="Duy")
	spy!(Duy')
	ax = Axis(fig[2,1:2], title="ico_mat")
	spy!(ico_mat')
	ax = Axis(fig[3,1:2], title="ico_mat")
	spy!(bound_mat')
	fig
end

# ╔═╡ 8dd51efd-e648-4935-ad9a-2e8049192412
begin
	λ = 1e-6  # Regularization parameter
	A_reg = [A; λ*I]  # Augment with regularization
	b_reg = [b; zeros(size(A, 2))]  # Augmented RHS
	u_res = A_reg \ b_reg  # Solve regularized system
end

# ╔═╡ 07f71376-a7a4-4e91-913c-9d8d92cb8663
@bind lgth Slider(0.001:0.001:0.01, show_value=true)

# ╔═╡ 17012da8-2804-4e6e-bb5e-abae29b40dff
let
	ux_res = u_res[1:N*M]
	uy_res = u_res[N*M+1:end]
	size = 0.01
	X = []
	Y = []
	U_field = zeros(N, M)
	Ux_field = zeros(N, M)
	Uy_field = zeros(N, M)
	
	for idx in 1:N*M
		i = idx_to_ij(idx, N, M, :x)[1]
		j = idx_to_ij(idx, N, M, :x)[2]
		push!(X, i)
		push!(Y, j)
		U_field[i, j] = sqrt(ux_res[idx]^2 + uy_res[idx]^2)
		U_field[i, j] = ux_res[idx]
		Uy_field[i, j] = uy_res[idx]
		Ux_field[i, j] = ux_res[idx] 
		
	end
	X = Float64.(X)
	Y = Float64.(Y)
	mags = sqrt.(ux_res.^2 .+ uy_res.^2)

	fig = Figure(resolution=(700, 700))
	ax = Axis(fig[1, 1], title = "Quiver Plot")
	arrows!(ax, X, Y, ux_res, uy_res, arrowsize = 10, lengthscale = lgth,
    arrowcolor = mags, linecolor = mags)
	ax = Axis(fig[2, 1], title = "U mag")
	heatmap!(U_field)
	ax = Axis(fig[1:2, 2], title = "Ux;Uy")
	uxuy = heatmap!([Ux_field Uy_field], colormap = :viridis)
	Colorbar(fig[1:2,3], uxuy)
	fig
end

# ╔═╡ f63929b4-6ce7-44c3-a00d-131fb57d12d4
@bind ΔT2 Slider(0:0.01:3, show_value=true)

# ╔═╡ 0b6985c2-cc3f-4888-9352-97b109f4e5f0
let
	Δ² = Lap((dx, dy); pad = :same, border = :smooth)
	T = imfilter(T, Kernel.gaussian(3))
    
    # Compute gradient components
    grad_T = real(▽(T))[:,:]  # x-component
	lap_T = real(Δ²(T))[:,:]

	T_real = simulate_thermal_cam(T, ΔT2)
	grad_T_real = real(▽(T_real))[:,:]
	lap_T_real = real(Δ²(T_real))[:,:]
	grad_x = zeros(N,M)
	grad_y = zeros(N,M)
	grad_x_real = zeros(N,M)
	grad_y_real = zeros(N,M)
	fig = Figure(resolution=(680, 680))
	ax = Axis(fig[1,1])
	for i in 1:N
		for j in 1:M
			grad_x[i,j] = grad_T[i,j][1]
			grad_y[i,j] = grad_T[i,j][2]
			grad_x_real[i,j] = grad_T_real[i,j][1]
			grad_y_real[i,j] = grad_T_real[i,j][2]
			
		end
	end
	heatmap!(ax, grad_x)
	ax = Axis(fig[2,1])
	heatmap!(ax, lap_T)
	ax = Axis(fig[1,2])
	heatmap!(ax, grad_x_real)
	ax = Axis(fig[2,2])
	heatmap!(ax, lap_T_real)
	fig
end

# ╔═╡ ca06947c-d6ff-4c54-9ee0-8e3990a51205
T_lap

# ╔═╡ Cell order:
# ╠═33877cce-873b-4c97-b87a-242a66abbb66
# ╠═5955b343-5af1-43a2-90c0-0ae17ea9a41f
# ╠═b8b450bb-0718-4868-964e-35b56d594c99
# ╠═c394313f-2e65-4287-a481-1b1d3f4216d9
# ╠═d05cea9b-a4e0-49ee-80e5-407e19c7d8b5
# ╠═89aa2199-75d6-44e4-83f8-3bd59c7a0e88
# ╠═59981483-90af-4586-b0f1-5a91e511bf45
# ╠═72128444-e1ac-4853-b902-efe8131738c1
# ╠═fa395f8a-0938-4e7c-a19a-1a26a7612d8f
# ╠═e0e8aa4f-4f05-4afe-86b3-08990d1a1366
# ╠═dac1b492-2831-4b53-a9d3-4db02b99ed48
# ╠═35e468e9-15f2-4860-ac72-66684c3005bf
# ╠═8dd51efd-e648-4935-ad9a-2e8049192412
# ╠═07f71376-a7a4-4e91-913c-9d8d92cb8663
# ╠═17012da8-2804-4e6e-bb5e-abae29b40dff
# ╠═f63929b4-6ce7-44c3-a00d-131fb57d12d4
# ╠═0b6985c2-cc3f-4888-9352-97b109f4e5f0
# ╠═ca06947c-d6ff-4c54-9ee0-8e3990a51205
