using LinearAlgebra
using SparseArrays # Required for Diagonal and I (identity matrix)
using EquivariantOperators # Assuming Del, Lap are defined here

# --- Helper functions (assuming they are defined as before) ---
# function ij_to_idx(i, j, N, M, axis)::Int ... end
# function idx_to_ij(idx, N, M, axis)::Tuple{Int, Int} ... end
# function build_divergence_matrix(N, M, dx, dy) ... end
# function build_boundary_mat_u(N, M) ... end
# --- (Include the definitions from the previous response if needed) ---

"""
    ij_to_idx(i, j, N, M, axis)::Int

Converts 2D grid indices (i, j) to a 1D index for a flattened vector.
Handles separate indexing for x and y components. Column-major order.
"""
function ij_to_idx(i, j, N, M, axis)::Int
    offset = 0
    if axis == :y
        offset = N*M # y-components start after all x-components
    end
    # Column-major flattening within each component block
    return (j-1)*N + i + offset
end

"""
    idx_to_ij(idx, N, M, axis)::Tuple{Int, Int}

Converts a 1D index back to 2D grid indices (i, j).
Handles separate indexing for x and y components. Column-major order.
"""
function idx_to_ij(idx, N, M, axis)::Tuple{Int, Int}
    if axis == :y
        idx -= N*M # Adjust index if it's for the y-component
    end
    # Inverse of column-major flattening
    j = div(idx - 1, N) + 1
    i = mod(idx - 1, N) + 1
    return i, j
end


"""
    build_divergence_matrix(N, M, dx, dy)

Builds the sparse matrix representing the discretized divergence operator (∇·u = 0 constraint).
Uses central differences (adjusting at boundaries). Assumes dx and dy are scalar cell dimensions.
"""
function build_divergence_matrix(N, M, dx, dy)
    num_cells = N * M
    num_vars = 2 * num_cells
    I_vec = Int[] # Renamed to avoid conflict with LinearAlgebra.I
    J_vec = Int[] # Renamed to avoid conflict
    V_vec = Float64[] # Renamed to avoid conflict
    cell_to_row(i, j) = (j-1)*N + i # Map cell (i,j) to a unique row index k

    for j in 1:M # Iterate column-wise
        for i in 1:N
            k = cell_to_row(i, j) # Row index for this cell's divergence equation

            # du/dx part
            idx_ux_curr = ij_to_idx(i, j, N, M, :x)
            if i == 1 # Forward difference at left boundary (i=1)
                 idx_ux_right = ij_to_idx(i+1, j, N, M, :x)
                 push!(I_vec, k); push!(J_vec, idx_ux_curr);  push!(V_vec, -1.0 / dx)
                 push!(I_vec, k); push!(J_vec, idx_ux_right); push!(V_vec,  1.0 / dx)
            elseif i == N # Backward difference at right boundary (i=N)
                 idx_ux_left = ij_to_idx(i-1, j, N, M, :x)
                 push!(I_vec, k); push!(J_vec, idx_ux_left); push!(V_vec, -1.0 / dx)
                 push!(I_vec, k); push!(J_vec, idx_ux_curr); push!(V_vec,  1.0 / dx)
            else # Central difference in interior
                 idx_ux_left = ij_to_idx(i-1, j, N, M, :x)
                 idx_ux_right = ij_to_idx(i+1, j, N, M, :x)
                 push!(I_vec, k); push!(J_vec, idx_ux_left);  push!(V_vec, -0.5 / dx)
                 push!(I_vec, k); push!(J_vec, idx_ux_right); push!(V_vec,  0.5 / dx)
            end

            # dv/dy part
            idx_uy_curr = ij_to_idx(i, j, N, M, :y)
            if j == 1 # Forward difference at bottom boundary (j=1)
                idx_uy_up = ij_to_idx(i, j+1, N, M, :y)
                push!(I_vec, k); push!(J_vec, idx_uy_curr); push!(V_vec, -1.0 / dy)
                push!(I_vec, k); push!(J_vec, idx_uy_up);   push!(V_vec,  1.0 / dy)
            elseif j == M # Backward difference at top boundary (j=M)
                 idx_uy_down = ij_to_idx(i, j-1, N, M, :y)
                 push!(I_vec, k); push!(J_vec, idx_uy_down); push!(V_vec, -1.0 / dy)
                 push!(I_vec, k); push!(J_vec, idx_uy_curr); push!(V_vec,  1.0 / dy)
            else # Central difference in interior
                 idx_uy_down = ij_to_idx(i, j-1, N, M, :y)
                 idx_uy_up   = ij_to_idx(i, j+1, N, M, :y)
                 push!(I_vec, k); push!(J_vec, idx_uy_down); push!(V_vec, -0.5 / dy)
                 push!(I_vec, k); push!(J_vec, idx_uy_up);   push!(V_vec,  0.5 / dy)
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec, num_cells, num_vars)
end

"""
    build_boundary_mat_u(N, M)

Builds the sparse matrix representing no-slip boundary conditions (u=0 or v=0 on walls).
Assumes u=0 on i=1, i=N boundaries and v=0 on j=1, j=M boundaries.
"""
function build_boundary_mat_u(N, M)
    num_vars = 2 * N * M
    num_bc_eqs = 2 * N + 2 * M # u=0 on top/bottom (2M), v=0 on left/right (2N)
    I_vec = Int[] # Renamed
    J_vec = Int[] # Renamed
    V_vec = Float64[] # Renamed
    eq_idx = 1 # Equation index

    # Enforce v = 0 on left (j=1) and right (j=M) walls
    for i in 1:N
        # Left wall (j=1)
        idx_v_left = ij_to_idx(i, 1, N, M, :y)
        push!(I_vec, eq_idx); push!(J_vec, idx_v_left); push!(V_vec, 1.0)
        eq_idx += 1
        # Right wall (j=M)
        idx_v_right = ij_to_idx(i, M, N, M, :y)
        push!(I_vec, eq_idx); push!(J_vec, idx_v_right); push!(V_vec, 1.0)
        eq_idx += 1
    end

     # Enforce u = 0 on bottom (i=1) and top (i=N) walls
    for j in 1:M
         # Bottom wall (i=1)
        idx_u_bottom = ij_to_idx(1, j, N, M, :x)
        push!(I_vec, eq_idx); push!(J_vec, idx_u_bottom); push!(V_vec, 1.0)
        eq_idx += 1
        # Top wall (i=N)
        idx_u_top = ij_to_idx(N, j, N, M, :x)
        push!(I_vec, eq_idx); push!(J_vec, idx_u_top); push!(V_vec, 1.0)
        eq_idx += 1
    end

    @assert eq_idx - 1 == num_bc_eqs "Mismatch in number of boundary condition equations"
    return sparse(I_vec, J_vec, V_vec, num_bc_eqs, num_vars)
end


"""
    calculate_velocity_from_temperature(
        T_field::AbstractMatrix{<:Real},
        dTdt_field::AbstractMatrix{<:Real},
        dx::Real,
        dy::Real;
        α::Real=1.0,
        λ::Real=1e-6
    )

Calculates the 2D velocity field (Ux, Uy) from a given 2D temperature field (T_field)
and its time derivative (dTdt_field), based on the unsteady heat equation:
∂T/∂t + u ⋅ ∇T = α ∇²T
This is rearranged to solve for u:
u ⋅ ∇T = α ∇²T - ∂T/∂t
subject to ∇ ⋅ u = 0 and no-slip boundary conditions.

Args:
    T_field: An N x M matrix representing the temperature field at a specific time.
    dTdt_field: An N x M matrix representing the partial derivative of temperature
                with respect to time (∂T/∂t) at that same time.
    dx: Grid spacing in the x-direction (first dimension).
    dy: Grid spacing in the y-direction (second dimension).
    α: Thermal diffusivity. Defaults to 1.0.
    λ: Regularization parameter for the linear system solver. Defaults to 1e-6.

Returns:
    A tuple (Ux_field, Uy_field) containing two N x M matrices for the x and y velocity components.

Raises:
    DimensionMismatch: If T_field and dTdt_field do not have the same dimensions.
"""
function calculate_velocity_from_temperature(
    T_field::AbstractMatrix{<:Real},
    dTdt_field::AbstractMatrix{<:Real},
    dx::Real,
    dy::Real;
    α::Real=1.0,
    λ::Real=1e-6
    )

    if size(T_field) != size(dTdt_field)
        throw(DimensionMismatch("T_field and dTdt_field must have the same dimensions"))
    end

    N, M = size(T_field)
    num_cells = N * M
    num_vars = 2 * num_cells

    # Ensure dx and dy are Float64 for operators
    dx_f::Float64 = Float64(dx)
    dy_f::Float64 = Float64(dy)

    # Define differential operators
    ▽ = Del((dx_f, dy_f)) # Gradient operator
    Δ² = Lap((dx, dy); pad = :same, border = :smooth)
    
    # 1. Calculate Laplacian and Gradient of the Temperature Field
    T_lap = real(Δ²(T_field))
    T_grad = real(▽(T_field)) # Array of 2-element vectors [∂T/∂x, ∂T/∂y]

    # 2. Build the right-hand side vector 'b' (representing α*∇²T - ∂T/∂t)
    # Flatten (α * T_lap - dTdt_field) into 'b' (column-major)
    b_source_term = α .* T_lap .- dTdt_field
    b = vec(b_source_term) # Equivalent to reshape(b_source_term, num_cells)

    # 3. Build the system matrix 'A'
    # 3a. Construct diagonal matrices from Temperature Gradient (∇T)
    diag_ux = zeros(num_cells)
    diag_uy = zeros(num_cells)
    idx = 1
    for j in 1:M # Column-major iteration
        for i in 1:N
            diag_ux[idx] = T_grad[i, j][1] # ∂T/∂x coefficient for ux
            diag_uy[idx] = T_grad[i, j][2] # ∂T/∂y coefficient for uy
            idx += 1
        end
    end
    Dux = Diagonal(diag_ux)
    Duy = Diagonal(diag_uy)

    # 3b. Build constraint matrices
    ico_mat = build_divergence_matrix(N, M, dx_f, dy_f) # ∇ ⋅ u = 0
    bound_mat = build_boundary_mat_u(N, M)             # u=0 or v=0 on boundaries

    # 3c. Assemble the 'A' matrix components
    A_eq = [Dux Duy] # Represents u ⋅ ∇T part (N*M equations)
    A = vcat(A_eq, ico_mat, bound_mat) # Stack equation, divergence, boundary constraints

    # 4. Regularize and Solve the Linear System
    b_augmented = [b; zeros(size(ico_mat, 1) + size(bound_mat, 1))] # Pad b for constraint eqns

    # --- Corrected Regularization Term ---
    # Create a sparse identity matrix of size num_vars x num_vars and scale by λ
    regularization_term = λ * sparse(I, num_vars, num_vars)
    A_reg = vcat(A, regularization_term) # Add regularization term rows
    # ------------------------------------

    b_reg = [b_augmented; zeros(num_vars)] # Pad b for regularization rows

    # Solve the regularized linear system A_reg * u = b_reg
    println("Solving system with size(A_reg) = ", size(A_reg), " and size(b_reg) = ", size(b_reg))
    u_res = A_reg \ b_reg

    # 5. Reshape the result into 2D velocity fields
    ux_flat = u_res[1:num_cells]
    uy_flat = u_res[num_cells+1 : end]

    Ux_field = reshape(ux_flat, N, M)
    Uy_field = reshape(uy_flat, N, M)

    return Ux_field, Uy_field
end

# == Example Usage ==
# N = 20
# M = 20
# dx = 0.1
# dy = 0.1
# α_physical = 1.5e-5 # Example thermal diffusivity for water
# λ_reg = 1e-6

# # Assume T_field is your temperature data [N x M] at time t
# T_field = reshape(X_data[:, 1], N, M) # Example

# # You need dTdt_field [N x M] at the same time t
# # This could come from finite differences (e.g., (T(t) - T(t-Δt))/Δt),
# # another model, or experimental data.
# dTdt_field = zeros(N, M) # Placeholder: If T is steady, dT/dt is zero

# # Call the updated function
# Ux, Uy = calculate_velocity_from_temperature(T_field, dTdt_field, dx, dy, α=α_physical, λ=λ_reg)

# println("Calculated Ux field size: ", size(Ux))
# println("Calculated Uy field size: ", size(Uy))

# # You can then use the Makie plotting code from the previous step with these Ux, Uy
