using LinearAlgebra

"""
    find_closest_outer_product(theta::AbstractVector{<:Real}, m::Int, n::Int)

Finds the vectors `u` (size `m`) and `v` (size `n`) such that their outer
product (`u * v'`) produces a matrix `Θ` whose flattened version is closest
to the target vector `theta` in the Euclidean norm sense.

This is equivalent to finding the best rank-1 approximation of the matrix `Θ`
formed by reshaping `theta`, using Singular Value Decomposition (SVD).

# Assumptions
- The input vector `theta` represents the elements of the target `m x n`
  matrix `Θ` flattened in **column-major** order (like Julia's `vec` function). # <--- CORRECTED

# Arguments
- `theta`: The target vector (AbstractVector{<:Real}). Accepts any 1D array type.
- `m`: The desired dimension (number of rows) of the matrix `Θ`, and the length of vector `u` (Int).
- `n`: The desired dimension (number of columns) of the matrix `Θ`, and the length of vector `v` (Int).

# Returns
- `Tuple{Vector, Vector}`: A tuple `(u_opt, v_opt)` containing the optimal vectors (as standard Vectors).

# Throws
- `ErrorException`: If `m * n` does not equal the length of `theta` or if m/n are not positive.
- `DimensionMismatch`: Propagated from `reshape` if dimensions are inconsistent despite passing the length check.
"""
function find_closest_outer_product(theta::AbstractVector{<:Real}, m::Int, n::Int)
    L = length(theta)
    # Input validation
    if m * n != L
        error("Dimensions m=$m, n=$n do not match theta length L=$L")
    end
    if m <= 0 || n <= 0
        error("Dimensions m and n must be positive integers.")
    end

    # Reshape theta into matrix Theta.
    # Assumes theta is column-major flattened. Julia's reshape fills column-by-column.
    local Theta
    try
        theta_vec = vec(theta) # Ensure it's treated as a 1D structure for reshape
        # Reshape directly to m x n for column-major input
        Theta = reshape(theta_vec, m, n) # <--- MODIFIED FOR CLARITY (Equivalent to transpose(reshape(..n,m)))
    catch e
        if isa(e, DimensionMismatch)
            error("Reshaping failed. Check if m=$m * n=$n matches length(theta)=$L. Original error: $e")
        else
            rethrow(e)
        end
    end

    # Compute Singular Value Decomposition (SVD)
    F = svd(Theta)
    U, S, V = F.U, F.S, F.V

    # Extract components for the best rank-1 approximation
    sigma1 = S[1]
    u1 = U[:, 1]
    v1 = V[:, 1]

    # Calculate the optimal u and v vectors
    sqrt_sigma1 = sqrt(max(0.0, sigma1))
    u_opt = sqrt_sigma1 * u1
    v_opt = sqrt_sigma1 * v1

    return Vector(u_opt), Vector(v_opt)
end

# # --- Example Usage ---
# println("--- Example 1: Noisy Data ---")
# # Example target vector theta (replace with your actual vector)
# # Should approximate outer product of u=[1, 2] and v=[3, 4, 5]
# theta_test1 = [3.1, 3.9, 5.2, 5.8, 8.1, 9.9]
# m1 = 2 # length of u
# n1 = 3 # length of v

# try
#     u_opt1, v_opt1 = find_closest_outer_product(theta_test1, m1, n1)

#     println("Input theta: ", theta_test1)
#     println("Optimal u: ", round.(u_opt1, digits=4))
#     println("Optimal v: ", round.(v_opt1, digits=4))

#     # Verification (optional) - Use outer product u * v'
#     Theta_reconstructed1 = u_opt1 * v_opt1'
#     theta_reconstructed1 = vec(transpose(Theta_reconstructed1)) # Flatten row-major

#     println("Reconstructed Matrix Θ:\n", round.(Theta_reconstructed1, digits=4))
#     println("Reconstructed theta: ", round.(theta_reconstructed1, digits=4))

#     # Calculate reconstruction error
#     Theta_original = transpose(reshape(theta_test1, n1, m1))
#     error_norm = norm(Theta_original - Theta_reconstructed1)
#     println("Reconstruction Error (Frobenius Norm): ", round(error_norm, digits=4))
# catch e
#     println("Error during Example 1: ", sprint(showerror, e)) # Print error details
# end


# println("\n--- Example 2: Perfect Case ---")
# # Example with a perfect outer product
# u_true = [1.0, 2.0, 3.0]
# v_true = [4.0, 5.0]
# Theta_true = u_true * v_true'
# theta_test2 = vec(transpose(Theta_true)) # Flatten row-major (result: [4.0, 5.0, 8.0, 10.0, 12.0, 15.0])
# m2 = 3
# n2 = 2

# try
#     # Call should now work with AbstractVector signature
#     u_opt2, v_opt2 = find_closest_outer_product(theta_test2, m2, n2)

#     println("Input theta: ", theta_test2)
#     # Note: SVD solution might have sign flips in u and v compared to original,
#     # e.g., (-u_true) * (-v_true)' gives the same matrix.
#     println("Optimal u: ", round.(u_opt2, digits=4))
#     println("Optimal v: ", round.(v_opt2, digits=4))

#     # Verification - Use outer product u * v'
#     Theta_reconstructed2 = u_opt2 * v_opt2'
#     error_norm2 = norm(Theta_true - Theta_reconstructed2)
#     println("Reconstruction Error (Frobenius Norm): ", round(error_norm2, digits=8)) # Should be near zero
# catch e
#      println("Error during Example 2: ", sprint(showerror, e)) # Print error details
# end

# println("\n--- Example 3: Invalid Dimensions ---")
# theta_test3 = [1.0, 2.0, 3.0, 4.0, 5.0] # Length 5 (prime)
# m3 = 2
# n3 = 3 # m3 * n3 = 6 != 5

# try
#     u_opt3, v_opt3 = find_closest_outer_product(theta_test3, m3, n3)
# catch e
#     println("Caught expected error for invalid dimensions: ")
#     println(sprint(showerror, e)) # Print the error message
# end
