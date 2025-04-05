using IterativeSolvers, LowRankApprox

function solve_pinv(X_curr, X_next)
    return pinv(X_curr) * X_next
end

function solve_explicit(X_curr, X_next)
    return (X_curr' \ X_next')'
end

function solve_iterative(X_curr, X_next; solver=lsqr, kwargs...)
    # Define the matrices for the standard Ax = B form derived from
    # Y * X_curr = X_next  =>  X_curr' * Y' = X_next'
    A = X_curr' # size n x m
    B = X_next' # size n x k

    # Solve A * Z = B using the specified iterative solver.
    # IterativeSolvers functions like lsmr/lsqr typically handle matrix B
    # by solving for each column.
    # Z will have size m x k
    println("Starting iterative solver ($(nameof(solver)))...")
    Z = solver(A, B; kwargs...)
    println("Iterative solver finished.")

    # The desired result is Y = Z'
    # Y will have size k x m
    Y = Z'

    return Y
end

function solve_lowrankapprox(X_curr, X_next; pinv_tol=sqrt(eps(real(eltype(X_curr)))), lra_kwargs...)
    # Input dimension check
    if size(X_curr, 2) != size(X_next, 2)
        error("Dimension mismatch: Number of columns in X_curr ($(size(X_curr, 2))) must match number of columns in X_next ($(size(X_next, 2))).")
    end

    m, n = size(X_curr)
    k = size(X_next, 1)

    # Step 1: Compute the principal SVD of X_curr using LowRankApprox
    # psvd returns U, S, V such that X_curr ≈ U * Diagonal(S) * V'
    # U is size m x r, S is vector of length r, V is size n x r (where r is the computed rank)
    println("Computing low-rank SVD of X_curr (size $m x $n)...")
    U, S, V = psvd(X_curr; lra_kwargs...) # Pass rank/tolerance via lra_kwargs
    r = length(S)
    println("Low-rank SVD computed with actual rank r = $r")
    # size(U) = m x r
    # size(S) = r
    # size(V) = n x r

    # Step 2: Compute the reciprocal of singular values with tolerance
    # Create diagonal matrix InvSigma (r x r) where diagonal elements are 1/S[i] if S[i] > pinv_tol, else 0
    inv_S_diag = zeros(eltype(S), r)
    non_zero_count = 0
    for i in 1:r
        if S[i] > pinv_tol
            inv_S_diag[i] = 1 / S[i]
            non_zero_count += 1
        end
        # Otherwise it remains 0.0
    end
    InvSigma = Diagonal(inv_S_diag) # r x r diagonal matrix
    println("Constructed InvSigma: $non_zero_count non-zero singular values above tolerance $pinv_tol")


    # Step 3: Compute Y_approx = X_next * V * InvSigma * U'
    # This approximates Y ≈ X_next * pinv(X_curr)
    # Calculate efficiently: (X_next * V) * (InvSigma * U')
    # Dimensions:
    # X_next: k x n
    # V:      n x r
    # InvSigma: r x r
    # U':     r x m
    # Result Y_approx: k x m

    println("Calculating approximate solution Y (size $k x $m)...")

    # Compute tmp1 = X_next * V (size k x r)
    print("  Calculating X_next * V ... ")
    tmp1 = X_next * V
    println("done (size $(size(tmp1))).")

    # Compute tmp2 = InvSigma * U' (size r x m)
    # U' is already computed efficiently by LAPACK/BLAS from psvd result structure
    print("  Calculating InvSigma * U' ... ")
    tmp2 = InvSigma * U'
    println("done (size $(size(tmp2))).")

    # Compute Y_approx = tmp1 * tmp2 (size k x m)
    print("  Calculating final Y_approx = tmp1 * tmp2 ... ")
    Y_approx = tmp1 * tmp2
    println("done (size $(size(Y_approx))).")

    println("Calculation finished.")

    return Y_approx
end
