using LinearAlgebra
using Statistics


"""
    pod_pca(data_matrix::Matrix{Float64}; num_modes::Union{Nothing, Int}=nothing, center::Bool=true)

Performs Proper Orthogonal Decomposition (POD) or Principal Component Analysis (PCA)
on a data matrix using Singular Value Decomposition (SVD).

Assumes the input `data_matrix` has features/variables as rows and
snapshots/observations as columns (n_features x n_snapshots).

# Arguments
- `data_matrix::Matrix{Float64}`: The input data (features x snapshots).
- `num_modes::Union{Nothing, Int}`: (Optional) Number of dominant modes/components to retain.
                                     If `nothing`, all modes are returned. Defaults to `nothing`.
- `center::Bool`: (Optional) Whether to subtract the mean across snapshots (temporal mean)
                  before performing SVD. Defaults to `true`.

# Returns
- `modes::Matrix{Float64}`: The spatial modes or principal component vectors (columns of U from SVD, possibly truncated).
- `singular_values::Vector{Float64}`: The singular values (possibly truncated). Square relates to energy/variance.
- `temporal_coeffs::Matrix{Float64}`: The temporal coefficients or principal component scores
                                      (projections of data onto modes, calculated as `Diagonal(S) * V'`, possibly truncated).
- `explained_variance_ratio::Vector{Float64}`: The fraction of total variance/energy captured by each mode (calculated from singular values).
- `mean_vector::Union{Nothing, Vector{Float64}}`: The mean vector subtracted if `center=true`, otherwise `nothing`.
"""
function pod_pca(data_matrix::Matrix{<:Real}; num_modes::Union{Nothing, Int}=nothing, center::Bool=true)
    n_features, n_snapshots = size(data_matrix)

    # --- 1. Center Data (Optional) ---
    X = copy(data_matrix) # Work on a copy
    mean_vector = nothing
    if center
        # Calculate mean across snapshots (temporal mean)
        mean_vector = vec(mean(X, dims=2))
        # Subtract mean from each snapshot column
        X .-= mean_vector
        @info "Data centered by subtracting the temporal mean vector."
    else
        @info "Data not centered."
    end

    # --- 2. Perform SVD ---
    # U: Left singular vectors (spatial modes) - n_features x rank
    # S: Singular values (vector) - rank
    # V: Right singular vectors (temporal structure) - n_snapshots x rank
    # rank = min(n_features, n_snapshots)
    @info "Performing SVD..."
    U, S, V = svd(X)
    @info "SVD complete."

    # --- 3. Calculate Explained Variance ---
    total_variance = sum(S.^2)
    explained_variance_ratio = (S.^2) ./ total_variance

    # --- 4. Truncate (Optional) ---
    rank = length(S)
    modes_to_keep = rank # Default to keeping all

    if !isnothing(num_modes)
        if num_modes > 0 && num_modes <= rank
            modes_to_keep = num_modes
            @info "Truncating results to $modes_to_keep modes."
        elseif num_modes > rank
             @warn "Requested num_modes ($num_modes) is greater than the rank ($rank). Returning all $rank modes."
        else
             @warn "Requested num_modes ($num_modes) is non-positive. Returning all $rank modes."
        end
        # Truncate results
        U = U[:, 1:modes_to_keep]
        S = S[1:modes_to_keep]
        V = V[:, 1:modes_to_keep]
        explained_variance_ratio = explained_variance_ratio[1:modes_to_keep]
    else
        @info "Returning all $rank modes."
    end

    # --- 5. Calculate Temporal Coefficients ---
    # These represent the projection of the (centered) data onto the modes.
    # Equivalent to U' * X or Diagonal(S) * V'
    # Size: modes_to_keep x n_snapshots
    temporal_coeffs = Diagonal(S) * V'

    # --- 6. Prepare Outputs ---
    modes = U # Spatial modes are columns of U

    return modes, S, temporal_coeffs, explained_variance_ratio, mean_vector
end

