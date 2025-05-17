# ############# evaluate_interpolation.jl #############
# Script to evaluate the Inverse Distance Weighting (IDW) interpolation
# used in the force field generation for thesis.
# It calculates the force field on the target OpenFOAM mesh for all
# 16 combinations of 4D basis vectors for electric (u) and magnetic (w) weights,
# accounting for a coordinate system shift between the source JLD data and the target mesh.
#
# Generates:
# 1. A CSV file ('theta_combinations_force_data.csv') with coordinates and all 16
#    calculated force fields (components + magnitude) for visualization.
# 2. A CSV file ('theta_combinations_metrics.csv') with quantitative metrics
#    (MAE, RMSE, Max Error) comparing each combination to a reference combination.
# #####################################################

using NearestNeighbors
using LinearAlgebra
# using WriteVTK # Removed dependency
# using CairoMakie # Not used in this version
using CSV
using DataFrames
using Printf
using Statistics # For mean
using JLD2 # For loading elMagData

# --- Configuration ---

# !! ADJUST THESE PATHS !!
# Base directory for output files (CSV)
const OUTPUT_DIR = "./figures"
# Directory containing the source .jld files (E*.jld, B*.jld)
const JLD_DATA_DIR = "./data/forceFields" # Adjust to your JLD directory
# Directory containing the target OpenFOAM case
const TARGET_CASE_DIR = "../cyl"
# Relative path to the cell centers file within TARGET_CASE_DIR
const TARGET_CENTERS_PATH = "0/C"

# Coordinate Shift: Vector to subtract from OpenFOAM coordinates
# to get the equivalent coordinates in the JLD data's system.
# Example: If OpenFOAM mesh is [0, 0.1] and JLD is [-0.05, 0.05],
# the shift is [0.05, 0.05, 0.0].
const COORD_SHIFT_VECTOR = [0.05, 0.05, 0.0] # Shift applied as: foam_coord - SHIFT = jld_coord

# IDW Parameters (Used for force calculation)
const INTERPOLATION_K = 5       # Default number of neighbors for interpolation
const INTERPOLATION_P = 2.0     # Default power for IDW (assuming 1/dist^2)

# Define 4D Basis Vectors
const BASIS_VECTORS = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0]
]
# Helper to create string representation for column names
basis_to_string(v) = join(Int.(v), "")

# Define the reference combination for metrics calculation
const REF_U_VEC = BASIS_VECTORS[1] # e.g., [1,0,0,0]
const REF_W_VEC = BASIS_VECTORS[1] # e.g., [1,0,0,0]
const REF_THETA = vcat(REF_U_VEC, REF_W_VEC)
const REF_U_STR = basis_to_string(REF_U_VEC)
const REF_W_STR = basis_to_string(REF_W_VEC)
const REF_COL_PREFIX = "F_u$(REF_U_STR)_w$(REF_W_STR)"

# --- Include User Code ---
# Make sure these files are accessible from where you run this script
try
    include("../openfoam/real_field.jl")
    include("../openfoam/field_utils.jl") # Assuming read_foam_vector_field is here
# If MHD_control_cooling module wraps these, adjust accordingly
# using .MHD_control_cooling
catch e
    @error "Failed to include necessary user scripts (real_field.jl, field_utils.jl). Ensure they are in the correct path. Error: $e"
    exit()
end


# --- Utility Functions ---

# Basic function to find the start of the data list in OpenFOAM files
function find_foam_data_start(io::IO)
    line_num = 0
    while !eof(io)
        line = readline(io)
        line_num += 1
        stripped_line = strip(line)
        if length(stripped_line) > 0 && (startswith(stripped_line, "(") || isdigit(stripped_line[1]))
            if isdigit(stripped_line[1])
                n_items = parse(Int, stripped_line)
                line = readline(io) # Read the '('
                line_num += 1
                if strip(line) == "("
                    return line_num, n_items
                else
                    error("Expected '(' after item count, found '$line'")
                end
            elseif startswith(stripped_line, "(")
                # Size might be omitted, return -1
                return line_num, -1
            end
        end
        # Limit search depth to avoid infinite loops on malformed files
        if line_num > 50
            error("Could not find data start marker '(' within 50 lines.")
        end
    end
    error("Reached end of file without finding data start marker '('")
end

# Read OpenFOAM vector field (like C or F) -> Matrix{Float64} (3 x N)
function read_foam_vector_field(field_path::String)::Matrix{Float64}
    if !isfile(field_path)
        error("Field file not found: $field_path")
    end
    vectors = Vector{Float64}[]
    open(field_path, "r") do io
        try
            start_line, n_items = find_foam_data_start(io)

            while !eof(io)
                line = strip(readline(io))
                if line == ")"
                    break
                end
                if !isempty(line) && !startswith(line, "//")
                    # Remove parentheses and split
                    parts = split(replace(line, r"[()]" => ""), keepempty=false)
                    if length(parts) == 3
                        try
                            push!(vectors, [parse(Float64, p) for p in parts])
                        catch e
                            @warn "Skipping invalid vector line: '$line' ($e)"
                        end
                    else
                        @warn "Skipping line with unexpected number of parts: '$line'"
                    end
                end
            end

            if n_items != -1 && length(vectors) != n_items
                @warn "Expected $n_items vectors in $field_path, found $(length(vectors))."
            end

        catch e
            @error "Error reading field file '$field_path': $e"
            rethrow(e)
        end
    end
    if isempty(vectors)
        @warn "No vectors read from $field_path. Returning empty matrix."
        return Matrix{Float64}(undef, 3, 0)
    end
    # Convert Vector{Vector{Float64}} to Matrix{Float64} (3 x N)
    return stack(vectors)
end



"""
    interpolate_tree_k(nn_tree, vecs, points, k, p)

Performs IDW interpolation using k neighbors and power p.
Adapted from the user's interpolate_tree to accept k and p.

Args:
    nn_tree (KDTree): KDTree built on source points (in JLD coordinate system).
    vecs (Matrix{Float64}): Source field data (3 x N_source).
    points (Matrix{Float64}): Target points (3 x N_target) *already shifted* into the JLD coordinate system.
    k (Int): Number of neighbors.
    p (Float64): IDW power.

Returns:
    Matrix{Float64}: Interpolated field data at target points (3 x N_target).
"""
function interpolate_tree_k(nn_tree, vecs::AbstractMatrix{Float64}, points::AbstractMatrix{Float64}, k::Int, p::Float64)
    n_target = size(points, 2)
    n_dim_field = size(vecs, 1)
    interpolated_vecs = zeros(Float64, n_dim_field, n_target)
    epsilon = 1e-12

    warning_count = 0
    warning_threshold = 5

    for i in 1:n_target
        target_point = points[:, i] # This point is already shifted
        idxs, dists = knn(nn_tree, target_point, k, true) # Get sorted distances

        sum_weights = 0.0
        weighted_sum = zeros(Float64, n_dim_field)
        valid_neighbors = 0

        # Handle exact match or very close point first
        if !isempty(dists) && dists[1] < epsilon
            weighted_sum = vecs[:, idxs[1]]
            sum_weights = 1.0
            valid_neighbors = 1
        else
            # Calculate weights for k neighbors
            for j in 1:length(idxs)
                idx = idxs[j]
                dist = dists[j]
                if dist < epsilon
                    continue
                end # Safety check

                weight = 1.0 / (dist^p)
                weighted_sum += weight * vecs[:, idx]
                sum_weights += weight
                valid_neighbors += 1
            end
        end

        # Assign value
        if valid_neighbors > 0 && sum_weights > epsilon
            interpolated_vecs[:, i] = weighted_sum / sum_weights
        else
            # Fallback: Use nearest neighbor
            idxs_nn, dists_nn = knn(nn_tree, target_point, 1)
            if !isempty(idxs_nn)
                interpolated_vecs[:, i] = vecs[:, idxs_nn[1]]
                if warning_count < warning_threshold
                    @warn "Weighting failed for target point $i (dist ~ $(dists_nn[1])), using nearest neighbor value (k=$k)."
                elseif warning_count == warning_threshold
                    @warn "Further weighting failure warnings suppressed..."
                end
                warning_count += 1
            else
                if warning_count < warning_threshold
                    @warn "Could not find any neighbors for target point $i. Setting to zero."
                elseif warning_count == warning_threshold
                    @warn "Further neighbor finding failure warnings suppressed..."
                end
                interpolated_vecs[:, i] .= 0.0 # Assign zero vector
                warning_count += 1
            end
        end
    end
    if warning_count > 0
        @warn "$warning_count points encountered issues during interpolation with k=$k."
    end
    return interpolated_vecs
end


"""
    get_force_at_points_k(elMagData, points_foam, θ, k, p, shift_vector)

Calculates the force field at multiple OpenFOAM points using specified k and p for IDW,
applying a coordinate shift before interpolation.

Args:
    elMagData (ElMagneticData): Loaded data from JLD files.
    points_foam (Matrix{Float64}): Target points in OpenFOAM coordinates (3 x N_target).
    θ (AbstractVector): Input parameters for combining E/B fields (length 8).
    k (Int): Number of neighbors for IDW.
    p (Float64): Power for IDW.
    shift_vector (Vector{Float64}): Vector to subtract from points_foam to get JLD coordinates.

Returns:
    Matrix{Float64}: Calculated force vectors at target points (3 x N_target).
"""
function get_force_at_points_k(elMagData::ElMagneticData, points_foam::Matrix{Float64}, θ::AbstractVector{<:Number}, k::Int, p::Float64, shift_vector::Vector{Float64})
    n_points = size(points_foam, 2)
    F_total = zeros(Float64, 3, n_points)

    if length(θ) != 8
        @error "Theta vector must have length 8 (4 for E, 4 for B). Got length $(length(θ))."
        return F_total # Return zeros
    end
    if length(shift_vector) != 3
         @error "Shift vector must have length 3. Got length $(length(shift_vector))."
        return F_total
    end

    # Shift target points to the JLD coordinate system *before* interpolation
    points_jld = points_foam .- shift_vector # Broadcasting subtraction

    # Pre-interpolate all E and B fields at all *shifted* target points for efficiency
    interpolated_fields = Dict{Int,Matrix{Float64}}() # Store interpolated E/B vectors for each input source
    num_el_sources = 0
    num_mag_sources = 0

    for (idx, input_data) in enumerate(elMagData.inputs)
        is_el = input_data.isEl
        local theta_idx::Int
        if is_el
            num_el_sources += 1
            if num_el_sources > 4
                @warn "More than 4 electric field sources found in elMagData. Check data structure and theta mapping."
                continue
            end
            theta_idx = num_el_sources
        else
            num_mag_sources += 1
            if num_mag_sources > 4
                @warn "More than 4 magnetic field sources found in elMagData. Check data structure and theta mapping."
                continue
            end
            theta_idx = num_mag_sources + 4
        end

        if abs(θ[theta_idx]) > 1e-9
            # Interpolate using the shifted points (points_jld)
            interpolated_fields[idx] = interpolate_tree_k(input_data.tree, input_data.vecs, points_jld, k, p)
        else
            interpolated_fields[idx] = zeros(Float64, 3, n_points)
        end
    end

    # Calculate forces by combining pre-interpolated fields
    el_idx_map = Dict{Int,Int}() # Map original index to theta index (1-4)
    mag_idx_map = Dict{Int,Int}() # Map original index to theta index (5-8)
    el_count = 0
    mag_count = 0
    for (idx, input_data) in enumerate(elMagData.inputs)
        if input_data.isEl
            el_count += 1
            if el_count <= 4
                el_idx_map[idx] = el_count
            end
        else
            mag_count += 1
            if mag_count <= 4
                mag_idx_map[idx] = mag_count + 4
            end
        end
    end


    for (Bi_orig, maginput) in enumerate(elMagData.inputs)
        if maginput.isEl continue end
        if !haskey(mag_idx_map, Bi_orig) continue end
        Bi_theta = mag_idx_map[Bi_orig]
        if abs(θ[Bi_theta]) < 1e-9 continue end

        for (Ei_orig, elinput) in enumerate(elMagData.inputs)
            if !elinput.isEl continue end
            if !haskey(el_idx_map, Ei_orig) continue end
            Ei_theta = el_idx_map[Ei_orig]
            if abs(θ[Ei_theta]) < 1e-9 continue end

            # Retrieve pre-interpolated fields using original index
            # These fields were calculated based on the shifted points (points_jld)
            # but the resulting force applies at the original OpenFOAM points.
            E_vectors = interpolated_fields[Ei_orig] .* θ[Ei_theta]
            B_vectors = interpolated_fields[Bi_orig] .* θ[Bi_theta]

            # Calculate cross product for all points element-wise
            F_total[1, :] .+= E_vectors[2, :] .* B_vectors[3, :] .- E_vectors[3, :] .* B_vectors[2, :]
            F_total[2, :] .+= E_vectors[3, :] .* B_vectors[1, :] .- E_vectors[1, :] .* B_vectors[3, :]
            F_total[3, :] .+= E_vectors[1, :] .* B_vectors[2, :] .- E_vectors[2, :] .* B_vectors[1, :]
        end
    end
    return F_total
end


# --- Error Calculation Function ---

function calculate_error_metrics(field_orig::Matrix{Float64}, field_mapped::Matrix{Float64})
    if size(field_orig) != size(field_mapped)
        error("Field dimensions do not match for error calculation: $(size(field_orig)) vs $(size(field_mapped))")
    end
    if isempty(field_orig)
        return (mae=NaN, rmse=NaN, max_err=NaN, error_vectors=Matrix{Float64}(undef, 3, 0))
    end

    error_vectors = field_orig .- field_mapped
    # Handle potential NaN/Inf in error vectors if interpolation failed badly
    error_magnitudes = [norm(ev) for ev in eachcol(error_vectors)]
    valid_indices = isfinite.(error_magnitudes)

    if !any(valid_indices)
        @warn "All error magnitudes are non-finite. Cannot calculate metrics."
        return (mae=NaN, rmse=NaN, max_err=NaN, error_vectors=error_vectors)
    end

    valid_error_magnitudes = error_magnitudes[valid_indices]

    mae = mean(valid_error_magnitudes)
    rmse = sqrt(mean(valid_error_magnitudes .^ 2))
    max_err = maximum(valid_error_magnitudes)

    return (mae=mae, rmse=rmse, max_err=max_err, error_vectors=error_vectors)
end

# --- Main Evaluation Script Logic ---

function run_evaluation()
    # Create output directory if it doesn't exist
    mkpath(OUTPUT_DIR)

    # --- Load Data ---
    @info "Loading electromagnetic data from $JLD_DATA_DIR..."
    local elMagData # Make it local
    try
        elMagData = create_trees(JLD_DATA_DIR)
    catch e
        @error "Failed to load electromagnetic data using create_trees: $e"
        return
    end

    @info "Loading target mesh cell centers from $TARGET_CASE_DIR..."
    local target_centers # Make it local
    local n_points
    try
        target_centers_path = joinpath(TARGET_CASE_DIR, TARGET_CENTERS_PATH)
        target_centers = read_foam_vector_field(target_centers_path) # These are in OpenFOAM coords
        n_points = size(target_centers, 2)
        @info "Loaded target mesh: $n_points cells."
        if n_points == 0
            @error "Target mesh centers file is empty or could not be read correctly."
            return
        end
    catch e
        @error "Failed to load target mesh centers: $e"
        return
    end

    # --- Experiment: Calculate Force Fields for Theta Combinations ---
    @info "\n--- Running Experiment: Calculate Force Fields for Theta Combinations ---"
    @info "Applying coordinate shift: $COORD_SHIFT_VECTOR"

    # Prepare DataFrame to store all visualization data
    # Store the ORIGINAL OpenFOAM coordinates
    df_viz = DataFrame(
        X = target_centers[1, :],
        Y = target_centers[2, :],
        Z = target_centers[3, :]
    )

    # Generate the 16 theta combinations
    theta_combinations = Vector{Float64}[]
    u_strings = String[]
    w_strings = String[]
    all_force_fields = Dict{String, Matrix{Float64}}() # Store calculated fields

    for u_vec in BASIS_VECTORS
        for w_vec in BASIS_VECTORS
            theta = vcat(u_vec, w_vec)
            u_str = basis_to_string(u_vec)
            w_str = basis_to_string(w_vec)
            push!(theta_combinations, theta)
            push!(u_strings, u_str)
            push!(w_strings, w_str)

            col_prefix = "F_u$(u_str)_w$(w_str)"

            @info "Calculating force field for theta = [$(join(theta, ','))] (u=$(u_str), w=$(w_str))..."
            # Pass the original target_centers and the shift vector
            F_theta = get_force_at_points_k(elMagData, target_centers, theta, INTERPOLATION_K, INTERPOLATION_P, COORD_SHIFT_VECTOR)
            all_force_fields[col_prefix] = F_theta # Store the calculated field

            # Add components and magnitude to the DataFrame
            df_viz[!, "$(col_prefix)_X"] = F_theta[1, :]
            df_viz[!, "$(col_prefix)_Y"] = F_theta[2, :]
            df_viz[!, "$(col_prefix)_Z"] = F_theta[3, :]
            df_viz[!, "$(col_prefix)_Mag"] = [norm(F_theta[:, j]) for j in 1:n_points]
        end
    end
    @info "Calculated $(length(theta_combinations)) force fields."

    # Save the combined visualization data to CSV
    viz_csv_path = joinpath(OUTPUT_DIR, "theta_combinations_force_data.csv")
    @info "Saving visualization data for all theta combinations to $viz_csv_path"
    CSV.write(viz_csv_path, df_viz)

    # --- Calculate and Save Quantitative Metrics ---
    @info "\n--- Calculating Quantitative Metrics (vs. u=$(REF_U_STR), w=$(REF_W_STR)) ---"
    metrics_results_list = []
    if haskey(all_force_fields, REF_COL_PREFIX)
        F_ref = all_force_fields[REF_COL_PREFIX]
        @info "Using reference field: $REF_COL_PREFIX"

        for i in 1:length(theta_combinations)
            u_str = u_strings[i]
            w_str = w_strings[i]
            col_prefix = "F_u$(u_str)_w$(w_str)"

            # Skip comparing the reference field to itself
            if col_prefix == REF_COL_PREFIX
                push!(metrics_results_list, Dict(
                    "u" => u_str, "w" => w_str,
                    "mae" => 0.0, "rmse" => 0.0, "max_err" => 0.0
                ))
                continue
            end

            if haskey(all_force_fields, col_prefix)
                F_current = all_force_fields[col_prefix]
                metrics = calculate_error_metrics(F_ref, F_current)
                push!(metrics_results_list, Dict(
                    "u" => u_str, "w" => w_str,
                    "mae" => metrics.mae, "rmse" => metrics.rmse, "max_err" => metrics.max_err
                ))
                @info @sprintf("Metrics for u=%s, w=%s (vs ref): MAE=%.4g, RMSE=%.4g, MaxErr=%.4g", u_str, w_str, metrics.mae, metrics.rmse, metrics.max_err)
            else
                 @warn "Could not find calculated field '$col_prefix' for metrics calculation."
            end
        end

        # Save metrics to CSV
        if !isempty(metrics_results_list)
            df_metrics = DataFrame(metrics_results_list)
            metrics_csv_path = joinpath(OUTPUT_DIR, "theta_combinations_metrics.csv")
            CSV.write(metrics_csv_path, df_metrics)
            @info "Quantitative metrics saved to $metrics_csv_path"
        end
    else
         @error "Reference field '$REF_COL_PREFIX' was not found in calculated fields. Cannot calculate metrics."
    end


    # --- Plot Generation (REMOVED) ---

    # --- Experiment 2: Smoothness Visualization (REMOVED - Use CSV) ---
    @info "\n--- Smoothness Visualization ---"
    @info "Visualization data saved in theta_combinations_force_data.csv"
    @info "Load this CSV into ParaView, use 'Table To Points' filter, and visualize desired F_u..._w..._Mag/vectors."


    @info "\nEvaluation script finished."
end

# --- Run the Evaluation ---
run_evaluation()
