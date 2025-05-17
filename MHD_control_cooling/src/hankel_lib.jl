#=
This script performs Hankel Alternative View of Koopman (HAVOK) analysis
on a single time-series observable from X_data, incorporating external
control inputs U_data, and produces augmented state-space matrices
A_full and B_full.

The model for the primary HAVOK coordinates v'_k is:
v'_{k+1} = sum(A_prime_i * v'_{k-i}) + sum(B_vr_i * v_{r,k-i}) + sum(C_ext_i * U_ext_{k-i})

A_full and B_full operate on an augmented state Z_k including delayed v', v_r, and U_ext.
The evolution of v_r itself is not dynamically modeled by A_full beyond shifting past values.

It also includes functionality to simulate the identified models and plot comparisons.
=#
using JLD2
using LinearAlgebra
using Statistics
using Plots # For plotting

# --- Configuration ---
const DATA_PATH = "./data.jld2" # Path to your JLD2 file
const RELEVANT_X_STATE_FOR_HANKEL = [6] # Index of the state in X_data for Hankel matrix
# const ALL_RELEVANT_X_STATES = [6] # Not directly used in this version, but kept for context

const HANKEL_Q = 100      # Number of rows for the Hankel matrix (embedding dimension q for the single observable)
const HANKEL_P_MAX = 3000 # Max number of columns for the Hankel matrix (snapshots p)
const HAVOK_R = 15        # Rank for SVD truncation / total number of HAVOK coordinates (r_H)
# Ensure HAVOK_R >= 2

const DELAY_COMMON = 60   # Common delay D for v', v_r, and U_ext in the DMDc-like regression
# and for constructing the augmented state Z_k.

const OUTPUT_FILENAME = "havok_dmdc_AB_matrices.jld2"
const NORMALIZE_X_FOR_HANKEL = true # Whether to center the single timeseries before Hankel
const NORMALIZE_OMEGA_INPUTS = true # Whether to center columns of Omega matrix in regression
const REGULARIZATION_LAMBDA = 1e-6 # Small regularization factor for pinv(Omega)

# --- Plotting Configuration ---
const START_INDICES_PLOT = [500, 1000, 1500, 2000] # Start indices for simulation plots (within havok_v_coords_cols)
const SIM_LEN_PLOT = 200          # Length of each simulation run for plotting
# Index (1-based) of the v' component to plot.
# E.g., if HAVOK_R=15, N_v_prime=14. PLOT_VP_COMPONENT_IDX=1 plots the first component of v'.
const PLOT_VP_COMPONENT_IDX = 1
const OUTPUT_DIR_PLOTS = "./figures/comparison_plots_julia" # Directory to save plots

# Plotting Backend and Style
gr()
default(
    fontfamily="Computer Modern",
    linewidth=1.5,
    markersize=3,
    legendfontsize=8,
    tickfontsize=8,
    guidefontsize=10,
    titlefontsize=10,
    dpi=300,
    grid=true,
    framestyle=:box
)

# --- Helper and Core Functions ---

"""
Loads X and U data, preprocesses U.
Extracts single time series for Hankel.
"""
function load_and_preprocess_data_havok(data_path::String)
    println("Loading data from $data_path...")
    local X_data_full, U_data_full
    try
        file = jldopen(data_path, "r")
        X_data_full = read(file, "X_data")
        U_data_full = read(file, "U_data")
        close(file)
    catch e
        println("Error loading data from $data_path: $e")
        rethrow(e)
    end

    if length(RELEVANT_X_STATE_FOR_HANKEL) != 1
        error("RELEVANT_X_STATE_FOR_HANKEL must specify exactly one state for Hankel construction.")
    end
    hankel_series_idx = RELEVANT_X_STATE_FOR_HANKEL[1]

    timeseries_for_hankel = vec(X_data_full[hankel_series_idx, 2:end])
    println("Extracted time series for Hankel (state $hankel_series_idx) with $(length(timeseries_for_hankel)) samples.")

    println("Preprocessing U_data...")
    U_ext_processed = let
        cols = []
        for u_col in eachcol(U_data_full)
            if size(u_col, 1) == 8
                push!(cols, vec(u_col[1:4] * u_col[5:8]'))
            elseif size(u_col, 1) == 16
                push!(cols, vec(u_col))
            else
                println("Warning: Unexpected size for column in U_data_full. Size: $(size(u_col)).")
                if ndims(u_col) == 1
                    push!(cols, u_col)
                else
                    error("Unhandled U_data column structure.")
                end
            end
        end
        hcat(cols...)
    end
    N_U_ext = size(U_ext_processed, 1)
    println("Processed U_ext_data. Size: $(size(U_ext_processed)) -> N_U_ext = $N_U_ext")

    min_total_samples = min(length(timeseries_for_hankel), size(U_ext_processed, 2))
    timeseries_for_hankel = timeseries_for_hankel[1:min_total_samples]
    U_ext_aligned = U_ext_processed[:, 1:min_total_samples]
    println("Aligned Hankel series and U_ext to $min_total_samples samples.")

    return timeseries_for_hankel, U_ext_aligned, N_U_ext, min_total_samples
end

"""
Constructs a Hankel matrix from a time series.
"""
function construct_hankel_matrix(series::AbstractVector, q::Int, p::Int; center::Bool=false)
    println("Constructing Hankel matrix ($q x $p)...")
    n_series = length(series)
    if q <= 0 || p <= 0
        error("Hankel dimensions q and p must be positive.")
    end
    if n_series < q + p - 1
        error("Time series too short. Needs $(q + p - 1) points, got $n_series.")
    end

    series_to_embed = copy(series)
    if center
        series_to_embed .-= mean(series_to_embed)
        println("Centered series for Hankel.")
    end

    H = similar(series_to_embed, q, p)
    for j in 1:p
        for i in 1:q
            H[i, j] = series_to_embed[i+j-1]
        end
    end
    println("Hankel matrix constructed.")
    return H
end

"""
Performs SVD on Hankel matrix, extracts r_H HAVOK coordinates.
"""
function get_havok_coordinates(H::AbstractMatrix, r_havok::Int)
    println("Performing SVD on Hankel matrix (size $(size(H)))...")
    if r_havok < 2
        error("HAVOK_R (r_havok) must be at least 2.")
    end

    U_h, S_h, V_h_full = svd(H)

    actual_r_havok = min(r_havok, size(V_h_full, 2), size(U_h, 2))
    if actual_r_havok != r_havok
        println("Warning: HAVOK_R ($r_havok) adjusted to $actual_r_havok due to SVD rank.")
    end
    if actual_r_havok < 2
        error("Effective HAVOK_R ($actual_r_havok) is less than 2.")
    end

    havok_coords_v_cols = V_h_full[:, 1:actual_r_havok]

    println("HAVOK coordinates (from V_h columns) extracted. Size: $(size(havok_coords_v_cols)), r_H = $actual_r_havok")
    return havok_coords_v_cols, S_h[1:actual_r_havok], actual_r_havok
end

"""
Creates a delay-embedded matrix from a multivariate time series matrix.
"""
function general_delay_embed(data::AbstractMatrix, delay::Int)
    n_features, n_samples = size(data)
    if delay < 0
        error("Delay must be non-negative.")
    end
    if delay == 0
        return data
    end # No delay, return original
    if n_samples <= delay
        error("Not enough samples ($n_samples) for delay ($delay).")
    end

    n_embedded_samples = n_samples - delay
    embedded_dim = n_features * (delay + 1)
    embedded_data = similar(data, embedded_dim, n_embedded_samples)

    for t in 1:n_embedded_samples
        for d in 0:delay
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features
            source_col = t + delay - d
            embedded_data[row_start:row_end, t] = data[:, source_col]
        end
    end
    return embedded_data
end

"""
Calculates HAVOK-DMDc components: A_prime_i, B_vr_i, C_ext_i
"""
function calculate_havok_dmdc_components(
    havok_v_coords_cols_arg::AbstractMatrix,
    U_ext_data_cols::AbstractMatrix,
    r_H::Int,
    N_U_ext::Int,
    delay_common::Int;
    lambda_reg::Float64=1e-6 # Regularization parameter
)
    println("Calculating HAVOK-DMDc components with common delay D = $delay_common, regularization lambda = $lambda_reg...")
    havok_v_coords_rows = permutedims(havok_v_coords_cols_arg, (2, 1))

    v_prime_series_rows = havok_v_coords_rows[1:r_H-1, :]
    v_r_series_rows = havok_v_coords_rows[r_H:r_H, :]

    N_v_prime = r_H - 1
    N_v_r = 1

    num_snapshots_available = size(havok_v_coords_rows, 2)

    if num_snapshots_available <= delay_common + 1
        error("Not enough snapshots in HAVOK coordinates ($num_snapshots_available) for common delay $delay_common and predicting next step.")
    end

    Y_target_v_prime = v_prime_series_rows[:, delay_common+2:num_snapshots_available]

    num_reg_samples = num_snapshots_available - (delay_common + 1)

    v_prime_for_embed = v_prime_series_rows[:, 1:num_snapshots_available-1]
    v_r_for_embed = v_r_series_rows[:, 1:num_snapshots_available-1]
    if size(U_ext_data_cols, 2) < (num_snapshots_available - 1)
        error("U_ext_data_cols has insufficient length ($(size(U_ext_data_cols,2))) for embedding. Needs at least $(num_snapshots_available-1).")
    end
    U_ext_for_embed = U_ext_data_cols[:, 1:num_snapshots_available-1]

    X_embed_v_prime = general_delay_embed(v_prime_for_embed, delay_common)
    X_embed_v_r = general_delay_embed(v_r_for_embed, delay_common)
    X_embed_U_ext = general_delay_embed(U_ext_for_embed, delay_common)

    Omega = [X_embed_v_prime; X_embed_v_r; X_embed_U_ext]
    local Omega_mean
    if NORMALIZE_OMEGA_INPUTS
        Omega_mean = mean(Omega, dims=2)
        Omega = Omega .- Omega_mean
        println("Omega matrix centered (mean subtracted column-wise).")
    else
        Omega_mean = zeros(size(Omega, 1), 1)
    end

    println("Performing regression for M_havok_dmdc...")
    println("Size of Y_target_v_prime: $(size(Y_target_v_prime))")
    println("Size of Omega: $(size(Omega))")

    if size(Y_target_v_prime, 2) != size(Omega, 2)
        error("Mismatch in number of samples for regression: Y_target has $(size(Y_target_v_prime, 2)), Omega has $(size(Omega, 2))")
    end

    # Regularized pseudo-inverse (Ridge Regression)
    # M = Y * Omega' * inv(Omega*Omega' + lambda*I)
    # This form is for when Omega has more columns than rows (n_features_omega < n_samples_omega)
    # Here, Omega rows = features, Omega cols = samples.
    # size(Omega) is (total_embedded_dim, num_reg_samples)
    # Your Omega is (1891, 2939), so more columns than rows.
    if size(Omega, 2) >= size(Omega, 1) # More samples than features (or equal)
        M_havok_dmdc = Y_target_v_prime * (Omega' * inv(Omega * Omega' + lambda_reg * I(size(Omega, 1))))
        println("Used Omega' * inv(Omega*Omega' + lambda*I) for M calculation.")
    else # More features than samples
        M_havok_dmdc = Y_target_v_prime * (inv(Omega' * Omega + lambda_reg * I(size(Omega, 2))) * Omega')
        println("Used inv(Omega'*Omega + lambda*I) * Omega' for M calculation.")
    end
    # M_havok_dmdc = Y_target_v_prime * pinv(Omega) # Original
    println("M_havok_dmdc matrix size: $(size(M_havok_dmdc))")

    A_prime_list = Vector{Matrix{Float64}}(undef, delay_common + 1)
    B_vr_list = Vector{Matrix{Float64}}(undef, delay_common + 1)
    C_ext_list = Vector{Matrix{Float64}}(undef, delay_common + 1)

    current_col = 1
    for i in 0:delay_common
        A_prime_list[i+1] = M_havok_dmdc[:, current_col:current_col+N_v_prime-1]
        current_col += N_v_prime
    end
    for i in 0:delay_common
        B_vr_list[i+1] = M_havok_dmdc[:, current_col:current_col+N_v_r-1]
        current_col += N_v_r
    end
    for i in 0:delay_common
        C_ext_list[i+1] = M_havok_dmdc[:, current_col:current_col+N_U_ext-1]
        current_col += N_U_ext
    end

    println("Extracted component matrices for HAVOK-DMDc.")
    return A_prime_list, B_vr_list, C_ext_list, M_havok_dmdc, Omega_mean
end

"""
Constructs augmented state-space matrices A_full, B_full for HAVOK-DMDc.
"""
function construct_havok_augmented_matrices(
    A_prime_list::Vector, B_vr_list::Vector, C_ext_list::Vector,
    r_H::Int, N_U_ext::Int, delay_common::Int
)
    println("Constructing augmented A_full, B_full for HAVOK-DMDc...")
    D = delay_common
    N_v_prime = r_H - 1
    N_v_r = 1

    dim_Z_v_prime = N_v_prime * (D + 1)
    dim_Z_v_r = N_v_r * (D + 1)
    dim_Z_U_ext = N_U_ext * D
    N_AUG_STATE = dim_Z_v_prime + dim_Z_v_r + dim_Z_U_ext

    A_full = zeros(N_AUG_STATE, N_AUG_STATE)
    B_full = zeros(N_AUG_STATE, N_U_ext)

    B_full[1:N_v_prime, :] = C_ext_list[1]

    if D > 0
        offset_U_in_Z = dim_Z_v_prime + dim_Z_v_r
        B_full[offset_U_in_Z+1:offset_U_in_Z+N_U_ext, :] = I(N_U_ext)
    end

    col_offset = 0
    for i in 0:D
        A_full[1:N_v_prime, col_offset+(i*N_v_prime)+1:col_offset+((i+1)*N_v_prime)] = A_prime_list[i+1]
    end
    col_offset += dim_Z_v_prime
    for i in 0:D
        A_full[1:N_v_prime, col_offset+(i*N_v_r)+1:col_offset+((i+1)*N_v_r)] = B_vr_list[i+1]
    end
    col_offset += dim_Z_v_r
    if D > 0
        for i in 1:D
            A_full[1:N_v_prime, col_offset+((i-1)*N_U_ext)+1:col_offset+(i*N_U_ext)] = C_ext_list[i+1]
        end
    end

    row_offset_v_prime = N_v_prime
    col_offset_v_prime = 0
    if D > 0
        A_full[row_offset_v_prime+1:dim_Z_v_prime, col_offset_v_prime+1:dim_Z_v_prime-N_v_prime] = I(N_v_prime * D)
    end

    row_offset_v_r = dim_Z_v_prime + N_v_r
    col_offset_v_r = dim_Z_v_prime
    if D > 0
        A_full[row_offset_v_r+1:dim_Z_v_prime+dim_Z_v_r, col_offset_v_r+1:dim_Z_v_prime+dim_Z_v_r-N_v_r] = I(N_v_r * D)
    end

    if D > 1
        dest_row_start_U_shift = dim_Z_v_prime + dim_Z_v_r + N_U_ext + 1
        dest_row_end_U_shift = N_AUG_STATE

        src_col_start_U_shift = dim_Z_v_prime + dim_Z_v_r + 1
        src_col_end_U_shift = dim_Z_v_prime + dim_Z_v_r + N_U_ext * (D - 1)

        A_full[dest_row_start_U_shift:dest_row_end_U_shift, src_col_start_U_shift:src_col_end_U_shift] = I(N_U_ext * (D - 1))
    end

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end

# --- Simulation and Plotting Functions ---

"""
Creates initial augmented state Z_0 for HAVOK-DMDc simulation.
"""
function create_havok_augmented_init_state(
    havok_v_coords_cols::AbstractMatrix,
    U_ext_data_cols::AbstractMatrix,
    sim_start_idx_in_v_coords::Int,
    delay_common::Int,
    r_H::Int,
    N_U_ext::Int
)
    D = delay_common
    N_v_prime = r_H - 1
    N_v_r = 1

    dim_Z_v_prime = N_v_prime * (D + 1)
    dim_Z_v_r = N_v_r * (D + 1)
    dim_Z_U_ext = N_U_ext * D
    N_AUG_STATE = dim_Z_v_prime + dim_Z_v_r + dim_Z_U_ext

    Z0 = zeros(N_AUG_STATE)

    min_v_idx = sim_start_idx_in_v_coords - D
    if min_v_idx < 1
        error("Cannot create initial state: sim_start_idx_in_v_coords ($sim_start_idx_in_v_coords) too early for delay_common ($D). Needs history back to index $min_v_idx.")
    end
    if sim_start_idx_in_v_coords > size(havok_v_coords_cols, 1)
        error("sim_start_idx_in_v_coords ($sim_start_idx_in_v_coords) exceeds available HAVOK coordinate samples ($(size(havok_v_coords_cols,1))).")
    end

    for i in 0:D
        row_start = i * N_v_prime + 1
        row_end = (i + 1) * N_v_prime
        Z0[row_start:row_end] = havok_v_coords_cols[sim_start_idx_in_v_coords-i, 1:N_v_prime]
    end

    offset_vr_in_Z = dim_Z_v_prime
    for i in 0:D
        row_start = offset_vr_in_Z + i * N_v_r + 1
        row_end = offset_vr_in_Z + (i + 1) * N_v_r
        Z0[row_start:row_end] = [havok_v_coords_cols[sim_start_idx_in_v_coords-i, r_H]]
    end

    if D > 0
        offset_U_in_Z = dim_Z_v_prime + dim_Z_v_r
        min_u_idx = sim_start_idx_in_v_coords - D # This is the index for U_ext_{k-D}
        if min_u_idx < 1 # sim_start_idx_in_v_coords is k. U_ext_{k-1} is at sim_start_idx_in_v_coords-1
            error("Cannot create initial state for U_ext: sim_start_idx_in_v_coords ($sim_start_idx_in_v_coords) too early for delay_common ($D) for U_ext. Needs U_ext data back to index $(sim_start_idx_in_v_coords-D).")
        end

        for i in 1:D # For U_ext_{k-1} down to U_ext_{k-D}
            row_start = offset_U_in_Z + (i - 1) * N_U_ext + 1
            row_end = offset_U_in_Z + i * N_U_ext
            # U_ext_data_cols is aligned with havok_v_coords_cols time
            # So U_ext_{k-i} is at index (sim_start_idx_in_v_coords - i)
            if (sim_start_idx_in_v_coords - i) < 1 || (sim_start_idx_in_v_coords - i) > size(U_ext_data_cols, 2)
                error("Index out of bounds for U_ext_data_cols: trying to access $(sim_start_idx_in_v_coords - i)")
            end
            Z0[row_start:row_end] = U_ext_data_cols[:, sim_start_idx_in_v_coords-i]
        end
    end
    return Z0
end

"""
Simulates the HAVOK augmented state-space model.
"""
function simulate_havok_augmented_model(
    A_full::AbstractMatrix, B_full::AbstractMatrix,
    Z0::AbstractVector,
    U_ext_segment::AbstractMatrix,
    sim_len::Int
)
    N_AUG_STATE = size(A_full, 1)
    states_Z_augmented = zeros(N_AUG_STATE, sim_len)
    current_Z = copy(Z0)

    if size(U_ext_segment, 2) < sim_len
        error("Not enough U_ext inputs for simulation length.")
    end

    for k in 1:sim_len
        try
            current_Z = A_full * current_Z + B_full * U_ext_segment[:, k]
            states_Z_augmented[:, k] = current_Z
        catch e
            println("Error during augmented simulation at step $k: $e")
            println("A_full size: $(size(A_full)), current_Z size: $(size(current_Z))")
            println("B_full size: $(size(B_full)), U_ext_segment[:,k] size: $(size(U_ext_segment[:, k]))")
            states_Z_augmented[:, k:end] .= NaN # Mark rest as NaN
            break
        end
        if any(isnan, current_Z) || any(isinf, current_Z)
            println("NaN or Inf detected in augmented simulation at step $k. Stopping.")
            states_Z_augmented[:, k:end] .= NaN
            break
        end
    end
    return states_Z_augmented
end

"""
Creates initial Omega vector for the M-Regression model simulation.
"""
function create_havok_M_model_init_omega(
    havok_v_coords_cols::AbstractMatrix,
    U_ext_data_cols::AbstractMatrix,
    sim_start_idx_in_v_coords::Int,
    delay_common::Int,
    r_H::Int, N_U_ext::Int,
    Omega_mean_val::AbstractVecOrMat
)
    D = delay_common
    N_v_prime = r_H - 1
    N_v_r = 1

    min_hist_idx = sim_start_idx_in_v_coords - D
    if min_hist_idx < 1
        error("Cannot create initial Omega: sim_start_idx ($sim_start_idx_in_v_coords) too early for delay ($D).")
    end

    v_prime_segment = permutedims(havok_v_coords_cols[min_hist_idx:sim_start_idx_in_v_coords, 1:N_v_prime], (2, 1))
    v_r_segment = permutedims(havok_v_coords_cols[min_hist_idx:sim_start_idx_in_v_coords, r_H:r_H], (2, 1))
    U_ext_segment_hist = U_ext_data_cols[:, min_hist_idx:sim_start_idx_in_v_coords]

    omega_v_prime = general_delay_embed(v_prime_segment, D) # Should be (N_v_prime*(D+1)) x 1
    omega_v_r = general_delay_embed(v_r_segment, D)     # Should be (N_v_r*(D+1)) x 1
    omega_U_ext = general_delay_embed(U_ext_segment_hist, D) # Should be (N_U_ext*(D+1)) x 1

    Omega0 = [omega_v_prime; omega_v_r; omega_U_ext]
    if NORMALIZE_OMEGA_INPUTS
        Omega0 .-= Omega_mean_val
    end
    return vec(Omega0)
end

"""
Simulates the M-Regression model.
"""
function simulate_havok_M_model(
    Omega0::AbstractVector,
    U_ext_segment_sim::AbstractMatrix,
    havok_v_r_truth_segment_sim::AbstractMatrix,
    M_regression::AbstractMatrix,
    sim_len::Int, delay_common::Int,
    r_H::Int, N_U_ext::Int,
    Omega_mean_val::AbstractVecOrMat
)
    D = delay_common
    N_v_prime = r_H - 1
    N_v_r = 1

    predicted_v_prime = zeros(N_v_prime, sim_len)
    current_Omega = copy(Omega0)

    history_vp = zeros(N_v_prime, D + 1)
    history_vr = zeros(N_v_r, D + 1)
    history_U = zeros(N_U_ext, D + 1)

    ptr = 1
    for i in 0:D
        history_vp[:, D-i+1] = current_Omega[ptr:ptr+N_v_prime-1]
        ptr += N_v_prime
    end
    for i in 0:D
        history_vr[:, D-i+1] = current_Omega[ptr:ptr+N_v_r-1]
        ptr += N_v_r
    end
    for i in 0:D
        history_U[:, D-i+1] = current_Omega[ptr:ptr+N_U_ext-1]
        ptr += N_U_ext
    end

    for k in 1:sim_len
        vp_next = M_regression * current_Omega
        predicted_v_prime[:, k] = vp_next

        if any(isnan, vp_next) || any(isinf, vp_next)
            println("NaN or Inf detected in M-model simulation at step $k. Stopping.")
            predicted_v_prime[:, k:end] .= NaN
            break
        end

        if k < sim_len
            history_vp = circshift(history_vp, (0, 1))
            history_vr = circshift(history_vr, (0, 1))
            history_U = circshift(history_U, (0, 1))

            history_vp[:, 1] = vp_next
            history_vr[:, 1] = havok_v_r_truth_segment_sim[:, k+1]
            history_U[:, 1] = U_ext_segment_sim[:, k+1]

            ptr = 1
            for i in 0:D
                current_Omega[ptr:ptr+N_v_prime-1] = history_vp[:, D-i+1]
                ptr += N_v_prime
            end
            for i in 0:D
                current_Omega[ptr:ptr+N_v_r-1] = history_vr[:, D-i+1]
                ptr += N_v_r
            end
            for i in 0:D
                current_Omega[ptr:ptr+N_U_ext-1] = history_U[:, D-i+1]
                ptr += N_U_ext
            end

            if NORMALIZE_OMEGA_INPUTS
                current_Omega .-= Omega_mean_val
            end
        end
    end
    return predicted_v_prime
end

"""
Plots a single comparison onto a specified subplot.
"""
function plot_single_havok_comparison_subplot!(
    p_obj::Plots.Plot, subplot_idx::Int,
    actual_vp_component_series::AbstractVector,
    sim_vp_component_augmented::AbstractVector,
    sim_vp_component_M_model::AbstractVector,
    sim_start_time_label::Int,
    plot_vp_comp_idx_label::Int
)
    sim_len = length(actual_vp_component_series)
    time_steps = 1:sim_len

    plot!(p_obj[subplot_idx], time_steps, actual_vp_component_series,
        label="Actual v'_{$(plot_vp_comp_idx_label)}", color=:black, linewidth=2.0)
    plot!(p_obj[subplot_idx], time_steps, sim_vp_component_augmented,
        label="Aug.Model v'_{$(plot_vp_comp_idx_label)}", color=:red, linestyle=:dash, linewidth=1.2)
    plot!(p_obj[subplot_idx], time_steps, sim_vp_component_M_model,
        label="M-Reg.Model v'_{$(plot_vp_comp_idx_label)}", color=:blue, linestyle=:dot, linewidth=1.2)

    title!(p_obj[subplot_idx], "Start Idx = $sim_start_time_label")
    xlabel!(p_obj[subplot_idx], "Time Step")
    if subplot_idx % Int(ceil(sqrt(length(p_obj.layout.grid)))) == 1 || length(p_obj.layout.grid) <= 2
        ylabel!(p_obj[subplot_idx], "v'_{$(plot_vp_comp_idx_label)} Value")
    end
    plot!(p_obj[subplot_idx], legend=:best)
end


# --- Main Execution ---
function main()
    println("Starting HAVOK-DMDc Analysis Script with Plotting...")

    timeseries_hankel, U_ext_data, N_U_ext, total_aligned_samples = load_and_preprocess_data_havok(DATA_PATH)

    current_hankel_p = min(HANKEL_P_MAX, total_aligned_samples - HANKEL_Q + 1)
    if current_hankel_p <= 0
        error("Not enough samples for Hankel matrix.")
    end
    if current_hankel_p < HAVOK_R
        error("Adjusted HANKEL_P < HAVOK_R.")
    end
    println("Using HANKEL_P = $current_hankel_p")

    H = construct_hankel_matrix(timeseries_hankel, HANKEL_Q, current_hankel_p, center=NORMALIZE_X_FOR_HANKEL)

    effective_r_H = min(HAVOK_R, HANKEL_Q, current_hankel_p)
    if HAVOK_R != effective_r_H
        println("HAVOK_R adjusted to $effective_r_H")
    end
    if effective_r_H < 2
        error("Effective HAVOK_R ($effective_r_H) < 2.")
    end

    havok_v_coords_cols, S_h_values, actual_r_H = get_havok_coordinates(H, effective_r_H)
    N_v_prime = actual_r_H - 1

    if size(U_ext_data, 2) < current_hankel_p
        error("U_ext_data has fewer samples than HANKEL_P.")
    end
    U_ext_for_reg_and_sim = U_ext_data[:, 1:current_hankel_p]

    A_prime_coeffs, B_vr_coeffs, C_ext_coeffs, M_reg, Omega_mean_val = calculate_havok_dmdc_components(
        havok_v_coords_cols,
        U_ext_for_reg_and_sim,
        actual_r_H,
        N_U_ext,
        DELAY_COMMON,
        lambda_reg=REGULARIZATION_LAMBDA
    )

    A_full_final, B_full_final = construct_havok_augmented_matrices(
        A_prime_coeffs, B_vr_coeffs, C_ext_coeffs,
        actual_r_H, N_U_ext, DELAY_COMMON
    )

    println("\n--- HAVOK-DMDc Analysis Complete ---")
    # Eigenvalue analysis of A_full_final
    try
        eigen_A_full = eigen(A_full_final)
        max_abs_eig = maximum(abs.(eigen_A_full.values))
        println("Maximum absolute eigenvalue of A_full_final: $max_abs_eig")
        if max_abs_eig > 1.0
            println("WARNING: A_full_final is unstable (max abs eigenvalue > 1). This likely explains simulation blow-up.")
        else
            println("A_full_final appears stable (max abs eigenvalue <= 1).")
        end
    catch e
        println("Could not compute eigenvalues of A_full_final: $e")
    end


    jldsave(OUTPUT_FILENAME;
        A_full=A_full_final, B_full=B_full_final,
        A_prime_list=A_prime_coeffs, B_vr_list=B_vr_coeffs, C_ext_list=C_ext_coeffs,
        M_regression=M_reg, Omega_mean=Omega_mean_val,
        HANKEL_Q_used=HANKEL_Q, HANKEL_P_used=current_hankel_p,
        HAVOK_R_used=actual_r_H, DELAY_COMMON_used=DELAY_COMMON,
        singular_values_hankel=S_h_values,
        N_U_ext_val=N_U_ext, N_v_prime_val=N_v_prime,
        RELEVANT_X_STATE_FOR_HANKEL_val=RELEVANT_X_STATE_FOR_HANKEL
    )
    println("\nResults saved to $OUTPUT_FILENAME")
    println("A_full size: $(size(A_full_final)), B_full size: $(size(B_full_final))")

    # --- Plotting Section ---
    println("\n--- Starting Plotting Section ---")
    if PLOT_VP_COMPONENT_IDX < 1 || PLOT_VP_COMPONENT_IDX > N_v_prime
        @error "Invalid PLOT_VP_COMPONENT_IDX ($PLOT_VP_COMPONENT_IDX). Must be between 1 and $N_v_prime. Skipping plots."
        return
    end

    num_plots = length(START_INDICES_PLOT)
    grid_cols = ceil(Int, sqrt(num_plots))
    grid_rows = ceil(Int, num_plots / grid_cols)
    plot_layout = (grid_rows, grid_cols)
    fig_width = 400 * grid_cols
    fig_height = 300 * grid_rows
    combined_plot_obj = plot(layout=plot_layout, size=(fig_width, fig_height), legend=false)
    subplot_counter = 0

    for sim_start_idx in START_INDICES_PLOT # These indices are into havok_v_coords_cols (1 to HANKEL_P)
        subplot_counter += 1
        println("\nPlotting for simulation start index in v_coords: $sim_start_idx")

        # Determine valid simulation length
        # Max length for actual_vp_data, and for U_ext_sim_segment, and for v_r_truth_for_M_sim
        max_len_from_v_coords = size(havok_v_coords_cols, 1) - sim_start_idx + 1
        max_len_from_U_ext = size(U_ext_for_reg_and_sim, 2) - sim_start_idx + 1

        current_sim_len = min(SIM_LEN_PLOT, max_len_from_v_coords, max_len_from_U_ext)

        # Check if enough history for initial states
        if (sim_start_idx - DELAY_COMMON < 1)
            println("Warning: Not enough history for initial state at sim_start_idx $sim_start_idx with delay $DELAY_COMMON. Skipping.")
            annotate!(combined_plot_obj[subplot_counter], (0.5, 0.5), text("Skipped (init hist)", :orange, :center, 8))
            title!(combined_plot_obj[subplot_counter], "Start = $sim_start_idx")
            continue
        end
        # Check if simulation length is positive
        if current_sim_len <= 0
            println("Warning: Non-positive simulation length ($current_sim_len) for start_idx $sim_start_idx. Skipping.")
            annotate!(combined_plot_obj[subplot_counter], (0.5, 0.5), text("Skipped (sim len <=0)", :orange, :center, 8))
            title!(combined_plot_obj[subplot_counter], "Start = $sim_start_idx")
            continue
        end
        # Check if M-model simulation needs data beyond v_coords for v_r_truth
        if (sim_start_idx + current_sim_len > size(havok_v_coords_cols, 1)) # M-model needs one step ahead for v_r and U_ext
            current_sim_len = size(havok_v_coords_cols, 1) - sim_start_idx # Adjust to not read past v_r for M-model
            if current_sim_len <= 0
                println("Warning: Adjusted sim_len for M-model is non-positive. Skipping.")
                annotate!(combined_plot_obj[subplot_counter], (0.5, 0.5), text("Skipped (M-model bound)", :orange, :center, 8))
                title!(combined_plot_obj[subplot_counter], "Start = $sim_start_idx")
                continue
            end
            println("Adjusted sim_len to $current_sim_len for M-model v_r truth bound.")
        end


        actual_vp_data = havok_v_coords_cols[sim_start_idx:sim_start_idx+current_sim_len-1, PLOT_VP_COMPONENT_IDX]
        U_ext_sim_segment = U_ext_for_reg_and_sim[:, sim_start_idx:sim_start_idx+current_sim_len-1]

        Z0_aug = create_havok_augmented_init_state(havok_v_coords_cols, U_ext_for_reg_and_sim, sim_start_idx, DELAY_COMMON, actual_r_H, N_U_ext)
        sim_Z_augmented = simulate_havok_augmented_model(A_full_final, B_full_final, Z0_aug, U_ext_sim_segment, current_sim_len)
        sim_vp_augmented = sim_Z_augmented[PLOT_VP_COMPONENT_IDX, :]

        Omega0_M = create_havok_M_model_init_omega(havok_v_coords_cols, U_ext_for_reg_and_sim, sim_start_idx, DELAY_COMMON, actual_r_H, N_U_ext, Omega_mean_val)
        # v_r truth for M model sim: needs v_r(k+1) when predicting v'(k+1) from Omega_k
        # So, if Omega_k uses data up to time t (sim_start_idx), M-model needs v_r from t up to t+sim_len
        # The loop in simulate_havok_M_model accesses havok_v_r_truth_segment_sim[:, k+1]
        # So it needs to be of length sim_len, corresponding to v_r at sim_start_idx+1 to sim_start_idx+sim_len
        v_r_truth_for_M_sim = permutedims(havok_v_coords_cols[sim_start_idx:sim_start_idx+current_sim_len, actual_r_H:actual_r_H], (2, 1)) # Length sim_len+1 for k+1 access

        sim_vp_M_model_all_comps = simulate_havok_M_model(Omega0_M, U_ext_sim_segment, v_r_truth_for_M_sim, M_reg, current_sim_len, DELAY_COMMON, actual_r_H, N_U_ext, Omega_mean_val)
        sim_vp_M_model = sim_vp_M_model_all_comps[PLOT_VP_COMPONENT_IDX, :]

        plot_single_havok_comparison_subplot!(
            combined_plot_obj, subplot_counter,
            actual_vp_data,
            sim_vp_augmented,
            sim_vp_M_model,
            sim_start_idx,
            PLOT_VP_COMPONENT_IDX
        )
    end

    plot!(combined_plot_obj, plot_title="HAVOK-DMDc Model Comparison (v'_{$(PLOT_VP_COMPONENT_IDX)}, r_H=$actual_r_H, D=$DELAY_COMMON)",
        plot_titlefontsize=12, top_margin=10Plots.mm)

    if !isdir(OUTPUT_DIR_PLOTS)
        println("Creating output directory: $OUTPUT_DIR_PLOTS")
        mkpath(OUTPUT_DIR_PLOTS)
    end
    figname = joinpath(OUTPUT_DIR_PLOTS, "havok_dmdc_compare_rH$(actual_r_H)_D$(DELAY_COMMON)_vp$(PLOT_VP_COMPONENT_IDX).png")
    try
        savefig(combined_plot_obj, figname)
        println("\nCombined plot saved to $figname")
    catch e
        println("\nError saving combined plot $figname: $e")
    end

    println("Script finished.")
end

# Run the main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
