#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays.
It constructs an augmented state-space model (A_full, B_full)
and compares its simulation results against the original DMDc formulation
and the actual data.

Key Functions:
- `load_and_preprocess_data`: Loads data from JLD2, preprocesses control input U.
- `delay_embed`: Creates delay-embedded matrices.
- `solve_fastdmd`: User-provided function to compute matrix using a specific pseudo-inverse method.
- `calculate_dmdc_matrices`: Computes DMDc system matrices (As, Bs, M) from data.
- `construct_augmented_matrices`: Builds the augmented A_full, B_full matrices.
- `create_abc_system_from_data`: Orchestrates data loading, DMDc, and ABC matrix construction.
- `create_abc_init_state`, `create_orig_init_state`: Generate initial states for simulations.
- `simulate_abc_model`, `simulate_orig_model`: Run simulations for both model types.
- `plot_comparison`: Generates and saves a comparison plot.
- `main`: Executes the full workflow.
=#
using JLD2
using LinearAlgebra
using Statistics
using Plots # Using Plots.jl for simplicity. Add Pkg.add("Plots") if needed.


# --- Core Data Processing and Model Creation Functions ---

"""
Loads X and U data, preprocesses U.
"""
function load_and_preprocess_data(data_path::String)
    println("Loading data from $data_path...")
    try
        X_data = jldopen(data_path, "r")["X_data"]
        X_data = X_data[RELEVANT_STATES, 2:end] .* 1e6
        # X_data = jldopen("./openfoam_surface_data_processed_N20xM20.jld2")["sampled_temperature_matrix_on_grid"]

        U_data = jldopen(data_path, "r")["U_data"]
        # U_data = jldopen("./openfoam_surface_data_processed_N20xM20.jld2")["experimental_input_matrix"]

        println("Preprocessing U data...")
        U = let
            cols = []
            for u_col in eachcol(U_data)
                if size(u_col, 1) == 8
                    push!(cols, vec(u_col[1:4] * u_col[5:8]'))
                else
                    println("Warning: Unexpected size for column in U_data. Skipping.")
                end
            end
            hcat(cols...)
        end # let

        N_X = size(X_data, 1)
        N_U = size(U, 1)
        TOTAL_SAMPLES = size(X_data, 2)
        println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")

        return X_data, U, N_X, N_U, TOTAL_SAMPLES
    catch e
        println("Error loading data from $data_path: $e")
        println("Please ensure the JLD2 file exists and the path is correct.")
        rethrow(e) # Propagate the error
    end # try-catch
end # function load_and_preprocess_data

"""
Creates a delay-embedded matrix from time series data.
"""
function delay_embed(data::AbstractMatrix, delay::Int)
    n_features, n_samples = size(data)
    if delay < 0
        error("Delay must be non-negative.")
    elseif delay == 0
        return data
    end # if delay

    if n_samples <= delay
        error("Not enough samples ($n_samples) to create embedding with delay ($delay).")
    end # if n_samples

    n_embedded_samples = n_samples - delay
    embedded_dim = n_features * (delay + 1)
    embedded_data = similar(data, embedded_dim, n_embedded_samples)

    for t in 1:n_embedded_samples
        for d in 0:delay
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features
            source_col = t + delay - d
            embedded_data[row_start:row_end, t] = data[:, source_col]
        end # for d
    end # for t
    return embedded_data
end # function delay_embed

function input_lift(U::AbstractMatrix)
    U_lifted = U
    return U_lifted
end # function input_lift

"""
User-provided function to solve for the system matrix M using a specific pseudo-inverse method.
Calculates M = (X_next * X_curr') * pinv(X_curr * X_curr').
The scaling factor 1/size(X_curr,2) in G and H (if used) cancels out in the final expression.
This method is computationally efficient when X_curr (Omega in DMDc) is a "fat" matrix
(i.e., size(X_curr,1) < size(X_curr,2)). If X_curr is "tall",
this method involves inverting a larger matrix (X_curr*X_curr') than
the standard pinv(X_curr) (which would typically use (X_curr'*X_curr) for tall matrices).
"""
function solve_fastdmd(X_curr, X_next)
    # num_snapshots = size(X_curr,2)
    # Original formulation with scaling factors:
    # G_scaled = (1/num_snapshots) * X_next * X_curr'
    # H_scaled = (1/num_snapshots) * X_curr * X_curr'
    # M = G_scaled * pinv(H_scaled)
    # This simplifies to M = (X_next * X_curr') * pinv(X_curr * X_curr')
    # as the scaling factors 1/num_snapshots cancel out.
    G = X_next * X_curr'
    H = X_curr * X_curr'
    return G*pinv(H)
end # function solve_fastdmd

"""
Calculates DMDc matrices (As, Bs, M) from data.
"""
function calculate_dmdc_matrices(X_data::AbstractMatrix, U::AbstractMatrix, delay::Int, N_X::Int, N_U::Int, TOTAL_SAMPLES::Int)
    println("Creating DMDc model with delay $delay...")
    max_delay = delay

    start_idx_curr = max_delay + 2
    end_idx_curr = TOTAL_SAMPLES
    start_idx_hist = 1
    end_idx_hist = TOTAL_SAMPLES - 1

    if start_idx_curr > end_idx_curr || end_idx_hist < start_idx_hist + max_delay
        error("Not enough data points to create DMDc model with delay $max_delay and data size $TOTAL_SAMPLES")
    end # if start_idx_curr

    U = input_lift(U)

    X_current = X_data[:, start_idx_curr:end_idx_curr]
    X_delayed = delay_embed(X_data[:, start_idx_hist:end_idx_hist], max_delay)
    U_delayed = delay_embed(U[:, start_idx_hist:end_idx_hist], max_delay)

    common_len = min(size(X_current, 2), size(X_delayed, 2), size(U_delayed, 2))
    if common_len <= 0
        error("Resulting common length for DMDc matrices is non-positive ($common_len).")
    end # if common_len

    X_next_dmd = X_current[:, 1:common_len]
    X_embed_dmd = X_delayed[:, 1:common_len]
    U_embed_dmd = U_delayed[:, 1:common_len]

    local mean_X
    if NORMALIZE
        mean_X = mean(X_data, dims=2)
    else
        mean_X = zeros(size(X_data, 1), 1)
    end # if NORMALIZE

    X_next_centered = X_next_dmd .- mean_X
    X_embed_centered = X_embed_dmd .- repeat(mean_X, max_delay + 1, 1)

    Omega = [X_embed_centered; U_embed_dmd]

    println("Performing DMDc regression using solve_fastdmd...")
    M = solve_fastdmd(Omega, X_next_centered)
    println("DMDc matrix M size: $(size(M))")

    A_matrices = M[:, 1:N_X*(max_delay+1)]
    B_matrices = M[:, N_X*(max_delay+1)+1:end]

    As = [A_matrices[:, i*N_X+1:(i+1)*N_X] for i in 0:max_delay]
    Bs = [B_matrices[:, i*N_U+1:(i+1)*N_U] for i in 0:max_delay]

    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, mean_X
end # function calculate_dmdc_matrices

"""
Constructs augmented state-space matrices A_full, B_full.
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, delay::Int)
    println("Constructing augmented state-space model (A_full, B_full)...")
    D = delay
    N_STATES_AUG = N_X * (1 + D) + N_U * D

    T = zeros(N_X, N_STATES_AUG)
    for i in 1:(D+1)
        T[:, (i-1)*N_X+1:i*N_X] = As[i]
    end # for i
    if D > 0
        u_del_offset = N_X * (1 + D)
        for i in 1:D
            T[:, u_del_offset+(i-1)*N_U+1:u_del_offset+i*N_U] = Bs[i+1]
        end # for i
    end # if D > 0 (for T matrix u part)

    Sx = zeros(N_X * D, N_X * (1 + D))
    if D > 0
        for i in 0:(D-1)
            Sx[i*N_X+1 : (i+1)*N_X, i*N_X+1 : (i+1)*N_X] = Matrix(I, N_X, N_X) # Or simply I
        end # for i
    end # if D > 0 (for Sx)

    Su = zeros(N_U * D, N_U * D)
    if D > 1 # Correctly zero if D=0 or D=1. Assigns only if D > 1.
        Su[N_U+1 : N_U*D, 1 : N_U*(D-1)] = Matrix(I, N_U*(D-1), N_U*(D-1))
    end # if D > 1 (for Su)

    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    A_full[1:N_X, :] = T
    if D > 0
        A_full[N_X+1:N_X*(1+D), 1:N_X*(1+D)] = Sx
        A_full[N_X*(1+D)+1:end, N_X*(1+D)+1:end] = Su
    end # if D > 0 (for A_full assembly)

    B_full = zeros(N_STATES_AUG, N_U)
    B_full[1:N_X, :] = Bs[1]
    if D > 0
        u_del_start_index_in_aug_state = N_X * (1 + D) + 1
        B_full[u_del_start_index_in_aug_state : u_del_start_index_in_aug_state+N_U-1, :] = Matrix(I, N_U, N_U)
    end # if D > 0 (for B_full assembly)

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end # function construct_augmented_matrices

"""
Main function to create the ABC system from data.
"""
function create_abc_system_from_data(data_path::String, delay::Int)
    X_data, U, N_X, N_U, TOTAL_SAMPLES = load_and_preprocess_data(data_path)
    As, Bs, M, mean_X = calculate_dmdc_matrices(X_data, U, delay, N_X, N_U, TOTAL_SAMPLES)
    A_full, B_full = construct_augmented_matrices(As, Bs, N_X, N_U, delay)
    return A_full, B_full, As, Bs, M, vec(mean_X), N_X, N_U, X_data, U, TOTAL_SAMPLES
end # function create_abc_system_from_data

# --- Simulation and  Reshape Functions ---

"""
Simulates the ABC model (augmented state-space).
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state::AbstractVector, inputs::AbstractMatrix, sim_length::Int)
    n_state_vars = size(A_full, 1)
    _ , n_input_steps = size(inputs)

    if n_input_steps < sim_length
        error("Not enough input data provided for simulation ($sim_length steps needed, $n_input_steps provided).")
    end # if n_input_steps

    states = similar(initial_state, n_state_vars, sim_length)
    current_state = copy(initial_state)
    abc_fn(state, input) = A_full * state + B_full * input

    for k in 1:sim_length
        states[:, k] = current_state
        current_state = abc_fn(current_state, inputs[:, k])
    end # for k
    return states
end # function simulate_abc_model

"""
Simulates the original DMDc model formulation.
"""
function simulate_orig_model(initial_state_orig_model::AbstractVector,
    inputs::AbstractMatrix,
    M_dmdc::AbstractMatrix,
    sim_length::Int,
    delay::Int,
    n_x::Int,
    n_u::Int)

    states_out_centered = zeros(n_x, sim_length)
    current_state_vec = copy(initial_state_orig_model)

    function model_fn(current_embedded_state_and_delayed_inputs, current_input_uk)
        x_embed_centered_from_state_vec = current_embedded_state_and_delayed_inputs[1 : n_x*(delay+1)]
        u_embed_for_M = zeros(n_u * (delay + 1))
        u_embed_for_M[1:n_u] = current_input_uk

        if delay > 0
            u_delayed_part_from_state_vec = current_embedded_state_and_delayed_inputs[n_x*(delay+1)+1 : end]
            u_embed_for_M[n_u+1 : end] = u_delayed_part_from_state_vec
        end # if delay

        omega_k_for_M = [x_embed_centered_from_state_vec; u_embed_for_M]
        return M_dmdc * omega_k_for_M
    end # function model_fn

    for k in 1:sim_length
        current_input = inputs[:, k]
        x_next_centered = model_fn(current_state_vec, current_input)
        states_out_centered[:, k] = x_next_centered

        new_x_embed_centered = zeros(n_x * (delay + 1))
        new_x_embed_centered[1:n_x] = x_next_centered
        if delay > 0
            new_x_embed_centered[n_x+1 : end] = current_state_vec[1 : n_x*delay]
        end # if delay (for x_embed)

        new_u_embed_delayed = zeros(n_u * delay)
        if delay > 0
            new_u_embed_delayed[1:n_u] = current_input
            if delay > 1
                old_u_delayed_part = current_state_vec[n_x*(delay+1)+1 : n_x*(delay+1) + n_u*(delay-1)]
                new_u_embed_delayed[n_u+1 : end] = old_u_delayed_part
            end # if delay > 1
        end # if delay (for u_embed)
        current_state_vec = [new_x_embed_centered; new_u_embed_delayed]
    end # for k
    return states_out_centered
end # function simulate_orig_model

"""
Creates the initial state vector for the ABC model.
"""
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVector, N_X::Int, N_U::Int)
    n_states_aug = N_X * (1 + delay) + N_U * delay
    init_state = zeros(n_states_aug)

    earliest_X_needed = start_time - delay
    earliest_U_needed = start_time - delay

    if earliest_X_needed < 1 || start_time > size(X,2) || (delay > 0 && (earliest_U_needed < 1 || (start_time-1) > size(U,2) || (start_time-1) < 1 ))
         error("Cannot create initial state at time $start_time with delay $delay: requires X data from index $earliest_X_needed to $start_time and U data from $earliest_U_needed to $(start_time-1). Available X=[1, $(size(X,2))], U=[1, $(size(U,2))]")
    end # if boundary check

    init_state[1:N_X] = X[:, start_time] .- mean_x

    for i in 1:delay
        state_idx_start = N_X * i + 1
        state_idx_end = N_X * (i + 1)
        init_state[state_idx_start:state_idx_end] = X[:, start_time-i] .- mean_x
    end # for i (delayed states)

    if delay > 0
        u_del_start_index_in_aug = N_X * (1 + delay) + 1
        for i in 1:delay
            input_idx_start = u_del_start_index_in_aug + N_U * (i - 1)
            input_idx_end = input_idx_start + N_U - 1
            init_state[input_idx_start:input_idx_end] = U[:, start_time-i]
        end # for i (delayed inputs)
    end # if delay > 0 (for U part)

    return init_state
end # function create_abc_init_state

"""
Creates the initial state vector for the original DMDc model simulation.
"""
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVector, N_X::Int, N_U::Int)
    state_len = N_X * (delay + 1) + N_U * delay
    init_state = zeros(state_len)

    earliest_X_needed = start_time - delay
    earliest_U_needed = start_time - delay
    if earliest_X_needed < 1 || start_time > size(X,2) || (delay > 0 && (earliest_U_needed < 1 || (start_time-1) > size(U,2) || (start_time-1) < 1))
        error("Cannot create initial state for original model at time $start_time with delay $delay: data indices out of bounds. Available X=[1, $(size(X,2))], U=[1, $(size(U,2))]")
    end # if boundary check

    for i in 0:delay
        idx_start = i * N_X + 1
        idx_end = (i + 1) * N_X
        init_state[idx_start:idx_end] = X[:, start_time-i] .- mean_x
    end # for i (state part)

    if delay > 0
        u_embed_delayed_start_idx = N_X * (delay + 1) + 1
        for i in 1:delay
            idx_start = u_embed_delayed_start_idx + (i - 1) * N_U
            idx_end = idx_start + N_U - 1
            init_state[idx_start:idx_end] = U[:, start_time-i]
        end # for i (input part)
    end # if delay > 0 (for U part)
    return init_state
end # function create_orig_init_state

"""
Plots the comparison between actual data and model simulations.
"""
function plot_comparison(actual_data::AbstractMatrix, states_abc::AbstractMatrix, states_orig::AbstractMatrix,
    plot_state_idx::Int, sim_start_time::Int, sim_len::Int, delay::Int)

    println("Plotting results...")
    time_axis = 1:sim_len

    p = plot(time_axis, actual_data[plot_state_idx, :],
        label="Actual Data (X$plot_state_idx)",
        linewidth=2, color=:black)
    plot!(p, time_axis, states_abc[plot_state_idx, :],
        label="ABC Model (X$plot_state_idx)",
        linewidth=1.5, linestyle=:dash, color=:red)
    plot!(p, time_axis, states_orig[plot_state_idx, :],
        label="Original DMDc Model (X$plot_state_idx)",
        linewidth=1.5, linestyle=:dot, color=:blue)

    title!("Model Comparison (Delay = $delay)")
    xlabel!("Time Step (relative to $sim_start_time)")
    ylabel!("State Value")
    plot!(p, legend=:outertopright)

    figname = "dmdc_comparison_delay$(delay)_start$(sim_start_time).png"
    try
        savefig(p, figname)
        println("Plot saved to $figname")
    catch e
        println("Error saving plot: $e. Make sure a plotting backend is working (e.g., GR).")
    end # try-catch
end # function plot_comparison


# --- Configuration ---
const DATA_PATH = "./data.jld2"
const DELAY = 100
const START_IDX = 2000
const SIM_LEN = 1000
const PLOT_STATE_IDX = 1
const RELEVANT_STATES = [4]
const NORMALIZE = true

# --- Main Execution ---
function main()
    A_full, B_full, As, Bs, M_dmdc, mean_X_vec, N_X_eff, N_U_eff, X_data_eff, U_eff, TOTAL_SAMPLES_eff = create_abc_system_from_data(DATA_PATH, DELAY)

    sim_start_time = max(START_IDX, DELAY + 1 + 1) # ensure DELAY steps before sim_start_time for U, and X data for X_current
                                                 # X needs up to start_time, U needs up to start_time-1
                                                 # initial state needs X[start_time-DELAY], U[start_time-DELAY]
                                                 # So start_time - DELAY >= 1 -> start_time >= DELAY + 1
    if sim_start_time > TOTAL_SAMPLES_eff
        println("Error: Simulation start time $sim_start_time is beyond total samples $TOTAL_SAMPLES_eff.")
        return
    end # if sim_start_time

    actual_sim_len = SIM_LEN
    if sim_start_time + actual_sim_len - 1 > TOTAL_SAMPLES_eff
        println("Warning: Simulation length exceeds data bounds for X_data. Adjusting sim_len.")
        actual_sim_len = TOTAL_SAMPLES_eff - sim_start_time + 1
    end # if X_data bounds
    if sim_start_time + actual_sim_len - 1 > size(U_eff,2)
         println("Warning: Simulation length exceeds data bounds for U_data. Adjusting sim_len further.")
         actual_sim_len = size(U_eff,2) - sim_start_time + 1
    end # if U_data bounds

    if actual_sim_len <= 0
        println("Error: Calculated simulation length is $actual_sim_len. Cannot simulate.")
        return
    end # if actual_sim_len

    sim_inputs = U_eff[:, sim_start_time : sim_start_time+actual_sim_len-1]

    if sim_start_time - DELAY < 1
        println("Error: sim_start_time ($sim_start_time) is too early for the given DELAY ($DELAY). Cannot create initial conditions.")
        return
    end # if initial condition check

    init_state_abc = create_abc_init_state(X_data_eff, U_eff, sim_start_time, DELAY, mean_X_vec, N_X_eff, N_U_eff)
    init_state_orig = create_orig_init_state(X_data_eff, U_eff, sim_start_time, DELAY, mean_X_vec, N_X_eff, N_U_eff)

    println("Running simulations for $actual_sim_len steps...")
    states_abc_full = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, actual_sim_len)
    states_abc = states_abc_full[1:N_X_eff, :] .+ mean_X_vec

    states_orig_centered = simulate_orig_model(init_state_orig, sim_inputs, M_dmdc, actual_sim_len, DELAY, N_X_eff, N_U_eff)
    states_orig = states_orig_centered .+ mean_X_vec

    actual_data = X_data_eff[:, sim_start_time : sim_start_time+actual_sim_len-1]
    println("Simulations complete.")

    plot_idx_to_use = PLOT_STATE_IDX
    if PLOT_STATE_IDX < 1 || PLOT_STATE_IDX > N_X_eff
        println("Warning: PLOT_STATE_IDX ($PLOT_STATE_IDX) is out of bounds for N_X=$N_X_eff. Plotting first state.")
        plot_idx_to_use = 1
    end # if PLOT_STATE_IDX
    plot_comparison(actual_data, states_abc, states_orig, plot_idx_to_use, sim_start_time, actual_sim_len, DELAY)

    println("Script finished.")
    jldsave("ABC.jld2";
        A=A_full,
        B=B_full,
        delay=DELAY,
        mean_vec=mean_X_vec,
        M_dmdc=M_dmdc,
        N_X_effective=N_X_eff,
        N_U_effective=N_U_eff,
        TOTAL_SAMPLES_X=TOTAL_SAMPLES_eff,
        RELEVANT_STATES_indices=RELEVANT_STATES
    )
    println("Results saved to ABC.jld2")

end # function main

# To run the script:
# main()
