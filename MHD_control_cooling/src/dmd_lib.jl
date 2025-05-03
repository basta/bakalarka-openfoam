#=
This script performs Dynamic Mode Decomposition with Control (DMDc)
on time-series data (X_data, U_data) with time delays.
It constructs an augmented state-space model (A_full, B_full)
and compares its simulation results against the original DMDc formulation
and the actual data.

Key Functions:
- `load_and_preprocess_data`: Loads data from JLD2, preprocesses control input U.
- `delay_embed`: Creates delay-embedded matrices.
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
        X_data = X_data[RELEVANT_STATES, 2:end]

        U_data = jldopen(data_path, "r")["U_data"]

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
        end

        N_X = size(X_data, 1)
        N_U = size(U, 1)
        TOTAL_SAMPLES = size(X_data, 2)
        println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")

        return X_data, U, N_X, N_U, TOTAL_SAMPLES
    catch e
        println("Error loading data from $data_path: $e")
        println("Please ensure the JLD2 file exists and the path is correct.")
        rethrow(e) # Propagate the error
    end
end

"""
Creates a delay-embedded matrix from time series data.
"""
function delay_embed(data::AbstractMatrix, delay::Int)
    n_features, n_samples = size(data)
    if delay < 0
        error("Delay must be non-negative.")
    elseif delay == 0
        return data
    end

    if n_samples <= delay
        error("Not enough samples ($n_samples) to create embedding with delay ($delay).")
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

function input_lift(U::AbstractMatrix)
    U_lifted = U
    return U_lifted
end

"""
Calculates DMDc matrices (As, Bs, M) from data.
"""
function calculate_dmdc_matrices(X_data::AbstractMatrix, U::AbstractMatrix, delay::Int, N_X::Int, N_U::Int, TOTAL_SAMPLES::Int)
    println("Creating DMDc model with delay $delay...")
    max_delay = delay

    # Prepare data indices
    start_idx_curr = max_delay + 2
    end_idx_curr = TOTAL_SAMPLES
    start_idx_hist = 1
    end_idx_hist = TOTAL_SAMPLES - 1

    if start_idx_curr > end_idx_curr || end_idx_hist < start_idx_hist + max_delay
        error("Not enough data points to create DMDc model with delay $max_delay and data size $TOTAL_SAMPLES")
    end

    U = input_lift(U)

    X_current = X_data[:, start_idx_curr:end_idx_curr]
    X_delayed = delay_embed(X_data[:, start_idx_hist:end_idx_hist], max_delay)
    U_delayed = delay_embed(U[:, start_idx_hist:end_idx_hist], max_delay)

    # Align dimensions
    common_len = min(size(X_current, 2), size(X_delayed, 2), size(U_delayed, 2))
    if common_len <= 0
        error("Resulting common length for DMDc matrices is non-positive ($common_len).")
    end

    X_next_dmd = X_current[:, 1:common_len]
    X_embed_dmd = X_delayed[:, 1:common_len]
    U_embed_dmd = U_delayed[:, 1:common_len]
    #==== LIFITING HERE ====#
    U_embed_dmd = abs.(U_embed_dmd)

    # Center X data
    local mean_X
    if NORMALIZE
        mean_X = mean(X_data, dims=2)
    else
        mean_X = zeros(size(X_data, 1), 1)
    end
    X_next_centered = X_next_dmd .- mean_X
    X_embed_centered = X_embed_dmd .- repeat(mean_X, max_delay + 1, 1)

    # Combined state for regression
    Omega = [X_embed_centered; U_embed_dmd]

    # Solve for the combined dynamics matrix M
    println("Performing DMDc regression...")
    M = X_next_centered * pinv(Omega)
    println("DMDc matrix M size: $(size(M))")

    # Extract A and B matrices
    A_matrices = M[:, 1:N_X*(max_delay+1)]
    B_matrices = M[:, N_X*(max_delay+1)+1:end]

    # Reshape into list of As and Bs
    As = [A_matrices[:, i*N_X+1:(i+1)*N_X] for i in 0:max_delay]
    Bs = [B_matrices[:, i*N_U+1:(i+1)*N_U] for i in 0:max_delay]

    println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")
    return As, Bs, M, mean_X
end

"""
Constructs augmented state-space matrices A_full, B_full.
"""
function construct_augmented_matrices(As::Vector, Bs::Vector, N_X::Int, N_U::Int, delay::Int)
    println("Constructing augmented state-space model (A_full, B_full)...")
    D = delay
    N_STATES_AUG = N_X * (1 + D) + N_U * D

    # Construct T matrix (top row of A_full)
    T = zeros(N_X, N_STATES_AUG)
    for i in 1:(D+1)
        T[:, (i-1)*N_X+1:i*N_X] = As[i]
    end
    if D > 0
        u_del_offset = N_X * (1 + D)
        for i in 1:D
            T[:, u_del_offset+(i-1)*N_U+1:u_del_offset+i*N_U] = Bs[i+1]
        end
    end

    # Construct Sx matrix (state shifting)
    Sx = zeros(N_X * D, N_X * (1 + D))
    if D > 0
        Sx[1:N_X, 1:N_X] = I(N_X)
        if D > 1
            Sx[N_X+1:N_X*D, N_X+1:N_X*D] = I(N_X * (D - 1))
        end
    end

    # Construct Su matrix (input shifting)
    Su = zeros(N_U * D, N_U * D)
    if D > 1
        Su[N_U+1:N_U*D, 1:N_U*(D-1)] = I(N_U * (D - 1))
    end

    # Assemble A_full
    A_full = zeros(N_STATES_AUG, N_STATES_AUG)
    A_full[1:N_X, :] = T
    if D > 0
        A_full[N_X+1:N_X*(1+D), 1:N_X*(1+D)] = Sx
        A_full[N_X*(1+D)+1:end, N_X*(1+D)+1:end] = Su
    end

    # Assemble B_full
    B_full = zeros(N_STATES_AUG, N_U)
    B_full[1:N_X, :] = Bs[1]
    if D > 0
        u_del_start_index = N_X * (1 + D) + 1
        B_full[u_del_start_index:u_del_start_index+N_U-1, :] = I(N_U)
    end

    println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")
    return A_full, B_full
end

"""
Main function to create the ABC system from data.
"""
function create_abc_system_from_data(data_path::String, delay::Int)
    X_data, U, N_X, N_U, TOTAL_SAMPLES = load_and_preprocess_data(data_path)
    As, Bs, M, mean_X = calculate_dmdc_matrices(X_data, U, delay, N_X, N_U, TOTAL_SAMPLES)
    A_full, B_full = construct_augmented_matrices(As, Bs, N_X, N_U, delay)

    # Return all relevant components
    return A_full, B_full, As, Bs, M, vec(mean_X), N_X, N_U, X_data, U, TOTAL_SAMPLES
end

# --- Simulation and  Reshape Functions ---

"""
Simulates the ABC model (augmented state-space).
"""
function simulate_abc_model(A_full::AbstractMatrix, B_full::AbstractMatrix, initial_state::AbstractVector, inputs::AbstractMatrix, sim_length::Int)
    n_state_vars = size(A_full, 1)
    n_inputs, n_input_steps = size(inputs)

    if n_input_steps < sim_length
        error("Not enough input data provided for simulation ($sim_length steps needed, $n_input_steps provided).")
    end

    states = similar(initial_state, n_state_vars, sim_length)
    current_state = copy(initial_state)
    abc_fn(state, input) = A_full * state + B_full * input

    for k in 1:sim_length
        states[:, k] = current_state
        current_state = abc_fn(current_state, inputs[:, k])
    end
    return states
end

"""
Simulates the original DMDc model formulation.
Requires careful state updates within the loop.
"""
function simulate_orig_model(initial_state_orig_model::AbstractVector,
    inputs::AbstractMatrix,
    M_dmdc::AbstractMatrix, # Pass the DMDc matrix M
    sim_length::Int,
    delay::Int,
    n_x::Int,
    n_u::Int)

    states_out_centered = zeros(n_x, sim_length)
    current_state_vec = copy(initial_state_orig_model) # [x(k)..x(k-D) centered; u(k-1)..u(k-D)]

    # Define the model function internally using M_dmdc
    function model_fn(current_embedded_state_u, current_input)
        x_embed_centered = current_embedded_state_u[1:n_x*(delay+1)]
        u_embed = zeros(n_u * (delay + 1))
        u_embed[1:n_u] = current_input # u(k)
        if delay > 0
            u_delayed_part_start_idx = n_x * (delay + 1) + 1
            u_delayed_part_end_idx = u_delayed_part_start_idx + n_u * delay - 1
            if length(current_embedded_state_u) >= u_delayed_part_end_idx
                u_delayed_part = current_embedded_state_u[u_delayed_part_start_idx:u_delayed_part_end_idx]
                u_embed[n_u+1:end] = u_delayed_part
            else
                println("Warning: Index out of bounds accessing delayed U in simulation.")
            end
        end
        omega_k = [x_embed_centered; u_embed]
        return M_dmdc * omega_k
    end

    for k in 1:sim_length
        current_input = inputs[:, k] # u(k)
        x_next_centered = model_fn(current_state_vec, current_input)
        states_out_centered[:, k] = x_next_centered

        # Update state vector for the next iteration
        new_x_embed_centered = zeros(n_x * (delay + 1))
        new_x_embed_centered[1:n_x] = x_next_centered
        if delay > 0
            new_x_embed_centered[n_x+1:end] = current_state_vec[1:n_x*delay]
        end

        new_u_embed_delayed = zeros(n_u * delay)
        if delay > 0
            new_u_embed_delayed[1:n_u] = current_input # u(k) becomes u(k-1) next step
            if delay > 1
                u_delayed_part_start_idx = n_x * (delay + 1) + 1
                old_u_delayed_len = n_u * (delay - 1)
                old_u_delayed_end_idx = u_delayed_part_start_idx + old_u_delayed_len - 1
                # Ensure indices are valid before accessing
                if length(current_state_vec) >= old_u_delayed_end_idx
                    old_u_delayed = current_state_vec[u_delayed_part_start_idx:old_u_delayed_end_idx] # u(k-1)...u(k-D+1)
                    new_u_embed_delayed[n_u+1:end] = old_u_delayed
                else
                    println("Warning: Index out of bounds updating delayed U in simulation.")
                end
            end
        end
        current_state_vec = [new_x_embed_centered; new_u_embed_delayed]
    end
    return states_out_centered
end


"""
Creates the initial state vector for the ABC model.
State: [x(k); x(k-1); ...; x(k-D); u(k-1); ...; u(k-D)]
"""
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVector, N_X::Int, N_U::Int)
    n_states_aug = N_X * (1 + delay) + N_U * delay
    init_state = zeros(n_states_aug)

    # Check boundaries
    max_required_time = start_time
    min_required_time = start_time - delay
    if min_required_time < 1 || max_required_time > size(X, 2) || max_required_time > size(U, 2)
        error("Cannot create initial state at time $start_time with delay $delay: requires data indices [$min_required_time, $max_required_time], available X=[1, $(size(X,2))], U=[1, $(size(U,2))]")
    end

    # Current state x(k) centered
    init_state[1:N_X] = X[:, start_time] .- mean_x

    # Delayed states x(k-1) ... x(k-D) centered
    for i in 1:delay
        state_idx_start = N_X * i + 1
        state_idx_end = N_X * (i + 1)
        init_state[state_idx_start:state_idx_end] = X[:, start_time-i] .- mean_x
    end

    # Delayed inputs u(k-1) ... u(k-D)
    if delay > 0
        u_del_start_index = N_X * (1 + delay) + 1
        for i in 1:delay
            input_idx_start = u_del_start_index + N_U * (i - 1)
            input_idx_end = input_idx_start + N_U - 1
            init_state[input_idx_start:input_idx_end] = U[:, start_time-i]
        end
    end

    return init_state
end

"""
Creates the initial state vector for the original DMDc model simulation.
State: [x(k); ... x(k-D) centered; u(k-1); ... u(k-D)]
"""
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVector, N_X::Int, N_U::Int)
    state_len = N_X * (delay + 1) + N_U * delay

    # Check boundaries
    max_required_time = start_time
    min_required_time = start_time - delay
    if min_required_time < 1 || max_required_time > size(X, 2) || max_required_time > size(U, 2)
        error("Cannot create initial state at time $start_time with delay $delay: requires data indices [$min_required_time, $max_required_time], available X=[1, $(size(X,2))], U=[1, $(size(U,2))]")
    end

    init_state = zeros(state_len)

    # State part: [x(k); x(k-1); ...; x(k-D)] centered
    for i in 0:delay
        idx_start = i * N_X + 1
        idx_end = (i + 1) * N_X
        init_state[idx_start:idx_end] = X[:, start_time-i] .- mean_x
    end

    # Input part: [u(k-1); ...; u(k-D)]
    if delay > 0
        u_embed_delayed_start_idx = N_X * (delay + 1) + 1
        for i in 1:delay
            idx_start = u_embed_delayed_start_idx + (i - 1) * N_U
            idx_end = idx_start + N_U - 1
            init_state[idx_start:idx_end] = U[:, start_time-i]
        end
    end
    return init_state
end

"""
Plots the comparison between actual data and model simulations.
"""
function plot_comparison(actual_data::AbstractMatrix, states_abc::AbstractMatrix, states_orig::AbstractMatrix,
    plot_state_idx::Int, sim_start_time::Int, sim_len::Int, delay::Int)

    println("Plotting results...")
    plot(1:sim_len, actual_data[plot_state_idx, :],
        label="Actual Data (X$plot_state_idx)",
        linewidth=2, color=:black)
    plot!(1:sim_len, states_abc[plot_state_idx, :],
        label="ABC Model (X$plot_state_idx)",
        linewidth=1.5, linestyle=:dash, color=:red)
    plot!(1:sim_len, states_orig[plot_state_idx, :],
        label="Original DMDc Model (X$plot_state_idx)",
        linewidth=1.5, linestyle=:dot, color=:blue)

    title!("Model Comparison (Delay = $delay)")
    xlabel!("Time Step (relative to $sim_start_time)")
    ylabel!("State Value")
    plot!(legend=:outertopright)

    # Save the plot
    figname = "dmdc_comparison_delay$(delay)_start$(sim_start_time).png"
    savefig(figname)
    println("Plot saved to $figname")
end


# --- Configuration ---
const DATA_PATH = "./data.jld2" # Adjust path if needed
const DELAY = 30       # Number of delay steps for states and inputs
const START_IDX = 800   # Starting index for simulation comparison
const SIM_LEN = 700     # Default length of the simulation
const PLOT_STATE_IDX = 1 # Which state component to plot
const RELEVANT_STATES = [6] # Relevant states for analysis
const NORMALIZE = false

# --- Main Execution ---
function main()
    # 1. Create ABC System
    A_full, B_full, As, Bs, M_dmdc, mean_X_vec, N_X, N_U, X_data, U, TOTAL_SAMPLES = create_abc_system_from_data(DATA_PATH, DELAY)

    # 2. Simulation Setup
    sim_start_time = max(START_IDX, DELAY + 1)
    sim_len = SIM_LEN # Use default or modify if needed
    if sim_start_time + sim_len - 1 > TOTAL_SAMPLES
        println("Warning: Simulation length exceeds data bounds. Adjusting sim_len.")
        sim_len = TOTAL_SAMPLES - sim_start_time + 1
    end
    sim_inputs = U[:, sim_start_time:sim_start_time+sim_len-1]

    # 3. Create Initial States
    init_state_abc = create_abc_init_state(X_data, U, sim_start_time, DELAY, mean_X_vec, N_X, N_U)
    init_state_orig = create_orig_init_state(X_data, U, sim_start_time, DELAY, mean_X_vec, N_X, N_U)

    # 4. Run Simulations
    println("Running simulations...")
    states_abc_full = simulate_abc_model(A_full, B_full, init_state_abc, sim_inputs, sim_len)
    states_abc = states_abc_full[1:N_X, :] .+ mean_X_vec # Add mean back

    states_orig_centered = simulate_orig_model(init_state_orig, sim_inputs, M_dmdc, sim_len, DELAY, N_X, N_U)
    states_orig = states_orig_centered .+ mean_X_vec # Add mean back

    # 5. Get Actual Data for Comparison
    actual_data = X_data[:, sim_start_time:sim_start_time+sim_len-1]
    println("Simulations complete.")

    # 6. Plot Results
    plot_comparison(actual_data, states_abc, states_orig, PLOT_STATE_IDX, sim_start_time, sim_len, DELAY)

    println("Script finished.")
    jldsave("ABC.jld2", A=A_full, B=B_full, delay=DELAY, mean_vec=mean_X_vec)
end
