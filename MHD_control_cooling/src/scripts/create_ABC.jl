# --- Julia DMDc Delay Comparison Script ---

using JLD2
using LinearAlgebra
using Statistics
using Plots # Using Plots.jl for simplicity. Add Pkg.add("Plots") if needed.

# --- Configuration ---
const DATA_PATH = "./data/dataset-live-fast2.jld2" # Adjust path if needed
const DELAY = 300      # Number of delay steps for states and inputs
const START_IDX = 1200 # Starting index for simulation comparison
const SIM_LEN = 700     # Length of the simulation
const PLOT_STATE_IDX = 1 # Which state component to plot

# --- Data Loading and Preprocessing ---
println("Loading data...")
try
    global X_data = jldopen(DATA_PATH, "r")["dataset"][2]
    global U_data = jldopen(DATA_PATH, "r")["dataset"][1]
catch e
    println("Error loading data from $DATA_PATH: $e")
    println("Please ensure the JLD2 file exists and the path is correct.")
    exit()
end

X_data = X_data[end:end, :]
@info "Loaded X_data with size $(size(X_data))"

# Preprocess U_data (outer product)
println("Preprocessing U data...")
U = let
    cols = []
    for u_col in eachcol(U_data)
        # Assuming u_col has 8 elements: first 4 for u_e, next 4 for u_m
        if size(u_col, 1) == 8
            push!(cols, vec(u_col[1:4] * u_col[5:8]'))
        else
            println("Warning: Unexpected size for column in U_data. Skipping.")
            # Handle error or different structure if necessary
        end
    end
    hcat(cols...) # Use hcat for efficiency
end

const N_X = size(X_data, 1)
const N_U = size(U, 1)
const TOTAL_SAMPLES = size(X_data, 2)

println("Data loaded: X_data ($(size(X_data))), U ($(size(U)))")
println("Parameters: DELAY=$DELAY, N_X=$N_X, N_U=$N_U")

# --- Helper Functions ---

"""
Creates a delay-embedded matrix from time series data.
Input: data (features x time), delay (number of steps)
Output: embedded_data (features*(delay+1) x time-(delay))
"""
function delay_embed(data::AbstractMatrix, delay::Int)
    n_features, n_samples = size(data)
    if delay == 0
        return data
    end
    n_embedded_samples = n_samples - delay
    embedded_dim = n_features * (delay + 1)
    embedded_data = similar(data, embedded_dim, n_embedded_samples)

    for t in 1:n_embedded_samples
        for d in 0:delay
            row_start = d * n_features + 1
            row_end = (d + 1) * n_features
            # Ensure indices are valid
            source_col = t + delay - d
            if source_col >= 1 && source_col <= n_samples
                embedded_data[row_start:row_end, t] = data[:, source_col]
            else
                # Handle out-of-bounds access if necessary, e.g., fill with zeros
                embedded_data[row_start:row_end, t] .= 0.0
                println("Warning: Out-of-bounds access in delay_embed at t=$t, d=$d")
            end
        end
    end
    return embedded_data
end

"""
Simulates a discrete-time system: x_{k+1} = f(x_k, u_k)
"""
function simulate_model(initial_state::AbstractVector,
    inputs::AbstractMatrix,
    model_fn::Function,
    sim_length::Int)

    n_state_vars = length(initial_state)
    n_inputs, n_input_steps = size(inputs)

    if n_input_steps < sim_length
        error("Not enough input data provided for the simulation length ($sim_length), need $n_input_steps.")
    end

    # Preallocate states matrix based on the expected output type of model_fn
    # We assume model_fn returns a vector of the same type as initial_state
    states = similar(initial_state, n_state_vars, sim_length)
    current_state = copy(initial_state)

    for k in 1:sim_length
        current_input = inputs[:, k]
        # Store the *current* state before updating
        states[:, k] = current_state
        # Update the state for the *next* iteration
        current_state = model_fn(current_state, current_input)
    end
    return states
end


# --- Model Creation: DMDc with Delays (Simplified) ---
println("Creating DMDc model with delays...")

# Prepare data for DMDc
# We need x(k+1) predicted from x(k)...x(k-D) and u(k)...u(k-D)
max_delay = DELAY
# Ensure indices are valid
start_idx_curr = max_delay + 2
end_idx_curr = TOTAL_SAMPLES

start_idx_hist = 1
end_idx_hist = TOTAL_SAMPLES - 1

if start_idx_curr > end_idx_curr
    error("Not enough data points to create DMDc model with delay $max_delay")
end

X_current = X_data[:, start_idx_curr:end_idx_curr]
X_delayed = delay_embed(X_data[:, start_idx_hist:end_idx_hist], max_delay) # Embeds x(k)...x(k-D)
U_delayed = delay_embed(U[:, start_idx_hist:end_idx_hist], max_delay)       # Embeds u(k)...u(k-D)

# Align dimensions
common_len = min(size(X_current, 2), size(X_delayed, 2), size(U_delayed, 2))
if common_len <= 0
    error("Resulting common length for DMDc matrices is non-positive ($common_len). Check data size and delay.")
end

X_next_dmd = X_current[:, 1:common_len]
X_embed_dmd = X_delayed[:, 1:common_len]
U_embed_dmd = U_delayed[:, 1:common_len]

# Center X data for DMDc
mean_X = mean(X_data, dims=2)
X_next_centered = X_next_dmd .- mean_X
X_embed_centered = X_embed_dmd .- repeat(mean_X, max_delay + 1, 1) # Center embedded states

# Combined state for regression: [x_embed_centered; u_embed]
Omega = [X_embed_centered; U_embed_dmd]

# Solve for the combined dynamics matrix M: X_next_centered = M * Omega
# M contains [A0 A1 ... AD B0 B1 ... BD]
println("Performing DMDc regression...")
M = X_next_centered * pinv(Omega)
println("DMDc matrix M size: $(size(M))")

# Extract A and B matrices
# M = [A_matrices B_matrices]
A_matrices = M[:, 1:N_X*(max_delay+1)]
B_matrices = M[:, N_X*(max_delay+1)+1:end]

# Reshape into list of As and Bs
As = [A_matrices[:, i*N_X+1:(i+1)*N_X] for i in 0:max_delay]
Bs = [B_matrices[:, i*N_U+1:(i+1)*N_U] for i in 0:max_delay]

println("Extracted $(length(As)) A matrices and $(length(Bs)) B matrices.")

# --- Define Original Style Model Function ---
function model_fn_orig(current_embedded_state_u::AbstractVector, current_input::AbstractVector)
    # current_embedded_state_u should be [x(k); ... x(k-D); u(k-1); ... u(k-D)] centered
    # We need to form the Omega vector for prediction: [x(k)...x(k-D); u(k)...u(k-D)]
    # Note: DMDc was trained with u(k)...u(k-D), but simulation step uses u(k) explicitly

    # Extract centered embedded state part
    x_embed_centered = current_embedded_state_u[1:N_X*(DELAY+1)]

    # Construct U_embed part for prediction [u(k); u(k-1); ...; u(k-D)]
    u_embed = zeros(N_U * (DELAY + 1))
    u_embed[1:N_U] = current_input # u(k)
    if DELAY > 0
        # Get u(k-1)...u(k-D) from the input state vector
        # Ensure indices are correct based on state vector structure
        u_delayed_part_start_idx = N_X * (DELAY + 1) + 1
        u_delayed_part_end_idx = u_delayed_part_start_idx + N_U * DELAY - 1
        if length(current_embedded_state_u) >= u_delayed_part_end_idx
            u_delayed_part = current_embedded_state_u[u_delayed_part_start_idx:u_delayed_part_end_idx]
            u_embed[N_U+1:end] = u_delayed_part
        else
            println("Warning: Index out of bounds when accessing delayed U in model_fn_orig")
            # Handle error or default values if needed
        end
    end

    omega_k = [x_embed_centered; u_embed]
    x_next_centered = M * omega_k
    return x_next_centered # Return the centered next state
end

# --- Define Initial State Function for Original Model ---
function create_orig_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVector)
    n_x = size(X, 1)
    n_u = size(U, 1)
    state_len = n_x * (delay + 1) + n_u * delay # Correct length

    # Check boundaries
    if start_time - delay < 1
        error("start_time ($start_time) is too early for the given delay ($delay)")
    end
    max_required_time = start_time
    min_required_time = start_time - delay
    if max_required_time > size(X, 2) || max_required_time > size(U, 2) || min_required_time < 1
        error("start_time ($start_time) with delay ($delay) requires data indices outside available range [1, $(size(X,2))]")
    end

    init_state = zeros(state_len)

    # State part: [x(k); x(k-1); ...; x(k-D)] centered
    for i in 0:delay
        idx_start = i * n_x + 1
        idx_end = (i + 1) * n_x
        init_state[idx_start:idx_end] = X[:, start_time-i] .- mean_x
    end

    # Input part: [u(k-1); ...; u(k-D)]
    if delay > 0
        u_embed_delayed_start_idx = n_x * (delay + 1) + 1
        for i in 1:delay
            idx_start = u_embed_delayed_start_idx + (i - 1) * n_u
            idx_end = idx_start + n_u - 1
            init_state[idx_start:idx_end] = U[:, start_time-i]
        end
    end

    return init_state
end


# --- Model Creation: Augmented State-Space (ABC Model) ---
println("Constructing augmented state-space model (A_full, B_full)...")

const D = DELAY # Alias for clarity
const N_STATES_AUG = N_X * (1 + D) + N_U * D # Size of the augmented state

# Construct T matrix (top row of A_full)
T = zeros(N_X, N_STATES_AUG)
# State part: As[1]*x(k) + As[2]*x(k-1) + ...
for i in 1:(D+1)
    T[:, (i-1)*N_X+1:i*N_X] = As[i] # As[1] corresponds to x(k), As[2] to x(k-1) etc.
end
# Delayed input part: Bs[2]*u(k-1) + Bs[3]*u(k-2) + ...
if D > 0
    u_del_offset = N_X * (1 + D)
    for i in 1:D
        T[:, u_del_offset+(i-1)*N_U+1:u_del_offset+i*N_U] = Bs[i+1] # Bs[2] corresponds to u(k-1)
    end
end

# Construct Sx matrix (state shifting)
Sx = zeros(N_X * D, N_X * (1 + D))
if D > 0
    # Copy x(k) to the first block of the next delayed state [x(k)]_{k+1}
    Sx[1:N_X, 1:N_X] = I(N_X)
    # Copy x(k-1)...x(k-D+1) to the remaining blocks [x(k-1);...;x(k-D+1)]_{k+1}
    if D > 1
        Sx[N_X+1:N_X*D, N_X+1:N_X*D] = I(N_X * (D - 1)) # Selects x(k-1)...x(k-D+1) from previous state
    end
end

# Construct Su matrix (input shifting)
Su = zeros(N_U * D, N_U * D)
if D > 1
    # Copy u(k-1)...u(k-D+1) for the next state's delayed inputs
    Su[N_U+1:N_U*D, 1:N_U*(D-1)] = I(N_U * (D - 1)) # Selects u(k-1)...u(k-D+1) from previous u_del state
end

# Assemble A_full
A_full = zeros(N_STATES_AUG, N_STATES_AUG)
A_full[1:N_X, :] = T # Dynamics for x(k+1)
if D > 0
    # State shifting part
    A_full[N_X+1:N_X*(1+D), 1:N_X*(1+D)] = Sx
    # Input shifting part
    A_full[N_X*(1+D)+1:end, N_X*(1+D)+1:end] = Su
end

# Assemble B_full
B_full = zeros(N_STATES_AUG, N_U)
B_full[1:N_X, :] = Bs[1] # Effect of u(k) on x(k+1)
if D > 0
    # Place u(k) into the start of the next u_del block [u(k); u(k-1); ...]_{k+1}
    u_del_start_index = N_X * (1 + D) + 1
    B_full[u_del_start_index:u_del_start_index+N_U-1, :] = I(N_U)
end

println("A_full size: $(size(A_full)), B_full size: $(size(B_full))")

# --- Define ABC Model Function ---
abc_fn(state::AbstractVector, input::AbstractVector) = A_full * state + B_full * input

# --- Define Initial State Function for ABC Model ---
function create_abc_init_state(X::AbstractMatrix, U::AbstractMatrix, start_time::Int, delay::Int, mean_x::AbstractVector)
    n_x = size(X, 1)
    n_u = size(U, 1)
    n_states_aug = n_x * (1 + delay) + n_u * delay
    init_state = zeros(n_states_aug)

    # Check boundaries
    if start_time - delay < 1
        error("start_time ($start_time) is too early for the given delay ($delay)")
    end
    max_required_time = start_time
    min_required_time = start_time - delay
    if max_required_time > size(X, 2) || max_required_time > size(U, 2) || min_required_time < 1
        error("start_time ($start_time) with delay ($delay) requires data indices outside available range [1, $(size(X,2))]")
    end

    # Current state x(k) centered
    init_state[1:n_x] = X[:, start_time] .- mean_x

    # Delayed states x(k-1) ... x(k-D) centered
    for i in 1:delay
        state_idx_start = n_x * i + 1
        state_idx_end = n_x * (i + 1)
        init_state[state_idx_start:state_idx_end] = X[:, start_time-i] .- mean_x
    end

    # Delayed inputs u(k-1) ... u(k-D)
    u_del_start_index = n_x * (1 + D) + 1
    for i in 1:delay
        input_idx_start = u_del_start_index + n_u * (i - 1)
        input_idx_end = input_idx_start + n_u - 1
        init_state[input_idx_start:input_idx_end] = U[:, start_time-i]
    end

    return init_state
end

# --- Simulation ---
println("Running simulations...")

# Ensure start index allows for delay history
sim_start_time = max(START_IDX, DELAY + 1)
global SIM_LEN # Allow modification
if sim_start_time + SIM_LEN - 1 > TOTAL_SAMPLES
    println("Warning: Simulation length exceeds data bounds. Adjusting SIM_LEN.")
    SIM_LEN = TOTAL_SAMPLES - sim_start_time + 1
end

# Inputs for simulation
sim_inputs = U[:, sim_start_time:sim_start_time+SIM_LEN-1]

# Initial state for ABC model
# Convert mean_X matrix to vector for the function call
init_state_abc = create_abc_init_state(X_data, U, sim_start_time, DELAY, vec(mean_X))

# Initial state for Original model
# Convert mean_X matrix to vector for the function call
init_state_orig = create_orig_init_state(X_data, U, sim_start_time, DELAY, vec(mean_X))


# Simulate ABC model
# The state vector includes delayed states/inputs, output is the next state vector
states_abc_full = simulate_model(init_state_abc, sim_inputs, abc_fn, SIM_LEN)
# Extract the actual state x(k) (which is centered) and add mean back
states_abc = states_abc_full[1:N_X, :] .+ mean_X

# Simulate Original model with correct state propagation
function simulate_model_orig(initial_state_orig_model::AbstractVector,
    inputs::AbstractMatrix,
    model_fn::Function,
    sim_length::Int,
    mean_x_vec::AbstractVector, # Expect vector here
    delay::Int,
    n_x::Int,
    n_u::Int)

    states_out_centered = zeros(n_x, sim_length)
    # State vector for original model: [x(k)..x(k-D) centered; u(k-1)..u(k-D)]
    current_state_vec = copy(initial_state_orig_model)

    for k in 1:sim_length
        current_input = inputs[:, k] # u(k)

        # Predict next centered state x(k+1) using the specific model function
        x_next_centered = model_fn(current_state_vec, current_input)
        states_out_centered[:, k] = x_next_centered

        # --- Update state vector for the *next* iteration ---
        # The state vector needs to be [x(k+1)..x(k+1-D) centered; u(k)..u(k+1-D)]

        # Shift x states: [x(k+1); x(k); ...; x(k+1-D)] all centered
        new_x_embed_centered = zeros(n_x * (delay + 1))
        new_x_embed_centered[1:n_x] = x_next_centered # x(k+1) centered
        if delay > 0
            # Copy x(k)...x(k+1-D) from previous state vector's x(k)...x(k-D+1) part
            new_x_embed_centered[n_x+1:end] = current_state_vec[1:n_x*delay]
        end

        # Shift u states: [u(k); u(k-1); ...; u(k+1-D)]
        new_u_embed_delayed = zeros(n_u * delay)
        if delay > 0
            new_u_embed_delayed[1:n_u] = current_input # u(k)
            if delay > 1
                # u(k-1)...u(k+1-D) come from old u(k-1)...u(k-D+1)
                u_delayed_part_start_idx = n_x * (delay + 1) + 1
                old_u_delayed = current_state_vec[u_delayed_part_start_idx:(u_delayed_part_start_idx+n_u*(delay-1)-1)] # u(k-1)...u(k-D+1)
                new_u_embed_delayed[n_u+1:end] = old_u_delayed
            end
        end
        current_state_vec = [new_x_embed_centered; new_u_embed_delayed]
        # --- End state update ---
    end
    return states_out_centered
end

# Call the redefined simulation for the original model
# Convert mean_X matrix to vector for the function call
states_orig_centered = simulate_model_orig(init_state_orig, sim_inputs, model_fn_orig, SIM_LEN, vec(mean_X), DELAY, N_X, N_U)
states_orig = states_orig_centered .+ mean_X # Add mean back

# Actual data for comparison
actual_data = X_data[:, sim_start_time:sim_start_time+SIM_LEN-1]

println("Simulations complete.")

# --- Plotting ---
println("Plotting results...")

time_vec = (sim_start_time:sim_start_time+SIM_LEN-1) # Or simply 1:SIM_LEN

plot(1:SIM_LEN, actual_data[PLOT_STATE_IDX, :],
    label="Actual Data (X$PLOT_STATE_IDX)",
    linewidth=2, color=:black)
plot!(1:SIM_LEN, states_abc[PLOT_STATE_IDX, :],
    label="ABC Model (X$PLOT_STATE_IDX)",
    linewidth=1.5, linestyle=:dash, color=:red)
plot!(1:SIM_LEN, states_orig[PLOT_STATE_IDX, :],
    label="Original DMDc Model (X$PLOT_STATE_IDX)",
    linewidth=1.5, linestyle=:dot, color=:blue)

title!("Model Comparison (Delay = $DELAY)")
xlabel!("Time Step (relative to $sim_start_time)")
ylabel!("State Value")
plot!(legend=:outertopright)

# Save the plot
savefig("dmdc_comparison_delay$(DELAY)_start$(START_IDX).png")
println("Plot saved to dmdc_comparison_delay$(DELAY)_start$(START_IDX).png")

println("Script finished.")
