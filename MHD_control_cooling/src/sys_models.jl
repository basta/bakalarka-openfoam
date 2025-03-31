using LinearAlgebra, Logging

function create_linmodel(X_data, U_data)
    X_combined = [X_data[:, 1:end-1]; U_data[:, 1:end-1]; ones(1,size(X_data,2)-1)]
	X_next = X_data[:, 2:end]

	solution = X_next*pinv(X_combined)
	return (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
		return solution*[x;u;1], recur
	end
end

function create_delay_linmodel(X_data, U_data, shift)
	X_shifted_data = zeros(0, size(X_data,2)-shift-1)
	for i in 0:shift
		X_shifted_data = [
			X_shifted_data;
			X_data[:, 1+(shift-i):end-i-1]
		]
	end
	X_next = X_data[:, 1+1+shift:end]
	solution = X_next*pinv([X_shifted_data; U_data[:, 1+shift:end-1]; ones(1,size(X_shifted_data,2))])
	return (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
		if !haskey(recur, "delay_vecs")
			recur["delay_vecs"] = repeat(x, 1, shift)
		end
		next_state = solution*[x;vec(recur["delay_vecs"]);u;1] 
		if shift > 0
			recur["delay_vecs"][:, 2:end] = recur["delay_vecs"][:, 1:end-1]
			recur["delay_vecs"][:, 1] = x
		end
		return next_state, recur
	end
end

"""
    create_delay_step_linmodel(X_data, U_data, shift::Int, N::Int)

Identifies a linear model using delay embeddings with a step N.
The model predicts x(k+1) based on:
x(k), x(k-N), x(k-2N), ..., x(k-shift*N), u(k), and a constant offset.

Args:
- X_data (Matrix{Float64}): State trajectory data (n_x rows, T columns).
- U_data (Matrix{Float64}): Input trajectory data (n_u rows, T columns).
- shift (Int): The number of delay steps to include (e.g., shift=2 means using x(k-N) and x(k-2N)).
- N (Int): The step size for delays (must be >= 1).

Returns:
- Function: A model function `(x, u, recur_dict) -> (next_state, updated_recur_dict)`.
- Dict: The `recur_dict` stores the necessary history (`delay_vecs`).
"""
function create_delay_step_linmodel(X_data, U_data, shift::Int, N::Int)
    @assert size(X_data, 2) == size(U_data, 2) "X_data and U_data must have the same number of columns (time steps)."
    @assert shift >= 0 "Shift must be non-negative."
    @assert N >= 1 "Step N must be >= 1."

    T = size(X_data, 2)
    n_x = size(X_data, 1)
    n_u = size(U_data, 1) # Get number of inputs

    # Maximum delay needed in the history (e.g., if shift=2, N=3, need x(k-6))
    max_delay = shift * N

    # Number of columns available for the training matrices
    num_cols = T - 1 - max_delay
    if num_cols <= 0
        error("Not enough data points (T = $T) for the given shift (shift = $shift) and step (N = $N). Need T > shift*N + 1.")
    end

    # Build the matrix of stacked selected delayed states [x(k); x(k-N); ...; x(k-shift*N)]
    # Size: (n_x * (shift+1)) x num_cols
    X_stepped_data = zeros(n_x * (shift + 1), num_cols)
    for i in 0:shift # i=0 is x(k), i=1 is x(k-N), ..., i=shift is x(k-shift*N)
        current_delay = i * N
        row_start = 1 + i * n_x
        row_end = (i + 1) * n_x

        # Time indices for x(k - current_delay): start_idx to end_idx
        # k runs from (1+max_delay) to (T-1)
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay

        X_stepped_data[row_start:row_end, :] = X_data[:, start_idx:end_idx]
    end

    # Target states x(k+1)
    # k runs from (1+max_delay) to (T-1), so k+1 runs from (2+max_delay) to T
    X_next = X_data[:, 2+max_delay : end] # Indices: 2+max_delay to T

    # Inputs u(k)
    # k runs from (1+max_delay) to (T-1)
    U_relevant = U_data[:, 1+max_delay : T-1] # Indices: 1+max_delay to T-1

    # Combined matrix for regression: [X_stepped_data; U(k); 1]
    combined_input = [X_stepped_data; U_relevant; ones(1, num_cols)]

    try
        solution = X_next * pinv(combined_input)

        # Return the prediction function
        return (x::Vector{Float64}, u::Vector{Float64}, recur_dict::Dict) -> begin
           # Initialize history if first call. Store full history up to max_delay.
           if !haskey(recur_dict, "delay_vecs")
              # Stores [x(k-1) x(k-2) ... x(k-max_delay)]
              # Size: n_x x max_delay
              recur_dict["delay_vecs"] = repeat(x, 1, max_delay)
           end

           # Select the necessary delayed states: x(k-N), x(k-2N), ..., x(k-shift*N)
           # These are columns N, 2N, ..., shift*N from the stored history.
           selected_delays = if max_delay > 0 && shift > 0
                               # Use range with step N. Handles edge case if shift=0 (range is empty)
                               recur_dict["delay_vecs"][:, N:N:max_delay] # Size: n_x x shift
                           else
                               zeros(n_x, 0) # Empty matrix if no delays needed
                           end

           # Flatten the selected delays for the feature vector
           delay_history_vec = vec(selected_delays) # Size: (n_x * shift) x 1

           # Construct feature vector [x(k); x(k-N); ...; x(k-shift*N); u(k); 1]
           feature_vec = [x; delay_history_vec; u; 1.0]

           # Predict next state
           next_state = solution * feature_vec

           # Update *full* history for the *next* step (if max_delay > 0)
           if max_delay > 0
              history_matrix = recur_dict["delay_vecs"]
              # Shift history: columns 1 to max_delay-1 move to 2 to max_delay
              history_matrix[:, 2:end] = history_matrix[:, 1:end-1]
              # Insert current state x(k) into the first column (becomes x(k-1) next time)
              history_matrix[:, 1] = x
              recur_dict["delay_vecs"] = history_matrix # Update dict (might not be needed if mutable)
           end

           return next_state, recur_dict
        end
    catch e
        @error "Linear regression failed in create_delay_step_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end



function simulate_model(x₀::Vector{Float64}, U::Matrix{Float64}, model_fn)
    init_recur = Dict() # Initialize empty recurrence dictionary
    n_x = length(x₀)
    T_sim = size(U, 2)
    states = zeros(n_x, T_sim)
    state = copy(x₀) # Use copy to avoid modifying x₀

    for i in 1:T_sim
       u = vec(U[:, i]) # Ensure u is a vector
       states[:, i] = state # Store current state *before* prediction

       # Predict next state and update recurrence dictionary
       state, init_recur = model_fn(state, u, init_recur)
    end
    # Note: Simulation runs for T_sim steps, generating T_sim state predictions.
    # If you need the state *after* the last input, you might run one more step
    # or adjust indexing depending on convention (predict x(k+1) from u(k) or u(k+1)).
    # This implementation stores x(k) when using u(k) to predict x(k+1).
    return states
end

function simulate_model(x₀::Vector{Float64}, U::Matrix{Float64}, model_fn)
	init_recur = Dict()
	states = zeros(length(x₀), size(U,2))
	state = x₀
	for (i) in 1:size(U,2)
		u = U[:, i]
		states[:, i] = state
		state, init_recur = model_fn(state, vec(u), init_recur)
	end
	return states
end