using LinearAlgebra, Logging, Statistics

include("pinverses.jl")
include("dim_reduction.jl")

"""
    create_linmodel(X_data, U_data)

Identifies a simple linear model x(k+1) = A*x(k) + B*u(k) + c.
Returns a tuple: (prediction_function, initialization_length).
Initialization length is 0 for this model.
"""
function create_linmodel(X_data, U_data)
    n_x = size(X_data, 1)
    n_u = size(U_data, 1)
    T = size(X_data, 2)

    if T < 2
        error("Need at least 2 data points for linear model identification.")
    end

    # Prepare matrices for linear regression: X_next = solution * [X_current; U_current; 1]
    X_combined = [X_data[:, 1:end-1]; U_data[:, 1:end-1]; ones(1, T - 1)]
    X_next = X_data[:, 2:end]

    try
        solution = X_next * pinv(X_combined) # Calculate model parameters (A, B, c combined)

        # Define the prediction function
        model_fn = (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
            # No recurrence/history needed for this model
            next_state = solution * [x; u; 1.0]
            return next_state, recur # Return next state and unchanged recurrence dict
        end

        initialization_length = 0 # No initialization needed

        return model_fn, initialization_length

    catch e
        @error "Linear regression failed in create_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end

function create_static_model(X_data, U_data, max_order)
    n_x = size(X_data, 1)
    n_u = size(U_data, 1)
    T = size(X_data, 2)

    if T < 2
        error("Need at least 2 data points for linear model identification.")
    end
    X_next = X_data[:, 2:end]
    X_data = X_data[:, :]
    X_data .= 0

    for order in 2:max_order
        U_data = [
            U_data;
            U_data[1:n_u, :].^order
        ]
    end

    # Prepare matrices for linear regression: X_next = solution * [X_current; U_current; 1]
    X_combined = [X_data[:, 1:end-1]; U_data[:, 1:end-1]; ones(1, T - 1)]

    try
        solution = X_next * pinv(X_combined) # Calculate model parameters (A, B, c combined)
        @info "Average error is MSE $(mean((abs.(X_next - solution * X_combined))))"
        # Define the prediction function
        model_fn = (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
            # No recurrence/history needed for this model
            for order in 2:max_order
                u = [
                    u;
                    u[1:n_u].^order
                ]
            end
            next_state = solution * [x; u; 1.0]
            return next_state, recur # Return next state and unchanged recurrence dict
        end

        initialization_length = 0 # No initialization needed

        return model_fn, initialization_length

    catch e
        @error "Linear regression failed in create_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end


"""
    create_delay_linmodel(X_data, U_data, shift)

Identifies a linear model using delay embeddings:
x(k+1) = A0*x(k) + A1*x(k-1) + ... + Ashift*x(k-shift) + B*u(k) + c.
Returns a tuple: (prediction_function, initialization_length).
Initialization length is `shift`.
"""
function create_delay_linmodel(X_data, U_data, shift::Int)
    @assert shift >= 0 "Shift must be non-negative."
    n_x = size(X_data, 1)
    n_u = size(U_data, 1)
    T = size(X_data, 2)

    # Check if enough data points are available
    min_data_points = shift + 2 # Need x(k-shift) up to x(k) to predict x(k+1)
    if T < min_data_points
        error("Not enough data points (T = $T) for the given shift (shift = $shift). Need T >= $(min_data_points).")
    end

    num_cols = T - 1 - shift # Number of training examples

    # Build the matrix of stacked current and delayed states [x(k); x(k-1); ...; x(k-shift)]
    # Size: (n_x * (shift+1)) x num_cols
    X_shifted_data = zeros(n_x * (shift + 1), num_cols)
    for i in 0:shift # i=0 is x(k), i=1 is x(k-1), ..., i=shift is x(k-shift)
        row_start = 1 + i * n_x
        row_end = (i + 1) * n_x
        # Time indices for x(k - i): k runs from (1+shift) to (T-1)
        # So, (k-i) runs from (1+shift-i) to (T-1-i)
        start_idx = 1 + shift - i
        end_idx = T - 1 - i
        X_shifted_data[row_start:row_end, :] = X_data[:, start_idx:end_idx]
    end

    # Target states x(k+1)
    # k runs from (1+shift) to (T-1), so k+1 runs from (2+shift) to T
    X_next = X_data[:, (2+shift):end]

    # Inputs u(k)
    # k runs from (1+shift) to (T-1)
    U_relevant = U_data[:, (1+shift):(T-1)]

    # Combined matrix for regression: [X_shifted_data; U(k); 1]
    combined_input = [X_shifted_data; U_relevant; ones(1, num_cols)]

    try
        solution = X_next * pinv(combined_input)

        # Define the prediction function
        model_fn = (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
            # Initialize history if first call (or if history is missing)
            if !haskey(recur, "delay_vecs") || size(recur["delay_vecs"], 2) != shift
                # Stores [x(k-1) x(k-2) ... x(k-shift)]
                # Size: n_x x shift
                # Initialize with the current state `x` assuming it persists backwards.
                recur["delay_vecs"] = repeat(x, 1, shift)
                @debug "Initialized delay_vecs history (size: $(size(recur["delay_vecs"])))"
            end

            history_matrix = recur["delay_vecs"] # Size: n_x x shift
            history_vec = vec(history_matrix)    # Flattened history: Size (n_x * shift) x 1

            # Construct feature vector [x(k); x(k-1); ...; x(k-shift); u(k); 1]
            feature_vec = [x; history_vec; u; 1.0]

            # Predict next state
            next_state = solution * feature_vec

            # Update history for the *next* step (if shift > 0)
            if shift > 0
                # Shift history: columns 1 to shift-1 move to 2 to shift
                history_matrix[:, 2:end] = history_matrix[:, 1:end-1]
                # Insert current state x(k) into the first column (becomes x(k-1) next time)
                history_matrix[:, 1] = x
                recur["delay_vecs"] = history_matrix # Update dict
            end

            return next_state, recur
        end

        initialization_length = shift # Need 'shift' steps to fill the history

        return model_fn, initialization_length

    catch e
        @error "Linear regression failed in create_delay_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end


"""
    create_delay_step_linmodel(X_data, U_data, shift::Int, N::Int)

Identifies a linear model using delay embeddings with a step N:
x(k+1) = A0*x(k) + A1*x(k-N) + ... + Ashift*x(k-shift*N) + B*u(k) + c.
Returns a tuple: (prediction_function, initialization_length).
Initialization length is `shift * N`.
"""
function create_delay_step_linmodel(X_data, U_data, shift::Int, N::Int; pinv_method=solve_fastdmd)
    @assert size(X_data, 2) == size(U_data, 2) "X_data and U_data must have the same number of columns (time steps)."
    @assert shift >= 0 "Shift must be non-negative."
    @assert N >= 1 "Step N must be >= 1."

    T = size(X_data, 2)
    n_x = size(X_data, 1)
    n_u = size(U_data, 1)

    # Maximum delay needed in the history (e.g., if shift=2, N=3, need x(k-6))
    max_delay = shift * N

    # Check if enough data points are available
    min_data_points = max_delay + 2 # Need data up to x(k) and x(k-max_delay) to predict x(k+1)
    if T < min_data_points
        error("Not enough data points (T = $T) for the given shift (shift = $shift) and step (N = $N). Need T >= $(min_data_points).")
    end

    num_cols = T - 1 - max_delay # Number of training examples

    # Build the matrix of stacked selected delayed states [x(k); x(k-N); ...; x(k-shift*N)]
    # Size: (n_x * (shift+1)) x num_cols
    X_stepped_data = zeros(n_x * (shift + 1), num_cols)
    for i in 0:shift # i=0 is x(k), i=1 is x(k-N), ..., i=shift is x(k-shift*N)
        current_delay = i * N
        row_start = 1 + i * n_x
        row_end = (i + 1) * n_x

        # Time indices for x(k - current_delay): k runs from (1+max_delay) to (T-1)
        # So, (k - current_delay) runs from (1+max_delay-current_delay) to (T-1-current_delay)
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay

        X_stepped_data[row_start:row_end, :] = X_data[:, start_idx:end_idx]
    end

    # Target states x(k+1)
    # k runs from (1+max_delay) to (T-1), so k+1 runs from (2+max_delay) to T
    X_next = X_data[:, (2+max_delay):end]

    # Inputs u(k)
    # k runs from (1+max_delay) to (T-1)
    U_relevant = U_data[:, (1+max_delay):(T-1)]

    # Combined matrix for regression: [X_stepped_data; U(k); 1]
    combined_input = [X_stepped_data; U_relevant; ones(1, num_cols)]

    try
        # solution = X_next * pinv(combined_input)
        @info "
        Calculating pseudoinverse for linear model with state size $(size(combined_input, 1))
        X dimension:$(size(X_stepped_data, 1)) U dimension: $(size(U_relevant, 1))       
        "
        solution = pinv_method(combined_input, X_next)
        @info "Done"
        X_pred = solution * combined_input
        Residuals = X_next - X_pred
        MSE = mean(Residuals.^2)
        RMSE = sqrt(MSE)
        RMSE_per_state = sqrt.(mean(Residuals.^2, dims=2))
        @info RMSE RMSE_per_state
        # Define the prediction function
        model_fn = (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
           # Initialize full history if first call or if history is missing/wrong size.
           if !haskey(recur, "delay_vecs") || size(recur["delay_vecs"], 2) != max_delay
              # Stores [x(k-1) x(k-2) ... x(k-max_delay)]
              # Size: n_x x max_delay
              # Initialize with the current state `x` assuming it persists backwards.
              recur["delay_vecs"] = repeat(x, 1, max_delay)
              @debug "Initialized delay_vecs history (size: $(size(recur["delay_vecs"])))"
           end

           # Select the necessary delayed states: x(k-N), x(k-2N), ..., x(k-shift*N)
           # These are columns N, 2N, ..., shift*N from the stored history.
           selected_delays = if max_delay > 0 && shift > 0 && N > 0
                               # Use range with step N. Handles edge case if shift=0 (range is empty)
                               indices = N:N:min(max_delay, size(recur["delay_vecs"], 2)) # Ensure indices are valid
                               if isempty(indices)
                                   zeros(n_x, 0)
                               else
                                   recur["delay_vecs"][:, indices] # Size: n_x x shift (or fewer if history not full yet)
                               end
                           else
                               zeros(n_x, 0) # Empty matrix if no delays needed
                           end

           # Flatten the selected delays for the feature vector
           delay_history_vec = vec(selected_delays) # Size: (n_x * shift) x 1 (potentially smaller during init)

           # Construct feature vector [x(k); x(k-N); ...; x(k-shift*N); u(k); 1]
           feature_vec = [x; delay_history_vec; u; 1.0]

           # Predict next state
           next_state = solution * feature_vec

           # Update *full* history for the *next* step (if max_delay > 0)
           if max_delay > 0
              history_matrix = recur["delay_vecs"]
              # Shift history: columns 1 to max_delay-1 move to 2 to max_delay
              if size(history_matrix, 2) > 1 # Check if there's anything to shift
                 history_matrix[:, 2:end] = history_matrix[:, 1:end-1]
              end
              # Insert current state x(k) into the first column (becomes x(k-1) next time)
              if size(history_matrix, 2) >= 1 # Check if there's a column to insert into
                 history_matrix[:, 1] = x
              end
              recur["delay_vecs"] = history_matrix # Update dict
           end

           return next_state, recur
        end

        initialization_length = max_delay # Need 'max_delay' steps to fill the history

        return model_fn, initialization_length

    catch e
        @error "Linear regression failed in create_delay_step_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end

function create_delay_step_linmodel_pod(X_data, U_data, shift::Int, N::Int, n_modes::Int; pinv_method=solve_explicit)
    @assert size(X_data, 2) == size(U_data, 2) "X_data and U_data must have the same number of columns (time steps)."
    @assert shift >= 0 "Shift must be non-negative."
    @assert N >= 1 "Step N must be >= 1."

    T = size(X_data, 2)
    n_x = size(X_data, 1)
    n_u = size(U_data, 1)

    # Maximum delay needed in the history (e.g., if shift=2, N=3, need x(k-6))
    max_delay = shift * N

    # Check if enough data points are available
    min_data_points = max_delay + 2 # Need data up to x(k) and x(k-max_delay) to predict x(k+1)
    if T < min_data_points
        error("Not enough data points (T = $T) for the given shift (shift = $shift) and step (N = $N). Need T >= $(min_data_points).")
    end

    num_cols = T - 1 - max_delay # Number of training examples

    # Build the matrix of stacked selected delayed states [x(k); x(k-N); ...; x(k-shift*N)]
    # Size: (n_x * (shift+1)) x num_cols
    X_stepped_data = zeros(n_x * (shift + 1), num_cols)
    @info "Constructing delay matrix with delay $shift and N $N"
    for i in 0:shift # i=0 is x(k), i=1 is x(k-N), ..., i=shift is x(k-shift*N)
        current_delay = i * N
        row_start = 1 + i * n_x
        row_end = (i + 1) * n_x

        # Time indices for x(k - current_delay): k runs from (1+max_delay) to (T-1)
        # So, (k - current_delay) runs from (1+max_delay-current_delay) to (T-1-current_delay)
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay

        X_stepped_data[row_start:row_end, :] = X_data[:, start_idx:end_idx]
    end

    modes, s, coordinates, _,  mean_vec = pod_pca(X_stepped_data, center=true, num_modes=n_modes)
    @info "Singular values are" s

    # Target states x(k+1)
    # k runs from (1+max_delay) to (T-1), so k+1 runs from (2+max_delay) to T
    X_next = X_data[:, (2+max_delay):end]

    # Inputs u(k)
    # k runs from (1+max_delay) to (T-1)
    U_relevant = U_data[:, (1+max_delay):(T-1)]

    # Combined matrix for regression: [X_stepped_data; U(k); 1]
    combined_input = [coordinates; U_relevant; ones(1, num_cols)]

    try
        # solution = X_next * pinv(combined_input)
        @info "
        Calculating pseudoinverse for linear model with state size $(size(combined_input, 1))
        X dimension:$(size(X_stepped_data, 1)) U dimension: $(size(U_relevant, 1))       
        "
        solution = pinv_method(combined_input, X_next)
        @info "Solution is" solution

        X_pred = solution * combined_input
        Residuals = X_next - X_pred
        MSE = mean(Residuals.^2)
        RMSE = sqrt(MSE)
        RMSE_per_state = sqrt.(mean(Residuals.^2, dims=2))
        @info RMSE RMSE_per_state
        # Define the prediction function
        model_fn = (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
           # Initialize full history if first call or if history is missing/wrong size.
           if !haskey(recur, "delay_vecs") || size(recur["delay_vecs"], 2) != max_delay
              # Stores [x(k-1) x(k-2) ... x(k-max_delay)]
              # Size: n_x x max_delay
              # Initialize with the current state `x` assuming it persists backwards.
              recur["delay_vecs"] = repeat(x, 1, max_delay)
              @debug "Initialized delay_vecs history (size: $(size(recur["delay_vecs"])))"
           end

           # Select the necessary delayed states: x(k-N), x(k-2N), ..., x(k-shift*N)
           # These are columns N, 2N, ..., shift*N from the stored history.
           selected_delays = if max_delay > 0 && shift > 0 && N > 0
                               # Use range with step N. Handles edge case if shift=0 (range is empty)
                               indices = N:N:min(max_delay, size(recur["delay_vecs"], 2)) # Ensure indices are valid
                               if isempty(indices)
                                   zeros(n_x, 0)
                               else
                                   recur["delay_vecs"][:, indices] # Size: n_x x shift (or fewer if history not full yet)
                               end
                           else
                               zeros(n_x, 0) # Empty matrix if no delays needed
                           end

           # Flatten the selected delays for the feature vector
           delay_history_vec = vec(selected_delays) # Size: (n_x * shift) x 1 (potentially smaller during init)

           # Construct feature vector [x(k); x(k-N); ...; x(k-shift*N); u(k); 1]
           feature_vec = [x; delay_history_vec;]
           feature_vec = modes' * (feature_vec - mean_vec)
           feature_vec = [feature_vec; u; 1.0]

           # Predict next state
           next_state = solution * feature_vec


           # Update *full* history for the *next* step (if max_delay > 0)
           if max_delay > 0
              history_matrix = recur["delay_vecs"]
              # Shift history: columns 1 to max_delay-1 move to 2 to max_delay
              if size(history_matrix, 2) > 1 # Check if there's anything to shift
                 history_matrix[:, 2:end] = history_matrix[:, 1:end-1]
              end
              # Insert current state x(k) into the first column (becomes x(k-1) next time)
              if size(history_matrix, 2) >= 1 # Check if there's a column to insert into
                 history_matrix[:, 1] = x
              end
              recur["delay_vecs"] = history_matrix # Update dict
           end

           return next_state, recur
        end

        initialization_length = max_delay # Need 'max_delay' steps to fill the history

        return model_fn, initialization_length

    catch e
        @error "Linear regression failed in create_delay_step_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end

function create_delay_step_linmodel_pod_first(X_data, U_data, shift::Int, N::Int, n_modes::Int; pinv_method=solve_fastdmd)
    @assert size(X_data, 2) == size(U_data, 2) "X_data and U_data must have the same number of columns (time steps)."
    @assert shift >= 0 "Shift must be non-negative."
    @assert N >= 1 "Step N must be >= 1."

    T = size(X_data, 2)
    original_n_x = size(X_data, 1) # Dimension of the original state space
    n_u = size(U_data, 1)

    # --- Apply POD to the full state data ---
    # modes: High-dim basis vectors (original_n_x x n_modes)
    # coordinates: Low-dim representation (n_modes x T)
    # mean_vec: Mean subtracted from X_data before POD (original_n_x x 1)
    modes, s, coordinates, _, mean_vec = pod_pca(X_data, center=true, num_modes=n_modes)
    @info "POD complete. Singular values:" s
    n_x = n_modes # Work with the dimension of the coordinate space internally

    # --- Setup for Delay Embedding ---
    max_delay = shift * N # Maximum delay needed in history

    # Check data availability
    min_data_points = max_delay + 2
    if T < min_data_points
        error("Not enough data points (T = $T) for shift=$shift, N=$N. Need T >= $min_data_points.")
    end

    num_cols = T - 1 - max_delay # Number of training examples

    # --- Build Delayed Coordinate Matrix ---
    # Stacked delayed POD coordinates [z(k); z(k-N); ...; z(k-shift*N)]
    # Size: (n_modes * (shift+1)) x num_cols
    X_stepped_coords_data = zeros(n_modes * (shift + 1), num_cols)
    @info "Constructing delay matrix for POD coordinates (delay $shift, N $N)"
    for i in 0:shift # i=0 is z(k), i=1 is z(k-N), ..., i=shift is z(k-shift*N)
        current_delay = i * N
        row_start = 1 + i * n_modes
        row_end = (i + 1) * n_modes

        # Time indices for z(k - current_delay): k runs from (1+max_delay) to (T-1)
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay

        # Use the POD coordinates calculated earlier
        X_stepped_coords_data[row_start:row_end, :] = coordinates[:, start_idx:end_idx]
    end

    # --- Target Coordinates ---
    # Target POD coordinates z(k+1)
    # k runs from (1+max_delay) to (T-1), so k+1 runs from (2+max_delay) to T
    X_next_coords = coordinates[:, (2+max_delay):end] # Target is now low-dimensional

    # --- Inputs ---
    # Inputs u(k) corresponding to the time steps where prediction starts
    # k runs from (1+max_delay) to (T-1)
    U_relevant = U_data[:, (1+max_delay):(T-1)]

    # --- Combined Input Matrix for Regression ---
    # [Delayed Coordinates; Relevant Inputs; Bias Term]
    combined_input = [X_stepped_coords_data; U_relevant; ones(1, num_cols)]

    try
        # --- Solve Linear System for Coordinate Dynamics ---
        # solution * combined_input ≈ X_next_coords
        @info """
        Calculating pseudoinverse for linear model in coordinate space.
        Input matrix size: $(size(combined_input))
        Target coordinate matrix size: $(size(X_next_coords))
        Coordinate dimension (n_modes): $n_modes
        Input dimension (n_u): $n_u
        Number of delay terms (shift+1): $(shift+1)
        """
        solution = pinv_method(combined_input, X_next_coords) # Predicts coordinates
        @info "Pseudoinverse calculation complete."

        # --- Evaluate Training Fit (in original high-dimensional space) ---
        X_pred_coords = solution * combined_input
        # Transform predicted coordinates back to high-dimensional space
        X_pred_full = modes * X_pred_coords .+ mean_vec
        # Get the corresponding true high-dimensional states
        X_next_full = X_data[:, (2+max_delay):end]
        Residuals = X_next_full - X_pred_full
        MSE = mean(Residuals.^2)
        RMSE = sqrt(MSE)
        RMSE_per_state = sqrt.(mean(Residuals.^2, dims=2))
        @info "Training RMSE (evaluated in original space): $RMSE"
        # @debug "Training RMSE per state (original space): $(vec(RMSE_per_state))"

        # --- Define the Prediction Function (Closure) ---
        model_fn = (x_high_dim::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
            # Project current high-dimensional state to POD coordinates
            # Ensure x_high_dim and mean_vec are column vectors for subtraction
            x_coords = modes' * (reshape(x_high_dim, :, 1) .- mean_vec)
            x_coords = vec(x_coords) # Make sure it's a vector

            # Initialize history buffer for *coordinates* if first call or wrong size
            if !haskey(recur, "delay_coords_vecs") || size(recur["delay_coords_vecs"], 1) != n_modes || size(recur["delay_coords_vecs"], 2) != max_delay
                # Stores [z(k-1) z(k-2) ... z(k-max_delay)]
                # Size: n_modes x max_delay
                recur["delay_coords_vecs"] = repeat(x_coords, 1, max_delay)
                @debug "Initialized delay_coords_vecs history (size: $(size(recur["delay_coords_vecs"])))"
            end

            # Select necessary delayed *coordinates* from history: z(k-N), ..., z(k-shift*N)
            selected_delay_coords = if max_delay > 0 && shift > 0 && N > 0
                indices = N:N:min(max_delay, size(recur["delay_coords_vecs"], 2))
                isempty(indices) ? zeros(n_modes, 0) : recur["delay_coords_vecs"][:, indices]
            else
                zeros(n_modes, 0) # Empty matrix if no delays needed
            end

            # Flatten selected delayed coordinates
            delay_history_coords_vec = vec(selected_delay_coords)

            # Construct feature vector for coordinate prediction
            # [z(k); z(k-N); ...; z(k-shift*N); u(k); 1]
            feature_vec = [
                x_coords;                # Current coordinates z(k)
                delay_history_coords_vec;# Delayed coordinates z(k-N)...
                u;                       # Current input u(k)
                1.0                      # Bias term
            ]

            # Predict *next coordinates*
            next_coords_pred = solution * feature_vec

            # --- Inverse Transform: Map predicted coordinates back to high-dimensional space ---
            next_state_high_dim = modes * next_coords_pred .+ vec(mean_vec) # Ensure mean_vec is added as a vector

            # --- Update History Buffer (with current coordinates) ---
            if max_delay > 0
                history_matrix = recur["delay_coords_vecs"]
                # Shift history
                if size(history_matrix, 2) > 1
                    history_matrix[:, 2:end] = history_matrix[:, 1:end-1]
                end
                # Insert current coordinates z(k)
                if size(history_matrix, 2) >= 1
                    history_matrix[:, 1] = x_coords
                end
                recur["delay_coords_vecs"] = history_matrix # Update dict
            end

            # Return predicted high-dimensional state and updated recurrence dict
            return vec(next_state_high_dim), recur
        end

        initialization_length = max_delay # Need 'max_delay' steps to fill history

        return model_fn, initialization_length

    catch e
        @error "Operation failed in create_delay_step_linmodel_pod_first." exception=(e, catch_backtrace())
        rethrow(e)
    end
end


function create_delay_step_linmodel_pod_first_with_input_delay(
    X_data, U_data, shift::Int, N::Int, n_modes::Int; 
    pinv_method=solve_fastdmd, skip_pod=false
    )
    @assert size(X_data, 2) == size(U_data, 2) "X_data and U_data must have the same number of columns (time steps)."
    @assert shift >= 0 "Shift must be non-negative."
    @assert N >= 1 "Step N must be >= 1."

    T = size(X_data, 2)
    original_n_x = size(X_data, 1)
    n_u = size(U_data, 1)

    # --- Apply POD only to the state data ---
    if skip_pod
        modes = I
        coordinates = X_data[:,:]
        mean_vec = mean(X_data, dims=2)
        coordinates .-= mean_vec
        @info "Skipping POD. Using identity matrix for modes."
    else
        modes, s, coordinates, _, mean_vec = pod_pca(X_data, center=true, num_modes=n_modes)
        @info "POD complete for state data. Singular values:" s

    end
    n_x_coords = n_modes # Dimension of the POD coordinate space

    # --- Setup for Delay Embedding ---
    max_delay = shift * N # Maximum delay needed in history for both state and input

    min_data_points = max_delay + 2
    if T < min_data_points
        error("Not enough data points (T = $T) for shift=$shift, N=$N. Need T >= $min_data_points.")
    end

    num_cols = T - 1 - max_delay # Number of training examples

    # --- Build Delayed Coordinate Matrix ---
    # Stacked delayed POD coordinates [z(k); z(k-N); ...; z(k-shift*N)]
    X_stepped_coords_data = zeros(n_x_coords * (shift + 1), num_cols)
    @info "Constructing delay matrix for POD coordinates (delay $shift, N $N)"
    for i in 0:shift
        current_delay = i * N
        row_start = 1 + i * n_x_coords
        row_end = (i + 1) * n_x_coords
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay
        X_stepped_coords_data[row_start:row_end, :] = coordinates[:, start_idx:end_idx]
    end

    # --- Build Delayed Input Matrix ---
    # Stacked delayed inputs [u(k); u(k-N); ...; u(k-shift*N)]
    U_stepped_data = zeros(n_u * (shift + 1), num_cols)
    @info "Constructing delay matrix for inputs (delay $shift, N $N)"
     for i in 0:shift
        current_delay = i * N
        row_start = 1 + i * n_u
        row_end = (i + 1) * n_u
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay
        U_stepped_data[row_start:row_end, :] = U_data[:, start_idx:end_idx]
    end

    # --- Target Coordinates ---
    # Target POD coordinates z(k+1)
    X_next_coords = coordinates[:, (2+max_delay):end]

    # --- Combined Input Matrix for Regression ---
    # [Delayed Coordinates; Delayed Inputs; Bias Term]
    combined_input = [X_stepped_coords_data; U_stepped_data; ones(1, num_cols)]

    try
        # --- Solve Linear System for Coordinate Dynamics ---
        @info """
        Calculating pseudoinverse for linear model in coordinate space.
        Input matrix size: $(size(combined_input))
        Target coordinate matrix size: $(size(X_next_coords))
        Coordinate dimension (n_modes): $n_x_coords
        Input dimension (n_u): $n_u
        Number of delay terms per variable (shift+1): $(shift+1)
        """
        solution = pinv_method(combined_input, X_next_coords) # Predicts next coordinates
        @info "Pseudoinverse calculation complete."

        # --- Evaluate Training Fit (optional, in original high-dimensional space) ---
        X_pred_coords = solution * combined_input
        X_pred_full = modes * X_pred_coords .+ mean_vec
        X_next_full = X_data[:, (2+max_delay):end]
        Residuals = X_next_full - X_pred_full
        RMSE = sqrt(mean(Residuals.^2))
        @info "Training RMSE (evaluated in original space): $RMSE"

        # --- Define the Prediction Function (Closure) ---
        model_fn = (x_high_dim::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
            # Project current high-dimensional state to POD coordinates
            x_coords = modes' * (reshape(x_high_dim, :, 1) .- mean_vec)
            x_coords = vec(x_coords) # Ensure it's a vector

            # Initialize history buffers if first call or wrong size
            # State coordinate history: [z(k-1) z(k-2) ... z(k-max_delay)]
            if !haskey(recur, "delay_coords_vecs") || size(recur["delay_coords_vecs"]) != (n_x_coords, max_delay)
                recur["delay_coords_vecs"] = repeat(x_coords, 1, max_delay)
                @debug "Initialized delay_coords_vecs history (size: $(size(recur["delay_coords_vecs"])))"
            end
            # Input history: [u(k-1) u(k-2) ... u(k-max_delay)]
             if !haskey(recur, "delay_u_vecs") || size(recur["delay_u_vecs"]) != (n_u, max_delay)
                recur["delay_u_vecs"] = repeat(reshape(u,:,1), 1, max_delay) # Use current u for init
                @debug "Initialized delay_u_vecs history (size: $(size(recur["delay_u_vecs"])))"
            end

            # Select necessary delayed states (coordinates) and inputs from history
            delay_indices = N:N:min(max_delay, size(recur["delay_coords_vecs"], 2)) # Same indices for both

            selected_delay_coords = if max_delay > 0 && shift > 0 && N > 0 && !isempty(delay_indices)
                 recur["delay_coords_vecs"][:, delay_indices]
            else
                zeros(n_x_coords, 0)
            end

            selected_delay_u = if max_delay > 0 && shift > 0 && N > 0 && !isempty(delay_indices)
                 recur["delay_u_vecs"][:, delay_indices]
            else
                zeros(n_u, 0)
            end

            # Flatten selected delayed variables
            delay_history_coords_vec = vec(selected_delay_coords)
            delay_history_u_vec = vec(selected_delay_u)

            # Construct feature vector for coordinate prediction
            # [z(k); z(k-N); ...; z(k-shift*N); u(k); u(k-N); ...; u(k-shift*N); 1]
            feature_vec = [
                x_coords;                 # Current coordinates z(k)
                delay_history_coords_vec; # Delayed coordinates z(k-N)...
                u;                        # Current input u(k)
                delay_history_u_vec;      # Delayed inputs u(k-N)...
                1.0                       # Bias term
            ]

            # Predict *next coordinates*
            next_coords_pred = solution * feature_vec

            # --- Inverse Transform: Map predicted coordinates back to high-dimensional space ---
            next_state_high_dim = modes * next_coords_pred .+ vec(mean_vec)

            # --- Update History Buffers (with current state coordinates and input) ---
            if max_delay > 0
                # Update state coordinate history
                history_coords = recur["delay_coords_vecs"]
                if size(history_coords, 2) > 1
                    history_coords[:, 2:end] = history_coords[:, 1:end-1]
                end
                if size(history_coords, 2) >= 1
                    history_coords[:, 1] = x_coords # Insert current z(k)
                end
                recur["delay_coords_vecs"] = history_coords

                # Update input history
                history_u = recur["delay_u_vecs"]
                 if size(history_u, 2) > 1
                    history_u[:, 2:end] = history_u[:, 1:end-1]
                end
                if size(history_u, 2) >= 1
                    history_u[:, 1] = u # Insert current u(k)
                end
                recur["delay_u_vecs"] = history_u
            end

            return vec(next_state_high_dim), recur
        end

        initialization_length = max_delay # Need 'max_delay' steps to fill history

        return model_fn, initialization_length, solution

    catch e
        @error "Operation failed in create_delay_step_linmodel_pod_first_with_input_delay." exception=(e, catch_backtrace())
        rethrow(e)
    end
end


function standardize_data(data::AbstractArray{Float64, 2})
    # Calculate mean for each feature (row)
    # Result is a column vector
    means = vec(mean(data, dims=2))

    # Calculate standard deviation for each feature (row)
    # Result is a column vector
    stds = vec(std(data, dims=2))

    # Handle features with zero standard deviation (constant features)
    # Avoid division by zero. Replace std=0 with 1.0.
    # The normalized value for these features will be (feature_value - mean) / 1.0 = 0.0 / 1.0 = 0.0
    stds_corrected = copy(stds)
    zero_std_indices = findall(isapprox.(stds, 0.0, atol=1e-10)) # Find indices where std is close to zero
    if !isempty(zero_std_indices)
        stds_corrected[zero_std_indices] .= 1.0
        @warn "Features (rows) with zero standard deviation found at indices: $zero_std_indices. Their normalized value will be 0. Std dev returned as 1.0 for these."
    end

    # Apply standardization: (data - mean) / std_dev
    # Broadcasting takes care of applying column vectors `means` and `stds_corrected` to each column of `data`
    normalized_data = (data .- means) ./ stds_corrected

    return normalized_data, means, stds_corrected
end

function apply_standardization(new_data::Matrix{Float64}, means::Vector{Float64}, stds_corrected::Vector{Float64})
    # Ensure dimensions match
    n_features_data = size(new_data, 1)
    n_features_params = length(means)
    @assert n_features_data == n_features_params "Number of features (rows) in new_data ($n_features_data) must match the length of the means vector ($n_features_params)."
    @assert n_features_data == length(stds_corrected) "Number of features (rows) in new_data ($n_features_data) must match the length of the stds_corrected vector."
    @assert all(stds_corrected .> 0) "Corrected standard deviations must all be positive (zero stds should have been replaced by 1.0)."

    # Apply the standardization using the provided means and corrected stds
    normalized_new_data = (new_data .- means) ./ stds_corrected

    return normalized_new_data
end

function invert_standardization(normalized_data::Matrix{Float64}, means::Vector{Float64}, stds_corrected::Vector{Float64})
    # Ensure dimensions match
    n_features_data = size(normalized_data, 1)
    n_features_params = length(means)
    @assert n_features_data == n_features_params "Number of features (rows) in normalized_data ($n_features_data) must match the length of the means vector ($n_features_params)."
    @assert n_features_data == length(stds_corrected) "Number of features (rows) in normalized_data ($n_features_data) must match the length of the stds_corrected vector."
    @assert all(stds_corrected .> 0) "Corrected standard deviations must all be positive."

    # Apply the inverse transformation: original = (normalized * std_dev) + mean
    # Broadcasting handles applying the vectors element-wise to the matrix columns
    original_scale_data = (normalized_data .* stds_corrected) .+ means

    return original_scale_data
end

function create_delay_step_linmodel_normalized(X_data, U_data, shift::Int, N::Int; pinv_method=solve_explicit)
    @assert size(X_data, 2) == size(U_data, 2) "X_data and U_data must have the same number of columns (time steps)."
    @assert shift >= 0 "Shift must be non-negative."
    @assert N >= 1 "Step N must be >= 1."

    X_data, means, stds = standardize_data(X_data)
    U_data, u_means, u_stds = standardize_data(U_data)

    T = size(X_data, 2)
    n_x = size(X_data, 1)
    n_u = size(U_data, 1)

    # Maximum delay needed in the history (e.g., if shift=2, N=3, need x(k-6))
    max_delay = shift * N

    # Check if enough data points are available
    min_data_points = max_delay + 2 # Need data up to x(k) and x(k-max_delay) to predict x(k+1)
    if T < min_data_points
        error("Not enough data points (T = $T) for the given shift (shift = $shift) and step (N = $N). Need T >= $(min_data_points).")
    end

    num_cols = T - 1 - max_delay # Number of training examples

    # Build the matrix of stacked selected delayed states [x(k); x(k-N); ...; x(k-shift*N)]
    # Size: (n_x * (shift+1)) x num_cols
    X_stepped_data = zeros(n_x * (shift + 1), num_cols)
    for i in 0:shift # i=0 is x(k), i=1 is x(k-N), ..., i=shift is x(k-shift*N)
        current_delay = i * N
        row_start = 1 + i * n_x
        row_end = (i + 1) * n_x

        # Time indices for x(k - current_delay): k runs from (1+max_delay) to (T-1)
        # So, (k - current_delay) runs from (1+max_delay-current_delay) to (T-1-current_delay)
        start_idx = 1 + max_delay - current_delay
        end_idx = T - 1 - current_delay

        X_stepped_data[row_start:row_end, :] = X_data[:, start_idx:end_idx]
    end

    # Target states x(k+1)
    # k runs from (1+max_delay) to (T-1), so k+1 runs from (2+max_delay) to T
    X_next = X_data[:, (2+max_delay):end]

    # Inputs u(k)
    # k runs from (1+max_delay) to (T-1)
    U_relevant = U_data[:, (1+max_delay):(T-1)]

    # Combined matrix for regression: [X_stepped_data; U(k); 1]
    combined_input = [X_stepped_data; U_relevant; ones(1, num_cols)]

    try
        # solution = X_next * pinv(combined_input)
        @info "
        Calculating pseudoinverse for linear model with state size $(size(combined_input, 1))
        X dimension:$(size(X_stepped_data, 1)) U dimension: $(size(U_relevant, 1))       
        "
        solution = pinv_method(combined_input, X_next)
        @info "Done"
        X_pred = solution * combined_input
        Residuals = X_next - X_pred
        MSE = mean(Residuals.^2)
        RMSE = sqrt(MSE)
        RMSE_per_state = sqrt.(mean(Residuals.^2, dims=2))
        @info RMSE RMSE_per_state
        # Define the prediction function
        model_fn = (x::Vector{Float64}, u::Vector{Float64}, recur::Dict) -> begin
           # Initialize full history if first call or if history is missing/wrong size.
           if !haskey(recur, "delay_vecs") || size(recur["delay_vecs"], 2) != max_delay
              # Stores [x(k-1) x(k-2) ... x(k-max_delay)]
              # Size: n_x x max_delay
              # Initialize with the current state `x` assuming it persists backwards.
              recur["delay_vecs"] = repeat(x, 1, max_delay)
              @debug "Initialized delay_vecs history (size: $(size(recur["delay_vecs"])))"
           end
           @info x
           x = apply_standardization(reshape(x, :, 1), means, stds)[:, 1]
           u = apply_standardization(reshape(u, :, 1), u_means, u_stds)[:, 1]
           @info x
           # Select the necessary delayed states: x(k-N), x(k-2N), ..., x(k-shift*N)
           # These are columns N, 2N, ..., shift*N from the stored history.
           selected_delays = if max_delay > 0 && shift > 0 && N > 0
                               # Use range with step N. Handles edge case if shift=0 (range is empty)
                               indices = N:N:min(max_delay, size(recur["delay_vecs"], 2)) # Ensure indices are valid
                               if isempty(indices)
                                   zeros(n_x, 0)
                               else
                                   recur["delay_vecs"][:, indices] # Size: n_x x shift (or fewer if history not full yet)
                               end
                           else
                               zeros(n_x, 0) # Empty matrix if no delays needed
                           end

           # Flatten the selected delays for the feature vector
           delay_history_vec = vec(selected_delays) # Size: (n_x * shift) x 1 (potentially smaller during init)

           # Construct feature vector [x(k); x(k-N); ...; x(k-shift*N); u(k); 1]
           feature_vec = [x; delay_history_vec; u; 1.0]

           # Predict next state
           next_state = solution * feature_vec

           # Update *full* history for the *next* step (if max_delay > 0)
           if max_delay > 0
              history_matrix = recur["delay_vecs"]
              # Shift history: columns 1 to max_delay-1 move to 2 to max_delay
              if size(history_matrix, 2) > 1 # Check if there's anything to shift
                 history_matrix[:, 2:end] = history_matrix[:, 1:end-1]
              end
              # Insert current state x(k) into the first column (becomes x(k-1) next time)
              if size(history_matrix, 2) >= 1 # Check if there's a column to insert into
                 history_matrix[:, 1] = x
              end
              recur["delay_vecs"] = history_matrix # Update dict
           end

           return invert_standardization(reshape(next_state, :, 1), means, stds)[:, 1], recur
        end

        initialization_length = max_delay # Need 'max_delay' steps to fill the history

        return model_fn, initialization_length

    catch e
        @error "Linear regression failed in create_delay_step_linmodel. Check data rank/conditioning." e
        rethrow(e)
    end
end


"""
    simulate_model(x₀::Vector{Float64}, U::Matrix{Float64}, model_fn_tuple)

Simulates a model forward in time, handling an initialization phase.

Args:
- x₀ (Vector{Float64}): Initial state.
- U (Matrix{Float64}): Input trajectory (n_u rows, T columns). T includes initialization steps.
- model_fn_tuple (Tuple{Function, Int}): A tuple containing:
    - The model prediction function `(x, u, recur) -> (next_state, updated_recur)`.
    - The required initialization length (`init_len`).

Returns:
- Matrix{Float64}: Matrix of predicted states (n_x rows, T - init_len columns).
                   Returns an empty matrix if T <= init_len.
"""
function simulate_model(x₀::Vector{Float64}, U::Matrix{Float64}, model_fn_tuple, x_pred = nothing)
    model_fn, init_len = model_fn_tuple
    n_x = length(x₀)
    T_total = size(U, 2)

    if T_total < init_len
        @warn "Total simulation time (T_total = $T_total) is less than the required initialization length (init_len = $init_len). Cannot perform prediction."
        return zeros(n_x, 0) # Return empty matrix
    end

    recur_dict = Dict() # Initialize empty recurrence dictionary
    state = copy(x₀)    # Start with the initial state
    predicted_states = zeros(n_x, T_total) # Allocate matrix for results
    predicted_states[:, 1] = state 

    # --- Initialization Phase ---
    # Run the model for init_len steps just to populate the history (recur_dict)
    @debug "Starting initialization phase (length = $init_len)..."
    for t in 1:init_len
        u = vec(U[:, t]) # Ensure u is a vector
        # Predict next state and update recurrence dictionary
        next_state, recur_dict = model_fn(state, u, recur_dict)
        state = next_state # Update state for the next initialization step
        if !isnothing(x_pred)
            state = x_pred[:, t]
        end
        predicted_states[:, t] = state
        # We don't store states during initialization
    end
    @debug "Initialization phase complete. Starting state for prediction: $state"
    @debug "Recurrence dict after init: $recur_dict"

    # --- Prediction Phase ---
    T_pred = T_total - init_len # Number of actual prediction steps

    @debug "Starting prediction phase (length = $T_pred)..."
    for t in 1:T_pred
        sim_step_idx = init_len + t # Index into the U matrix for the current step
        u = vec(U[:, sim_step_idx]) # Get the input for the current prediction step

        # Store the *current* state before predicting the next one
        predicted_states[:, sim_step_idx] = state

        # Predict the next state and update recurrence dictionary
        next_state, recur_dict = model_fn(state, u, recur_dict)
        state = next_state # Update state for the next prediction step
    end
    @debug "Prediction phase complete."

    # The loop stores the state *at* time step k (using input u(k)) before calculating state k+1.
    # The size of predicted_states is T_pred = T_total - init_len.
    return predicted_states
end

function extract_dmdc_matrices(solution::Matrix{Float64}, n_x_coords::Int, n_u::Int, shift::Int)
    num_A_matrices = shift + 1
    num_B_matrices = shift + 1

    # --- Calculate expected dimensions ---
    expected_cols = n_x_coords * num_A_matrices + n_u * num_B_matrices + 1
    expected_rows = n_x_coords

    if size(solution) != (expected_rows, expected_cols)
        error("Dimension mismatch: Input solution matrix has size $(size(solution)), but expected ($expected_rows, $expected_cols) based on n_x_coords=$n_x_coords, n_u=$n_u, shift=$shift.")
    end

    # --- Initialize storage ---
    As = Vector{Matrix{Float64}}(undef, num_A_matrices)
    Bs = Vector{Matrix{Float64}}(undef, num_B_matrices)

    # --- Extract A matrices ---
    current_col = 1
    for i in 0:shift # Corresponds to A₀, A₁, ..., A_shift
        start_col = current_col
        end_col = current_col + n_x_coords - 1
        As[i+1] = solution[:, start_col:end_col]
        current_col = end_col + 1
    end

    # --- Extract B matrices ---
    for i in 0:shift # Corresponds to B₀, B₁, ..., B_shift
        start_col = current_col
        end_col = current_col + n_u - 1
        # Handle case where n_u might be 0 (no control input)
        if n_u > 0
             Bs[i+1] = solution[:, start_col:end_col]
             current_col = end_col + 1
        else
            # If no inputs, create empty matrices of correct size
             Bs[i+1] = zeros(n_x_coords, 0)
             # current_col remains unchanged as no columns were consumed
        end
    end

    # --- Extract C vector ---
    # C is the last column
    C = solution[:, end] # This automatically becomes a Vector

    return As, Bs, C
end