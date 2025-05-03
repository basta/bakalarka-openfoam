# src/controllers/blackbox_online.jl
# Controller using BlackBoxOptim.jl's ask-tell interface for online optimization.

using BlackBoxOptim
using LinearAlgebra
using TOML
using Logging
using Statistics # For potential criterion calculation
using JSON

# If using CircularBuffer for state history:
# using DataStructures

# Assuming AbstractController and interface functions are defined in interface.jl
include("interface.jl")

# --- Helper Function for Criterion Calculation (Example) ---
# You MUST tailor this function to your specific needs.
# It calculates the performance score based on the state observed at the END of the evaluation period.
function default_end_period_criterion(state_at_end::Dict{Symbol, Any})::Float64
    fluxes = get(state_at_end, :wallheatflux, Float32[])
    @info fluxes
    if length(fluxes) >= 2
        # Example: minimize heat going INTO cold wall (usually negative, so maximize it towards 0)
        # return -Float64(fluxes[2]) # Maximize negative value means minimize magnitude
        # Example: Minimize total heat removed (sum is usually negative)
        return Float64(fluxes[1]) # Minimize the (negative) sum
    else
        @warn "Cannot calculate criterion: :wallheatflux data missing or insufficient."
        return Inf
    end
end

# --- Controller Definition ---

mutable struct BlackBoxOnlineController <: AbstractController
    search_range::Vector{Tuple{Float64, Float64}} # Lower/Upper bounds for each input dim
    T_eval::Float64                              # Evaluation duration for each candidate
    num_inputs::Int                              # Dimension of the control vector
    num_corners::Int                             # 2^num_inputs
    criterion_calculator::Function               # Function (state_dict) -> score
    minimize_criterion::Bool                     # True if lower score is better
    output_file_path::String                     # File path to save the best input

    best_input_so_far::Vector{Float64}           # Best corner input found over time
    best_score_so_far::Float64                   # Score associated with best_input_so_far

    current_corner_index::Int                    # Index (0 to num_corners-1) representing the corner being tested
    current_candidate_input::Union{Vector{Float64}, Nothing} # Corner input currently being evaluated
    current_eval_start_time::Float64             # Sim time when current evaluation began
    is_evaluating::Bool                          # Flag: currently evaluating a candidate?
    cycle_complete::Bool                         # Flag: Have we tested all corners at least once?

    function BlackBoxOnlineController(config::Dict)
        cso_params = config["controller"]["params"]["BlackBoxOnline"]

        num_inputs = get(cso_params, "num_inputs", 8)
        if num_inputs > 20 # Prevent very large number of corners
            @warn "num_inputs ($num_inputs) is large, resulting in $(2^num_inputs) corners. Consider reducing inputs or using random corner sampling."
        end
        num_corners = 2^num_inputs

        T_eval = Float64(get(cso_params, "T_eval", 100.0))
        minimize_criterion = get(cso_params, "minimize_criterion", true)

        # --- Get Search Range (Mandatory) ---
        search_range_tuple = get(cso_params, "SearchRange", nothing)
        if isnothing(search_range_tuple)
            error("Missing mandatory 'SearchRange' parameter in BlackBoxOnline config.")
        end
        # Ensure SearchRange is interpreted correctly for num_inputs
        if length(search_range_tuple) != num_inputs
            if length(search_range_tuple) == 2 && num_inputs > 1
                @warn "SearchRange appears to be a single range for all inputs. Applying [$(search_range_tuple[1]), $(search_range_tuple[2])] to all $num_inputs dimensions."
                 search_range = [(Float64(search_range_tuple[1]), Float64(search_range_tuple[2])) for _ in 1:num_inputs]
            else
                 error("Length of SearchRange ($(length(search_range_tuple))) does not match num_inputs ($num_inputs). Provide a list of tuples, one for each input dimension.")
             end
        else
             # Assuming search_range_tuple is already a list of [low, high] pairs or tuples
             search_range = [(Float64(t[1]), Float64(t[2])) for t in search_range_tuple]
        end


        # --- Criterion Calculator ---
        criterion_calculator = default_end_period_criterion
        @info "Using criterion calculator: default_end_period_criterion. Minimizing: $minimize_criterion."

        # --- Output File Path ---
        # Define a default path, potentially make it configurable later
        default_output_file = "best_blackbox_input.json"
        output_file_path = get(cso_params, "output_file", default_output_file)
        @info "Best input will be saved to: $output_file_path"


        # --- Initialize State ---
        # Start with corner 0 as the initial "best" guess until evaluated
        initial_input = generate_corner_candidate(search_range, num_inputs, 0) # Corner index 0
        initial_score = minimize_criterion ? Inf : -Inf

        new(search_range, T_eval, num_inputs, num_corners, criterion_calculator, minimize_criterion, output_file_path,
            copy(initial_input), initial_score,  # best_input, best_score
            0,         # current_corner_index (start with index 0)
            nothing,   # current_candidate_input (generated on first compute call)
            -Inf,      # current_eval_start_time
            false,     # is_evaluating
            false      # cycle_complete
           )
    end
end

# --- Helper to generate a candidate based on corner index ---
function generate_corner_candidate(search_range::Vector{Tuple{Float64, Float64}}, num_inputs::Int, corner_index::Int)
    candidate = zeros(num_inputs)
    for i in 1:num_inputs
        low, high = search_range[i]
        # Check the (i-1)-th bit of corner_index (0 = low, 1 = high)
        # Right shift corner_index by (i-1) positions, then check the last bit using & 1
        if ((corner_index >> (i - 1)) & 1) == 0
            candidate[i] = low
        else
            candidate[i] = high
        end
    end
    return candidate
end

# --- Implement the Controller Interface ---

function compute_control_action(controller::BlackBoxOnlineController, time::Float64, system_state::Dict{Symbol, Any})

    # If we are currently evaluating a specific corner, keep applying its input
    if controller.is_evaluating
        @debug "Continuing evaluation for corner index $(controller.current_corner_index) at time $time."
        # Ensure current_candidate_input is not nothing if is_evaluating is true
        if controller.current_candidate_input === nothing
            @error "State inconsistency: is_evaluating is true, but current_candidate_input is nothing. Generating candidate for index $(controller.current_corner_index)."
            controller.current_candidate_input = generate_corner_candidate(controller.search_range, controller.num_inputs, controller.current_corner_index)
            controller.current_eval_start_time = time # Reset start time as we recovered state
        end
        return controller.current_candidate_input

    # If not evaluating, check if we have completed a full cycle and should stop/repeat
    elseif controller.cycle_complete
        @info "All $controller.num_corners corners tested. Applying best found input: $(round.(controller.best_input_so_far, digits=3))"
        # Keep applying the overall best input found after the cycle finished
        return controller.best_input_so_far

    # If not evaluating and cycle not complete, generate and start evaluating the next corner
    else
        # Generate the candidate for the *current* corner_index
        next_candidate = generate_corner_candidate(controller.search_range, controller.num_inputs, controller.current_corner_index)

        controller.current_candidate_input = next_candidate
        controller.current_eval_start_time = time
        controller.is_evaluating = true

        @info "Starting evaluation for Corner Index: $(controller.current_corner_index)/$(controller.num_corners-1), Input: $(round.(next_candidate, digits=3))"
        return next_candidate
    end
end


function update_controller_state!(controller::BlackBoxOnlineController, time::Float64, new_state_data::Dict{Symbol, Any})
    # Check if we were evaluating and if the evaluation period has ended
    if controller.is_evaluating && time >= controller.current_eval_start_time + controller.T_eval

        evaluated_candidate = controller.current_candidate_input
        evaluated_index = controller.current_corner_index
        @info "Evaluation period ended for Corner Index $evaluated_index: $(round.(evaluated_candidate, digits=3)) at time $time."

        # 1. Calculate the criterion for the corner that just finished
        local current_score::Float64
        try
            current_score = controller.criterion_calculator(new_state_data)
        catch e
            @error "Error calculating criterion for corner index $evaluated_index." exception=(e, catch_backtrace())
            current_score = controller.minimize_criterion ? Inf : -Inf
        end

        if !isfinite(current_score)
            @warn "Criterion calculation resulted in non-finite value ($current_score). Discarding result."
            current_score = controller.minimize_criterion ? Inf : -Inf
        end

        # 2. Compare with the best score found so far
        is_better = false
        if controller.minimize_criterion
            is_better = current_score < controller.best_score_so_far
        else # Maximize
            is_better = current_score > controller.best_score_so_far
        end

        if is_better
            @info "New best corner found! Index: $evaluated_index, Score: $current_score (Previous best: $(controller.best_score_so_far))"
            controller.best_score_so_far = current_score
            # Make sure to store a copy
            controller.best_input_so_far = copy(evaluated_candidate)

            # --- Save the new best input to a file ---
            try
                # Create a dictionary to store the best input along with its score for context
                output_data = Dict(
                    "best_input" => controller.best_input_so_far,
                    "best_score" => controller.best_score_so_far,
                    "time_found" => time,
                    "corner_index" => evaluated_index
                )
                # Open the file in write mode (overwrites previous best)
                open(controller.output_file_path, "w") do f
                    JSON.print(f, output_data, 4) # Use indentation for readability
                end
                @info "Saved new best input to $(controller.output_file_path)"
            catch e
                @error "Failed to save best input to $(controller.output_file_path)." exception=(e, catch_backtrace())
            end
            # --- End of saving logic ---

        else
            @info "Corner $evaluated_index score ($current_score) did not improve over best ($(controller.best_score_so_far))."
        end

        # 3. Advance to the next corner index
        controller.current_corner_index += 1

        # 4. Check if we have completed a full cycle
        if controller.current_corner_index >= controller.num_corners
            @info "Completed testing all $controller.num_corners corners."
            controller.cycle_complete = true
            controller.current_corner_index = 0 # Reset index if we wanted to cycle again later
            # If we only want to run once, cycle_complete flag handles stopping new evaluations
        end

        # 5. Reset evaluation state
        controller.current_candidate_input = nothing
        controller.current_eval_start_time = -Inf
        controller.is_evaluating = false # Ready for the next compute_control_action call

    elseif controller.is_evaluating
        @debug "Received state update during evaluation period for corner $controller.current_corner_index at time $time."
    end

    return nothing
end
