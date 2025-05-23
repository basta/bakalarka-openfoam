# src/controllers/identification.jl
# Implements the Identification Signal controller.

# Ensure the interface is available
# include("interface.jl") # Usually handled by the main script's include order or module structure

"""
    IdentificationController

A simple controller that generates a piecewise constant random signal,
often used for system identification experiments.
"""
mutable struct IdentificationController <: AbstractController
    holding_time::Float64       # Current duration to hold the output constant
    last_change_time::Float64   # Simulation time when the output last changed
    last_output::Vector{Float64} # The current constant output vector
    T_stabilize::Float64        # Initial period during which output is zero
    num_inputs::Int             # Number of control inputs

    # Constructor
    function IdentificationController(T_stabilize::Real=100.0, num_inputs::Int=8)
        new(0.0, -Inf, zeros(Float64, num_inputs), Float64(T_stabilize), num_inputs)
    end
end

# Implement the interface function
function compute_control_action(controller::IdentificationController, time::Float64, state::Dict{Symbol, Any})
    # Reset last_change_time if time goes backwards (e.g., simulation restart)
    if time < controller.last_change_time
        @warn "[$(typeof(controller))] Time $time is less than last_change_time $(controller.last_change_time). Resetting last_change_time."
        controller.last_change_time = time
    end

    # Check if holding time has elapsed
    if time >= controller.last_change_time + controller.holding_time
        # Determine next holding duration
        controller.holding_time = rand(50:200)
        controller.last_change_time = time

        # Generate new random output after stabilization period
        if time < controller.T_stabilize
            controller.last_output = zeros(controller.num_inputs)
            # @info "[$(typeof(controller))] Time $time < T_stabilize $(controller.T_stabilize). Outputting zeros."
        else
            # Generate random values, e.g., in [-1.5, 1.5]
            controller.last_output = (rand(controller.num_inputs) .- 0.) .* 1.
            # @info "[$(typeof(controller))] Time $time >= T_stabilize $(controller.T_stabilize). New random output: $(round.(controller.last_output, digits=3))"
        end
    end

    # Return the currently held output value
    return controller.last_output
end

# Identification controller is stateless regarding system observations,
# so the default update_controller_state! implementation (does nothing) is sufficient.
# function update_controller_state!(controller::IdentificationController, time::Float64, new_state_data::Dict{Symbol, Any})
#     return nothing
# end
