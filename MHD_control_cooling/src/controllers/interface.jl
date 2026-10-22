# src/controllers/interface.jl
# Defines the abstract type and interface functions for all controllers.

# --- Abstract Type Definition ---
abstract type AbstractController end

# --- Interface Functions ---

"""
    compute_control_action(controller::AbstractController, time::Float64, state::Dict{Symbol, Any}) -> Vector{Float64}

Computes the control action based on the current time and system state.

# Arguments
- `controller`: The specific controller instance.
- `time`: Current simulation time.
- `state`: A dictionary containing the latest known system state (e.g., :time, :temperature, :wallheatflux).

# Returns
- `Vector{Float64}`: The computed control input vector (e.g., 8 actuator values).
"""
function compute_control_action(controller::AbstractController, time::Float64, state::Dict{Symbol, Any})
    # Ensure concrete subtypes implement this method
    error("compute_control_action not implemented for controller type $(typeof(controller))")
end

"""
    update_controller_state!(controller::AbstractController, time::Float64, new_state_data::Dict{Symbol, Any})

Updates the internal state of the controller, if necessary, based on new system data.
Useful for controllers with internal states like observers, filters, or integral terms.

# Arguments
- `controller`: The specific controller instance.
- `time`: Current simulation time corresponding to the `new_state_data`.
- `new_state_data`: A dictionary containing the latest received system state data.

# Returns
- `nothing`
"""
function update_controller_state!(controller::AbstractController, time::Float64, new_state_data::Dict{Symbol, Any})
    # Default implementation does nothing. Controllers requiring state updates should override this.
    return nothing
end

# --- Helper/Placeholder Functions (Potentially move to a dedicated utils file) ---

"""
Placeholder function for state estimation.
This needs to be implemented based on how you process raw temperature/flux data
into a state vector suitable for your controller model.
"""
function estimate_state(state_data::Dict{Symbol, Any})
    # Highly simplified example: Use average temperature and total flux
    # Replace this with your actual state estimation logic (e.g., using specific sensor points, POD modes, etc.)
    temp_data = get(state_data, :temperature, Float32[])
    flux_data = get(state_data, :wallheatflux, Float32[])

    avg_temp = isempty(temp_data) ? 0.0f0 : mean(temp_data)
    total_flux = isempty(flux_data) ? 0.0f0 : sum(flux_data)

    # Return a state vector. The structure depends entirely on your MPC model needs.
    return [avg_temp, total_flux]
end

"""
Placeholder function to get the target state for the controller.
This could be constant, time-varying, or based on external inputs.
"""
function get_target_state(time::Float64)
    # Example: Constant target state [target_avg_temp, target_total_flux]
    # Replace with your actual target definition logic.
    return [500.0, -1000.0] # Example target values
end
