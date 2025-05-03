# verify_mpc_linear.jl
using LinearAlgebra
using JLD2
using TOML
using Plots
using Logging

# Adjust the path if your controller file is located elsewhere
include("../controllers/mpc_dense.jl")

# --- Configuration ---
# Create a configuration dictionary similar to your TOML file
# **** MUST BE ADJUSTED BY USER ****
reference = -1e6
config = Dict(
    "controller" => Dict(
        "type" => "DenseMPC",
        "params" => Dict(
            "DenseMPC" => Dict(
                # --- CRITICAL: Set the path to your matrix file ---
                "matrix_file" => "./ABC.jld2",

                # --- MPC Tuning Parameters (Example Values - Adjust!) ---
                "Np" => 50,  # Prediction horizon
                "Nc" => 10,   # Control horizon
                "Q" => 1,  # Base weight for state error (applied to diagonal)
                             # Q[1,1] might be implicitly higher if xref focuses on x[1]
                "R_diag" => fill(0.1, 16), # Weight for control input theta (16 elements)
                "x0ref" => reference, # Base value for reference state (see below)
                "input_delays" => 100,

                # --- Constraints on theta ---
                "umin" => fill(-1.5, 16), # Lower bound for theta
                "umax" => fill(1.5, 16),  # Upper bound for theta

                # --- Optional OSQP Settings ---
                "osqp_settings" => Dict(
                    # :verbose => true, # Uncomment for solver details
                    :eps_abs => 1e-4,
                    :eps_rel => 1e-4,
                    :max_iter => 10000
                )
                # Add other necessary params like "input_delays" if your constructor uses them
                # "input_delays" => 2, # Example if needed by constructor/state logic
            )
        )
    )
)
# --- End Configuration ---

# --- Simulation Parameters ---
Tsim = 5000.       # Simulation time
dt = 5.       # Time step (SHOULD MATCH model discretization dt)
Nsim = Int(round(Tsim / dt))

# --- Load Model and Instantiate Controller ---
local A, B, mean_vec, controller
try
    # Instantiate the controller using the config
    global controller = DenseMPCController(config)

    # Extract A, B, mean_vec (for simulation state updates and reference)
    # We assume the controller loaded these correctly internally
    # Access them via the controller's mpc_data field
    global A = controller.mpc_data.A
    global B = controller.mpc_data.B
    global mean_vec = controller.mean_vec # Needed for simulating measurement
    global nx = controller.mpc_data.nx
    global nu = controller.mpc_data.nu # Should be 16

    @info "Controller instantiated successfully."
    @info "Model dimensions: A is $(size(A)), B is $(size(B)), nx=$nx, nu=$nu"
catch e
    @error "Failed to load model or instantiate controller:" exception=(e, catch_backtrace())
    exit()
end

# --- Reference Trajectory ---
# Define the target state trajectory *in centered coordinates*
xref_sim = zeros(nx, Nsim + 1)
target_y_value = reference
step_start_time = 1.0
step_start_idx = Int(round(step_start_time / dt)) + 1

# Set the reference for the first state element after the step time
# Keep the rest of the reference state at zero (assuming target is steady state at origin)
# A more sophisticated reference might be needed if x0ref != 0 or target theta != 0
xref_sim[1, step_start_idx:end] .= target_y_value

# --- Initial State ---
# Start at the origin (zero deviation from mean_vec)
x0 = zeros(nx)
xk = copy(x0) # Current state

# --- Data Logging ---
time_hist = zeros(Nsim + 1)
x_hist = zeros(nx, Nsim + 1)
theta_hist = zeros(nu, Nsim) # Store calculated theta
x_hist[:, 1] = xk

# --- Simulation Loop ---
@info "Starting linear simulation..."
for k = 1:Nsim
    current_time = (k - 1) * dt
    time_hist[k] = current_time

    # 1. Prepare system_state for the controller
    # Controller expects the "raw" measurement, which it centers internally.
    # Here, the raw measurement is xk[1] (centered) + mean_vec[1] (offset)
    # We assume the primary measured output corresponds to xk[1].
    measurement = xk[1]
    if length(mean_vec) >= 1
         measurement += mean_vec[1]
    end
    # The controller expects a Dict, potentially with specific keys
    # We mimic providing the 'wallheatflux' based on the original controller code
    system_state = Dict{Symbol, Any}(:wallheatflux => [measurement])

    # 2. Update controller reference if necessary (optional, depends on strategy)
    # For this test, the controller uses the constant xref from its config.
    # If you wanted a time-varying reference *for the MPC internal target*,
    # you would update controller.mpc_data.xref here.
    # We will compare x_hist to xref_sim later.

    # 3. Call the controller's main function
    # We ignore the return value (8D physical input), as we need theta_k
    try
        compute_control_action(controller, current_time, system_state)
    catch e
         @error "Error during compute_control_action at step $k:" exception=(e, catch_backtrace())
         # Decide how to handle: break, continue with zero input?
         break
    end


    # 4. Get the calculated internal control theta_k
    # The controller stores the result in 'last_uk' internally
    theta_k = controller.last_uk
    if theta_k === nothing
        @warn "Controller did not produce a control input at step $k. Using zeros."
        theta_k = zeros(nu)
    end
    theta_hist[:, k] = theta_k

    # 5. Simulate the next state using the LINEAR MODEL
    try
        global xk = A * xk + B * theta_k
    catch e
        @error "Error during state update (A*xk + B*theta_k) at step $k:" exception=(e, catch_backtrace())
        break
    end

    # 6. Store history
    x_hist[:, k+1] = xk
    time_hist[k+1] = current_time + dt

    # Optional: Logging progress
    if k % 10 == 0
        @info "Sim step $k/$Nsim, x[1]=$(round(xk[1], digits=4)), theta[1]=$(round(theta_k[1], digits=4))"
    end
end
@info "Simulation finished."

# --- Plotting Results ---
@info "Plotting results..."
# Plot the first state variable (output) vs. reference
p1 = plot(time_hist, x_hist[1, :], label="x[1] (Simulated Output)", xlabel="Time (s)", ylabel="State Value", title="MPC Tracking Verification (Linear Model)")
plot!(p1, time_hist, xref_sim[1, :], label="x[1] Reference", linestyle=:dash)

# Plot the first element of the control input theta
p2 = plot(time_hist[1:Nsim], theta_hist[1, :], label="theta[1] (Control Input)", xlabel="Time (s)", ylabel="Control Value", title="Internal Control Input")
# You could plot more elements of theta if desired

# Combine plots
plot(p1, p2, layout=(2, 1), legend=true)

# Save the plot (optional)
savefig("mpc_linear_verification.png")
@info "Plot displayed. Close plot window to exit."

# Keep plot window open until closed
gui() # May be needed in some environments
