# src/scripts/sock_field_server.jl
# Main server script for OpenFOAM control experiments

# --- Core Julia Libraries ---
using Sockets
using Dates
using Statistics
using TOML # For reading config files (add Pkg.add("TOML") if needed)
using ArgParse # For command-line arguments (add Pkg.add("ArgParse") if needed)

# --- Project-Specific Includes ---
# Note: Adjust paths if this script is moved or run from a different directory
try
    # Assuming MHD_control_cooling.jl defines a module with necessary functions
    include("../MHD_control_cooling.jl")
    using .MHD_control_cooling # Or specific function imports

    # OpenFOAM utilities
    include("../openfoam/field_utils.jl")

    # Controller definitions
    include("../controllers/interface.jl") # Defines AbstractController etc.
    include("../controllers/identification.jl") # Defines IdentificationController
    include("../controllers/mpc.jl") # Defines MPCController
    include("../controllers/mpc_dmdc.jl") # Defines DMDC_MPCController

catch e
    @error "Failed to include necessary project files. Ensure paths are correct relative to script location." exception=(e, catch_backtrace())
    exit(1)
end


# --- Global Variables ---
global centers = Matrix{Float64}(undef, 3, 0) # Populated from OpenFOAM case
global current_controller::Union{AbstractController, Nothing} = nothing # Active controller instance
global system_state = Dict{Symbol, Any}( # Dictionary to hold the latest known state
    :time => 0.0,
    :temperature => Float32[],
    :wallheatflux => Float32[],
    # Add other state variables if needed (e.g., velocity, pressure)
)
global config::Dict = Dict() # Loaded experiment configuration
global csv_logger::Union{IOStream, Nothing} = nothing # File handle for primary CSV logging
global state_logger::Union{IOStream, Nothing} = nothing # Optional: File handle for detailed state logging
global last_state_log_time::Float64 = -Inf # Track when state was last logged

# --- Configuration & Setup Functions ---

function parse_commandline()
    s = ArgParseSettings(description="Run the OpenFOAM control server with a specified configuration.")

    @add_arg_table! s begin
        "--config", "-c"
            help = "Path to the experiment configuration TOML file."
            arg_type = String
            default = "configs/default.toml" # Default config path
        # Add other command-line args if needed (e.g., overriding port)
    end

    return parse_args(s)
end

function load_configuration(filepath::String)
    global config
    try
        config = TOML.parsefile(filepath)
        @info "Configuration loaded successfully from '$filepath'"
        # Basic validation (add more specific checks as needed)
        if !haskey(config, "controller") || !haskey(config["controller"], "type")
            error("Configuration missing [controller] section or 'type' key.")
        end
        if !haskey(config, "simulation") || !haskey(config["simulation"], "case_dir")
             @warn "Configuration missing [simulation] section or 'case_dir'. Using default '.'."
             config["simulation"]["case_dir"] = "." # Provide a default if missing
        end
         if !haskey(config, "logging") || !haskey(config["logging"], "log_dir")
             @warn "Configuration missing [logging] section or 'log_dir'. Using default 'experiment_runs'."
             config["logging"]["log_dir"] = "experiment_runs" # Provide a default
        end
    catch e
        @error "Failed to load or parse configuration file '$filepath'." exception=(e, catch_backtrace())
        rethrow(e) # Propagate error to stop server start
    end
end

function setup_logging(config_path::String)
    global csv_logger, state_logger, config
    log_dir_base = get(config["logging"], "log_dir", "experiment_runs")
    base_filename = get(config["logging"], "base_filename", "run_data")
    controller_type = get(config["controller"], "type", "unknown_controller")
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    experiment_dir = joinpath(log_dir_base, "$(timestamp)_$(controller_type)_$(base_filename)")

    try
        mkpath(experiment_dir)
        @info "Created experiment log directory: $experiment_dir"

        # Copy config file for reproducibility
        try
            cp(config_path, joinpath(experiment_dir, "config.toml"), force=true)
        catch e
             @warn "Could not copy config file '$config_path' to log directory." exception=e
        end


        # --- Setup Primary CSV Logger (Control Actions) ---
        csv_filename = joinpath(experiment_dir, "$(base_filename)_control.csv")
        csv_logger = open(csv_filename, "w")
        # Define header based on expected logged data
        # Assuming 8 control inputs (u1-u8)
        header_control = "Time," * join(["u$i" for i in 1:8], ",")
        println(csv_logger, header_control)
        flush(csv_logger)
        @info "Opened CSV log for control inputs: $csv_filename"

        # --- Optional: Setup State CSV Logger ---
        log_state_detail = get(config["logging"], "log_state_detail", false) # Check if detailed state logging is enabled
        if log_state_detail
             state_csv_filename = joinpath(experiment_dir, "$(base_filename)_state.csv")
             state_logger = open(state_csv_filename, "w")
             # Define header for state log (example)
             header_state = "Time,AvgTemp,MinTemp,MaxTemp,TotalFlux,HotFlux,ColdFlux" # Customize as needed
             println(state_logger, header_state)
             flush(state_logger)
             @info "Opened CSV log for detailed state: $state_csv_filename"
        end

    catch e
        @error "Failed to set up logging directory or file(s)." exception=(e, catch_backtrace())
        rethrow(e) # Propagate error
    end
end


function initialize_controller()
    global current_controller, config
    controller_type = config["controller"]["type"]
    params = get(config, "controller", Dict()) # Get controller section
    specific_params = get(params, "params", Dict()) # Get nested params section
    controller_specific_params = get(specific_params, controller_type, Dict()) # Params for this specific controller type

    @info "Initializing controller of type: $controller_type"

    try
        if controller_type == "Identification"
            params = get(config, "controller", Dict()) # Get controller section for common params
            T_stabilize = get(params, "T_stabilize", 100.0) # Common param might be outside nested params
            current_controller = IdentificationController(T_stabilize)

            elseif controller_type == "DMDC_MPC"
                # The constructor handles extracting parameters from the config dictionary
                current_controller = DMDC_MPCController(config)
        # Add elseif blocks for other controller types (e.g., PID)
        # elseif controller_type == "PID"
        #     Kp = get(...)
        #     Ki = get(...)
        #     Kd = get(...)
        #     current_controller = PIDController(Kp, Ki, Kd)

        else
            error("Unsupported controller type specified in configuration: '$controller_type'")
        end
        @info "Controller initialized successfully: $(typeof(current_controller))"
    catch e
         @error "Failed to initialize controller '$controller_type'." exception=(e, catch_backtrace())
         rethrow(e) # Propagate error
    end
end

function read_cell_centers()
    global centers, config
    case_dir = get(config["simulation"], "case_dir", ".") # Get from config
    centers_path = joinpath(case_dir, "0/C") # Standard location for cell centers
    @info "Attempting to read cell centers from: $centers_path"
    try
        # read_field_vector should be defined in openfoam/field_utils.jl
        centers = read_field_vector(centers_path)
        if isempty(centers)
             error("Read cell centers but the resulting matrix is empty.")
        end
        @info "Successfully read cell centers. Size: $(size(centers))"
    catch e
        @error "Failed to read cell centers from '$centers_path'." exception=(e, catch_backtrace())
        @error "Server cannot start without cell centers."
        rethrow(e) # Propagate error
    end
end


# --- Core Server Logic ---

function log_data(time::Float64, input::Vector{Float64})
    global csv_logger, state_logger, system_state, config, last_state_log_time

    # --- Log Control Input (Always) ---
    if csv_logger !== nothing && isopen(csv_logger)
        try
            # Format: Time, u1, u2, ..., u8
            println(csv_logger, join([time; round.(input, digits=6)...], ","))
            # flush(csv_logger) # Flushing frequently can impact performance, consider removing if logging rate is high
        catch e
            @error "Failed to write control input to CSV log at time $time: $e"
        end
    end

    # --- Log System State (Periodically) ---
    log_interval = get(config["logging"], "log_state_interval", 1.0)
    log_state_detail = get(config["logging"], "log_state_detail", false)

    if log_state_detail && state_logger !== nothing && isopen(state_logger) && time >= last_state_log_time + log_interval
         try
            # Extract state metrics (handle empty data)
            temp_data = system_state[:temperature]
            flux_data = system_state[:wallheatflux]

            avg_temp = isempty(temp_data) ? NaN32 : mean(temp_data)
            min_temp = isempty(temp_data) ? NaN32 : minimum(temp_data)
            max_temp = isempty(temp_data) ? NaN32 : maximum(temp_data)
            total_flux = isempty(flux_data) ? NaN32 : sum(flux_data)
            hot_flux = isempty(flux_data) ? NaN32 : flux_data[1]
            cold_flux = isempty(flux_data) ? NaN32 : flux_data[2]

            # Format matches header defined in setup_logging
            println(state_logger, join([time, avg_temp, min_temp, max_temp, total_flux, hot_flux, cold_flux], ","))
            flush(state_logger)
            last_state_log_time = time
         catch e
             @error "Failed to write state data to log at time $time: $e"
         end
    end
end

function handle_client(client::TCPSocket)
    global system_state, current_controller, centers # Use globals

    client_ip, client_port = getpeername(client)
    @info "Handling client connection: $client_ip:$client_port"

    try
        while isopen(client)
            local request
            try
                # Read request line by line (assuming text-based commands from OpenFOAM)
                request = readline(client)
                if isempty(request)
                    @info "[$client_ip:$client_port] Received empty line, assuming client closed connection."
                    break # Exit loop for this client
                end
                # @debug "[$client_ip:$client_port] Received request: $request" # Use debug for less noise
            catch e
                if e isa EOFError
                    @info "[$client_ip:$client_port] Client closed connection (EOF)."
                else
                    @error "[$client_ip:$client_port] Error reading request: $e"
                end
                break # Exit loop on read error
            end

            # --- Process Request ---
            if startswith(request, "REQUEST_FORCES_")
                 m = match(r"REQUEST_FORCES_(\d+(\.?\d*([eE][-+]?\d+)?))", request) # Match float/scientific notation time
                 if m !== nothing
                     time_value = parse(Float64, m.captures[1])
                     system_state[:time] = time_value # Update time in global state

                     # *** Core Control Logic Call ***
                     if current_controller !== nothing
                         local input # Control input vector
                         try
                            # Compute control action using the active controller and current state
                            # input = compute_control_action(current_controller, time_value, system_state)
                            input = [-2.6134549598668113, 1.5088922853126288, 1.5088901387951934, 1.508845830268384, -2.61346857727252, -1.5088389626108187, -1.5088604541887394, 1.5089052510684984]
                            input = [-2.6134549598668113, 1.5088922853126288, 1.5088901387951934, 1.508845830268384, 2.61346857727252, 1.5088389626108187, 1.5088604541887394, -1.5089052510684984]
                         catch e
                             @error "[$client_ip:$client_port] Error during compute_control_action at time $time_value" exception=(e, catch_backtrace())
                             # Fallback: Send zero forces if controller fails
                             input = zeros(8) # Assuming 8 control inputs
                         end


                         # Log time and computed input (and potentially state)
                         log_data(time_value, input)

                         # Calculate forces based on controller output
                         local force_data_vectors
                         try
                             # calculate_force_fields should be defined elsewhere (e.g., MHD_control_cooling.jl)
                             force_data_vectors = calculate_force_fields(input, centers)
                         catch e
                              @error "[$client_ip:$client_port] Error during calculate_force_fields at time $time_value" exception=(e, catch_backtrace())
                              # Fallback: Send zero forces if calculation fails
                              force_data_vectors = [zeros(3) for _ in 1:size(centers, 2)] # Zero vector for each cell
                         end


                         # Flatten, convert to Float32, and send data back to OpenFOAM
                         flat_force_data = Float32[]
                         for vec in force_data_vectors
                             append!(flat_force_data, Float32.(vec)) # Append components
                         end
                         num_floats = Int32(length(flat_force_data)) # Must be Int32 for OpenFOAM

                         try
                             write(client, num_floats) # Send size first
                             write(client, flat_force_data) # Then send data
                             flush(client) # Ensure data is sent immediately
                             # @debug "[$client_ip:$client_port] Sent $num_floats force floats for time $time_value."
                         catch e
                             @error "[$client_ip:$client_port] Error sending force data: $e"
                             break # Exit loop on write error
                         end
                     else
                         @error "[$client_ip:$client_port] No controller initialized! Cannot compute forces."
                         # Send zero forces as fallback if controller is somehow null
                         try
                             write(client, Int32(0))
                             flush(client)
                         catch e; break; end
                     end
                 else
                     @warn "[$client_ip:$client_port] Could not parse time value from force request: $request"
                     # Send zero forces if request format is wrong
                     try
                         write(client, Int32(0))
                         flush(client)
                     catch e; break; end
                 end

            elseif startswith(request, "SEND_TEMP_")
                # --- Handle Temperature Reception ---
                m = match(r"SEND_TEMP_(\d+(\.?\d*([eE][-+]?\d+)?))", request)
                time_value_str = (m !== nothing) ? m.captures[1] : "unknown"
                # @debug "[$client_ip:$client_port] Receiving temperature data for time ≈ $time_value_str"
                try
                    # 1. Read the number of floats (Int32)
                    num_temp_floats = read(client, Int32)
                    # @debug "[$client_ip:$client_port] Expecting $num_temp_floats temperature values."

                    if num_temp_floats < 0
                        @warn "[$client_ip:$client_port] Received negative size ($num_temp_floats) for temperature data. Ignoring."
                        continue # Skip reading data for this request
                    elseif num_temp_floats == 0
                         # @debug "[$client_ip:$client_port] Received zero size for temperature data."
                         system_state[:temperature] = Float32[] # Clear previous data
                    else
                        # 2. Read the Float32 temperature data
                        temp_data = Vector{Float32}(undef, num_temp_floats)
                        read!(client, temp_data) # Reads Float32 directly into the vector
                        system_state[:temperature] = temp_data # *** Update global state ***
                        # @debug "[$client_ip:$client_port] Received Temp Min: $(minimum(temp_data)), Max: $(maximum(temp_data))"
                    end

                    # 3. Update controller state if needed (e.g., for observers)
                    if current_controller !== nothing
                        update_controller_state!(current_controller, system_state[:time], system_state)
                    end

                catch e
                    if e isa EOFError
                         @info "[$client_ip:$client_port] Client closed connection (EOF) while receiving temperature data."
                    else
                         @error "[$client_ip:$client_port] Error receiving temperature data: $e"
                    end
                    break # Exit loop on read error
                end

            elseif startswith(request, "SEND_WALLHEATFLUX_")
                 # --- Handle Wall Heat Flux Reception ---
                 m = match(r"SEND_WALLHEATFLUX_(\d+(\.?\d*([eE][-+]?\d+)?))", request)
                 time_value_str = (m !== nothing) ? m.captures[1] : "unknown"
                 # @debug "[$client_ip:$client_port] Receiving wall heat flux data for time ≈ $time_value_str"
                 try
                     # 1. Read the number of floats (Int32)
                     num_flux_values = read(client, Int32)
                     # @debug "[$client_ip:$client_port] Expecting $num_flux_values wall heat flux value(s)."

                     if num_flux_values < 0
                        @warn "[$client_ip:$client_port] Received negative size ($num_flux_values) for flux data. Ignoring."
                        continue # Skip reading data
                     elseif num_flux_values == 0
                        # @debug "[$client_ip:$client_port] Received zero size for wall heat flux data."
                        system_state[:wallheatflux] = Float32[] # Clear previous data
                     else
                        # 2. Read the Float32 heat flux data
                        flux_data = Vector{Float32}(undef, num_flux_values)
                        read!(client, flux_data) # Reads Float32 directly
                        system_state[:wallheatflux] = flux_data # *** Update global state ***
                        # @debug "[$client_ip:$client_port] Received Wall Heat Fluxes: $flux_data"
                     end

                     # 3. Update controller state if needed
                     if current_controller !== nothing
                        update_controller_state!(current_controller, system_state[:time], system_state)
                     end

                 catch e
                     if e isa EOFError
                         @info "[$client_ip:$client_port] Client closed connection (EOF) while receiving wall heat flux data."
                     else
                         @error "[$client_ip:$client_port] Error receiving wall heat flux data: $e"
                     end
                     break # Exit loop on read error
                 end

            else
                # --- Handle Unknown Request ---
                @warn "[$client_ip:$client_port] Received unknown or malformed request: $request"
                # Optionally send an error message back, or just ignore
            end
        end # End while isopen(client)

    catch e
        # Catch errors occurring outside the inner read/write try blocks
        # but within the client handling scope.
        @error "Unhandled error handling client $client_ip:$client_port: $e"
        showerror(stderr, e, catch_backtrace()) # Print stack trace for debugging
    finally
        # Ensure the client socket is closed when the loop exits or an error occurs
        if isopen(client)
            close(client)
        end
        @info "Client $client_ip:$client_port disconnected."
    end
end


function start_data_server(config_path::String)
    global config, csv_logger, state_logger # Use globals

    # --- Load Configuration ---
    try
        load_configuration(config_path)
    catch e
        @error "Exiting server due to configuration error."
        return # Stop execution
    end

    # --- Setup Logging ---
    try
        setup_logging(config_path)
    catch e
        @error "Exiting server due to logging setup error."
        # Ensure any partially opened log files are closed if possible
        if csv_logger !== nothing && isopen(csv_logger) close(csv_logger) end
        if state_logger !== nothing && isopen(state_logger) close(state_logger) end
        return # Stop execution
    end

    # --- Read Cell Centers ---
    try
        read_cell_centers()
    catch e
        @error "Exiting server because cell centers could not be read."
        if csv_logger !== nothing && isopen(csv_logger) close(csv_logger) end
        if state_logger !== nothing && isopen(state_logger) close(state_logger) end
        return # Stop execution
    end

    # --- Initialize Controller ---
    try
        initialize_controller()
    catch e
        @error "Exiting server due to controller initialization error."
        if csv_logger !== nothing && isopen(csv_logger) close(csv_logger) end
        if state_logger !== nothing && isopen(state_logger) close(state_logger) end
        return # Stop execution
    end

    # --- Start Server Listening ---
    port = get(config["simulation"], "port", 8080) # Get port from config, default 8080
    server_socket = nothing # Initialize variable

    try
        server_socket = listen(port)
        println("\n" * "="^40)
        println("Server listening on port $port...")
        println("Using controller: $(typeof(current_controller))")
        println("Logging to directory: $(dirname(csv_logger.name))") # Show log dir
        println("Configuration file: $config_path")
        println("="^40 * "\n")

        # Main server loop to accept client connections
        while true
            client_socket = accept(server_socket)
            println("Client connected from $(getpeername(client_socket))")
            # Handle each client connection asynchronously using Tasks
            # This allows handling potential multiple connections, though OpenFOAM usually uses one.
            @async begin
                handle_client(client_socket)
            end
        end
    catch e
        # Handle server errors (e.g., port already in use, network issues)
        if e isa Base.IOError && occursin("address already in use", e.msg)
             @error "Server error: Port $port is already in use."
        elseif e isa InterruptException
             println("\nServer interrupted by user (Ctrl+C). Shutting down...")
        else
             @error "Server error encountered:" exception=(e, catch_backtrace())
        end
    finally
        # --- Cleanup ---
        println("Shutting down server...")
        if server_socket !== nothing && isopen(server_socket)
            close(server_socket)
            println("Server socket closed.")
        end
        if csv_logger !== nothing && isopen(csv_logger)
            close(csv_logger)
            println("Control log file closed.")
        end
         if state_logger !== nothing && isopen(state_logger)
            close(state_logger)
            println("State log file closed.")
        end
        println("Server shut down complete.")
    end
end

# --- Main Execution Block ---
function main()
    parsed_args = parse_commandline()
    config_file = parsed_args["config"]

    if !isfile(config_file)
        @error "Configuration file not found at specified path: '$config_file'. Please provide a valid path using --config or -c."
        exit(1)
    end

    start_data_server(config_file)
end

# --- Script Entry Point ---
# Ensures main() is called only when the script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
