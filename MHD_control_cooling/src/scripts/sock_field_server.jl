using Sockets, Pkg
Pkg.activate(".") # Activate the project environment

# Assuming these files exist and contain necessary definitions
include("../MHD_control_cooling.jl")
include("../openfoam/field_utils.jl")

@info "Script launched in pwd = $(pwd())"

# Use the module defined in MHD_control_cooling.jl
using .MHD_control_cooling
using CSV
using Dates
using Statistics # Added for potential use with received data (e.g., mean)

# --- Global variable for cell centers ---
# Initialize as an empty array, will be populated when the server starts
centers = Matrix{Float64}(undef, 3, 0)

# --- Function to calculate force fields (remains unchanged) ---
function calculate_force_fields(inputs, centers_matrix)
    # Ensure centers_matrix is correctly formatted if needed by real_force_generator
    force = MHD_control_cooling.real_force_generator(inputs)
    # Assuming real_force_generator works with columns of the matrix
    forces = [force(cell) for cell in eachcol(centers_matrix)]
    return forces # Returns Vector{Vector{Float64}}
end

# --- Function to create a timestamped CSV file (remains unchanged) ---
function create_csv_file_with_timestamp()
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "force_data_$timestamp.csv"
    # Create or open the file and write the header
    open(filename, "w") do file
        println(file, "Time,u1,u2,u3,u4,u5,u6,u7,u8") # Assuming 8 control inputs
    end
    @info "Created CSV log file: $filename"
    return filename
end

# Create the CSV file at the start
csv_filename = create_csv_file_with_timestamp()

# --- Function to start the server ---
function start_data_server(port::Int=8080)
    global centers # Ensure we modify the global variable
    try
        # Read cell centers once when the server starts
        centers = read_field_vector(joinpath(MHD_control_cooling.CASE_DIR, "0/C"))
        @info "Successfully read cell centers. Size: $(size(centers))"
    catch e
        @error "Failed to read cell centers from $(joinpath(MHD_control_cooling.CASE_DIR, "0/C"))" exception = (e, catch_backtrace())
        @error "Server cannot start without cell centers. Exiting."
        return # Exit if centers cannot be read
    end

    # Create a server socket
    server = listen(port)
    println("Server listening on port $port...")

    try
        while true
            # Accept incoming connection
            client = accept(server)
            println("Client connected from $(getpeername(client))")

            # Handle client connection asynchronously using Tasks
            # This allows handling multiple clients concurrently if needed,
            # although the OpenFOAM setup likely uses one persistent connection.
            @async begin
                handle_client(client)
            end
        end
    catch e
        # Handle potential server errors (e.g., port already in use)
        @error "Server error: $e"
    finally
        close(server)
        println("Server shut down")
    end
end

# --- Function to handle each client connection (Modified) ---
function handle_client(client::TCPSocket)
    client_ip, client_port = getpeername(client)
    @info "Handling client: $client_ip:$client_port"
    try
        # Loop to handle multiple requests from the same client
        while isopen(client)
            local request # Ensure request is local to the loop iteration
            try
                # Read the client's request line by line
                request = readline(client)
                if isempty(request)
                    # Client might have closed connection gracefully after sending data
                    @info "Received empty line, potentially client closed connection."
                    break # Exit loop if empty line received
                end
                println("[$client_ip:$client_port] Received request: $request")

            catch e
                if e isa EOFError
                    @info "[$client_ip:$client_port] Client closed connection (EOF)."
                else
                    @error "[$client_ip:$client_port] Error reading request: $e"
                end
                break # Exit loop on any read error
            end

            # --- Process Request ---
            if startswith(request, "REQUEST_FORCES_")
                # --- Handle Force Request ---
                m = match(r"REQUEST_FORCES_(\d+(\.?\d*))", request) # Match integer or float time
                if m !== nothing
                    time_value = parse(Float64, m.captures[1])
                    println("[$client_ip:$client_port] Extracted time value: $time_value")

                    # Generate control input signal
                    input = identification_signal(time_value, 0) # Assuming T_stabilize is 0
                    @info "[$client_ip:$client_port] Generated control input: $input"

                    # Calculate force fields
                    force_data_vectors = calculate_force_fields(input, centers) # Vector{Vector{Float64}}

                    # Flatten the data and convert to Float32
                    flat_force_data = Float32[]
                    for vec in force_data_vectors
                        append!(flat_force_data, Float32.(vec)) # Append components as Float32
                    end

                    # Get the total number of floats
                    num_floats = Int32(length(flat_force_data)) # Ensure Int32

                    try
                        # 1. Send the number of floats (as Int32)
                        write(client, num_floats)

                        # 2. Send the flattened Float32 data
                        write(client, flat_force_data)

                        println("[$client_ip:$client_port] Sent $num_floats floats (representing $(length(force_data_vectors)) vectors) for forces.")
                        flush(client) # Ensure data is sent immediately

                    catch e
                        @error "[$client_ip:$client_port] Error sending force data: $e"
                        break # Exit loop on write error
                    end
                else
                    @warn "[$client_ip:$client_port] No valid time value found in force request: $request"
                    # Decide how to respond - maybe send size 0? Or just close?
                    # Sending size 0 might be safer if the client expects a response.
                    try
                        write(client, Int32(0)) # Send size 0
                        flush(client)
                    catch e
                        @error "[$client_ip:$client_port] Error sending zero size for invalid request: $e"
                        break
                    end
                end

            elseif startswith(request, "SEND_TEMP_")
                # --- Handle Temperature Reception ---
                m = match(r"SEND_TEMP_(\d+(\.?\d*))", request) # Match integer or float time
                time_value_str = (m !== nothing) ? m.captures[1] : "unknown"
                @info "[$client_ip:$client_port] Receiving temperature data for time ≈ $time_value_str"

                try
                    # 1. Read the number of floats (as Int32)
                    num_temp_floats = read(client, Int32)
                    @info "[$client_ip:$client_port] Expecting $num_temp_floats temperature values."

                    if num_temp_floats < 0
                        @warn "[$client_ip:$client_port] Received negative size for temperature data ($num_temp_floats). Aborting read."
                        continue # Skip to next request read attempt
                    elseif num_temp_floats == 0
                        @info "[$client_ip:$client_port] Received zero size for temperature data. No data to read."
                        continue
                    end


                    # 2. Read the Float32 temperature data
                    temp_data = Vector{Float32}(undef, num_temp_floats)
                    # Use read! for potentially better performance reading into existing vector
                    read!(client, temp_data) # Reads Float32 directly

                    @info "[$client_ip:$client_port] Successfully received $num_temp_floats temperature values."
                    # TODO: Process the received temp_data (e.g., log, use in control)
                    # Example: Print min/max temp received
                    if !isempty(temp_data)
                        println("[$client_ip:$client_port] Received Temp Min: $(minimum(temp_data)), Max: $(maximum(temp_data))")
                    end

                catch e
                    if e isa EOFError
                        @info "[$client_ip:$client_port] Client closed connection (EOF) while receiving temperature data."
                    else
                        @error "[$client_ip:$client_port] Error receiving temperature data: $e"
                    end
                    break # Exit loop on read error
                end

                # --- ADDED BLOCK: Handle Wall Heat Flux Reception ---
            elseif startswith(request, "SEND_WALLHEATFLUX_")
                m = match(r"SEND_WALLHEATFLUX_(\d+(\.?\d*))", request) # Match integer or float time
                time_value_str = (m !== nothing) ? m.captures[1] : "unknown"
                @info "[$client_ip:$client_port] Receiving wall heat flux data for time ≈ $time_value_str"

                try
                    # 1. Read the number of floats (as Int32) - this is the number of patches
                    num_flux_values = read(client, Int32)
                    @info "[$client_ip:$client_port] Expecting $num_flux_values wall heat flux value(s)."

                    if num_flux_values < 0
                        @warn "[$client_ip:$client_port] Received negative size for wall heat flux data ($num_flux_values). Aborting read."
                        continue # Skip to next request read attempt
                    elseif num_flux_values == 0
                        @info "[$client_ip:$client_port] Received zero size for wall heat flux data. No data to read."
                        continue
                    end

                    # 2. Read the Float32 heat flux data
                    flux_data = Vector{Float32}(undef, num_flux_values)
                    read!(client, flux_data) # Reads Float32 directly

                    @info "[$client_ip:$client_port] Successfully received $num_flux_values wall heat flux value(s)."
                    # TODO: Process the received flux_data
                    # Example: Print the received values
                    println("[$client_ip:$client_port] Received Wall Heat Fluxes (W/m^2): $flux_data")

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
                @warn "[$client_ip:$client_port] Received unknown request: $request"
                # No response needed as the client doesn't expect one for unknown requests
            end
        end # End while isopen(client)

    catch e
        # Catch errors occurring outside the inner read/write try blocks
        # but within the client handling scope.
        @error "Error handling client $client_ip:$client_port: $e"
        # Print stack trace for debugging
        showerror(stderr, e, catch_backtrace())
    finally
        # Ensure the client socket is closed when the loop exits or an error occurs
        if isopen(client)
            close(client)
        end
        @info "Client $client_ip:$client_port disconnected."
    end
end

# --- Function for generating identification signal (remains unchanged) ---
holding_time = 0.0 # Initialize holding time
last_change = -Inf # Initialize last change time
last_out = zeros(8) # Initialize last output
last_outs = [] # Initialize buffer for moving average
N_MA = 5 # Moving average window size
T_stabilize = 100.0 # Example stabilization time

function identification_signal(t::Real, T_stabilize_local::Real) # Renamed T_stabilize to avoid conflict
    global holding_time, last_change, last_out, last_outs # Ensure modification of globals

    # Reset last_change if time goes backwards (e.g., simulation restart)
    if t < last_change
        @warn "Time $t is less than last_change $last_change. Resetting last_change."
        last_change = t
    end

    # Check if holding time has elapsed
    if t >= last_change + holding_time
        holding_time = 150.0 # Set next holding duration
        last_change = t      # Record time of change

        # Generate new random output after stabilization period
        if t < T_stabilize # Use the global T_stabilize here
            last_out = zeros(8)
            @info "Time $t < T_stabilize $T_stabilize. Outputting zeros."
        else
            last_out = (rand(8) .- 0.5) .* 3.0 # Generate random values in [-1.5, 1.5]
            @info "Time $t >= T_stabilize $T_stabilize. Generating new random output: $last_out"
        end
    end

    # --- Moving Average (Optional - kept from original) ---
    # push!(last_outs, last_out)
    # if length(last_outs) > N_MA
    #     # Correct way to remove the first element
    #     popfirst!(last_outs)
    # end
    # Calculate moving average if needed - currently not used, just stores history
    # avg_out = mean(last_outs) # Requires Statistics pkg if you want to use mean
    # For now, just return the last generated output
    out = last_out
    # --- End Moving Average ---


    # Log the output to CSV
    try
        open(csv_filename, "a") do file
            # Format output for CSV: time, u1, u2, ...
            println(file, join([t; round.(out, digits=5)...], ","))
        end
    catch e
        @error "Failed to write to CSV file $csv_filename: $e"
    end

    return out # Return the last generated (non-averaged) output
end

# --- Define CASE_DIR (remains unchanged) ---
# Ensure this path is correct relative to where the script is run
CASE_DIR = "../2d-example/"

# --- Main function to start the server ---
function main()
    start_data_server(8080) # Start server on default port 8080
end

# --- Entry point: Run main() if script is executed directly ---
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
