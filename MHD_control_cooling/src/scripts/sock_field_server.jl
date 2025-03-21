using Sockets, Pkg
Pkg.activate(".")

include("../MHD_control_cooling.jl")
include("../openfoam/field_utils.jl")

@info "Script launched in pwd = $(pwd())"

using .MHD_control_cooling

function calculate_force_fields(inputs, centers)
    force = MHD_control_cooling.real_force_generator(inputs)
    cells = centers
    forces = [force(cell) for cell in eachcol(cells)]
    return forces
end

using CSV

using Dates

function create_csv_file_with_timestamp()
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    filename = "force_data_$timestamp.csv"
    file = open(filename, "w")
    println(file, "Time,u1,u2,u3,u4,u5,u6,u7,u8")
    close(file)
    return filename
end

csv_filename = create_csv_file_with_timestamp()

centers = []
# Function to start a server that sends vector data on request
function start_data_server(port::Int=8080)
    global centers = read_field_vector(joinpath(MHD_control_cooling.CASE_DIR, "0/C"))
    # Create a server socket that listens on the specified port
    server = listen(port)
    println("Server listening on port $port...")
    
    try
        while true
            # Accept incoming connection
            client = accept(server)
            println("Client connected from")
            
            # Handle client connection in an asynchronous task
            handle_client(client)
        end
    # catch e
    #     println("Server error: $e")
    finally
        close(server)
        println("Server shut down")
    end
end

# Function to handle each client connection
function handle_client(client::TCPSocket)

    try
        # Read the client's request
        @info "Reading request"
        request = readline(client)
        println("Received request: $request")
        
        if occursin("REQUEST_ARRAY", request)
            # Extract the float value from the request string
            m = match(r"REQUEST_ARRAY_(\d+\.\d+)", request)
            if m !== nothing
                float_value = parse(Float64, m.captures[1])
                println("Extracted float value: $float_value")
                
                # Generate some sample data based on the extracted float value
                data = calculate_force_fields(identification_signal(float_value, 0), centers)
            else
                println("No valid float value found in the request")
                data = []
            end
            
            # First send the array size
            array_size = length(data)
            write(client, array_size)
            
            # Then send each vector in the array
            for vec in data
                # Send each component of the 3D vector
                for component in vec
                    write(client, Float32(component))
                end
            end
            
            println("Sent array with $array_size vectors")
        else
            # Handle unknown request
            write(client, "Unknown request")
        end
    finally
        close(client)
        println("Client disconnected")
    end
end

holding_time = 0
last_change = -holding_time
last_out = zeros(8)
last_outs = []
N_MA = 5
function identification_signal(t::Real, T_stabilize::Real)
    if t < last_change
        global last_change = 0
    end
    if t > last_change + holding_time
        global holding_time = rand()*40+10
        global last_change = t
        if t < T_stabilize
            global last_out = zeros(8)
        else
            global last_out = (rand(8).-0.5).*2
        end
    end
    push!(last_outs, last_out)
    if length(last_outs) > N_MA
        global last_outs = last_outs[2:end]
    end
    out = sum(last_outs)/N_MA

    open(csv_filename, "a") do file
        println(file, join([t; out...], ","))
    end

    return out 
end

CASE_DIR = "../2d-example-live/"

function main()

    start_data_server(8080)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end