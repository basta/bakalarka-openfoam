# backend/src/server.jl
using HTTP
using JSON3
using Makie
using CairoMakie # For rendering plots to PNG
using Base64 # For encoding the image

# Define a function to handle WebSocket connections
function websocket_handler(ws::HTTP.WebSockets.WebSocket) # Explicitly type ws
    println("WebSocket client connected: ", ws)

    # Send a welcome message to the client
    try
        welcome_message = Dict(
            "type" => "STATUS_MESSAGE",
            "payload" => Dict(
                "level" => "INFO",
                "message" => "Welcome! Connected to Makie Plot Editor backend."
            )
        )
        HTTP.WebSockets.send(ws, JSON3.write(welcome_message))
    catch e
        println("Error sending welcome message: ", e)
    end

    while !eof(ws) # Check if the WebSocket is still open before trying to read
        try
            msg_bytes = HTTP.WebSockets.receive(ws) # Reads a single message
            if isempty(msg_bytes)
                continue
            end

            msg_string = String(msg_bytes)
            println("Received message from client: ", msg_string)

            json_msg = JSON3.read(msg_string)
            println("Parsed JSON message: ", json_msg)

            message_type = get(json_msg, :type, "UNKNOWN_TYPE")

            if message_type == "REQUEST_INITIAL_PLOT"
                println("Processing REQUEST_INITIAL_PLOT...")
                try
                    Makie.set_theme!(Makie.theme_light())
                    fig = Figure(resolution = (600, 400))
                    ax = Axis(fig[1, 1], xlabel = "X Label", ylabel = "Y Label", title = "Initial Plot")
                    xs = 0:0.1:10
                    ys = sin.(xs)
                    lines!(ax, xs, ys, color = :blue, linewidth = 2)

                    io_buffer = IOBuffer()
                    Makie.save(io_buffer, fig, format="png")
                    img_data_raw = take!(io_buffer)

                    img_base64 = base64encode(img_data_raw)
                    plot_data_uri = "data:image/png;base64," * img_base64

                    initial_plot_state_msg = Dict(
                        "type" => "INITIAL_PLOT_STATE",
                        "payload" => Dict(
                            "plotData" => plot_data_uri,
                            "elements" => [
                                Dict(
                                    "id" => "figure_main", "type" => "Figure", "name" => "Main Figure",
                                    "attributes" => [
                                        Dict("name" => "backgroundcolor", "value" => "white", "type" => "color", "label" => "Background"),
                                        Dict("name" => "resolution", "value" => [600,400], "type" => "dimensions_wh", "label" => "Resolution")
                                    ],
                                    "children" => ["axis_1"]
                                ),
                                Dict(
                                    "id" => "axis_1", "type" => "Axis", "name" => "Primary Axis",
                                    "attributes" => [
                                        Dict("name" => "title", "value" => "Initial Plot", "type" => "string", "label" => "Title"),
                                        Dict("name" => "xlabel", "value" => "X Label", "type" => "string", "label" => "X-Axis Label"),
                                        Dict("name" => "ylabel", "value" => "Y Label", "type" => "string", "label" => "Y-Axis Label")
                                    ],
                                    "children" => ["line_plot_1"]
                                ),
                                Dict(
                                    "id" => "line_plot_1", "type" => "Lines", "name" => "Sine Wave",
                                    "attributes" => [
                                        Dict("name" => "color", "value" => "blue", "type" => "color", "label" => "Line Color"),
                                        Dict("name" => "linewidth", "value" => 2, "type" => "number", "label" => "Line Width")
                                    ]
                                )
                            ],
                            "availableThemes" => ["theme_light", "theme_dark"],
                            "availableColormaps" => ["viridis", "plasma"]
                        )
                    )
                    HTTP.WebSockets.send(ws, JSON3.write(initial_plot_state_msg))
                    println("Sent INITIAL_PLOT_STATE to client.")

                catch plot_err
                    println("Error generating or sending initial plot: ", plot_err)
                    showerror(stdout, plot_err, catch_backtrace())
                    error_response = Dict(
                        "type" => "ERROR_MESSAGE",
                        "payload" => Dict("message" => "Failed to generate initial plot.", "details" => sprint(showerror, plot_err))
                    )
                    try HTTP.WebSockets.send(ws, JSON3.write(error_response)) catch end
                end
            elseif message_type != "UNKNOWN_TYPE"
                 response_payload = Dict(
                    "type" => "STATUS_MESSAGE",
                    "payload" => Dict("level" => "INFO", "message" => "Received your message of type: $(get(json_msg, :type, "UNKNOWN"))")
                )
                HTTP.WebSockets.send(ws, JSON3.write(response_payload))
            else
                println("Received message with unknown type or no type field.")
                 error_response = Dict(
                    "type" => "ERROR_MESSAGE",
                    "payload" => Dict("message" => "Received message with unknown type: $message_type", "details" => "Ensure your message has a 'type' field.")
                )
                try HTTP.WebSockets.send(ws, JSON3.write(error_response)) catch end
            end

        catch e
            if e isa HTTP.WebSockets.WebSocketError || e isa EOFError
                println("WebSocket read error or closed by client: ", e)
                break
            elseif e isa JSON3.JSONParseException
                println("JSON parsing error: ", e)
                error_response = Dict(
                    "type" => "ERROR_MESSAGE",
                    "payload" => Dict("message" => "Could not parse your message. Ensure it is valid JSON.", "details" => sprint(showerror, e))
                )
                try HTTP.WebSockets.send(ws, JSON3.write(error_response)) catch end
            else
                println("Error processing message or sending response: ", e)
                # Ensure ws is captured by the closure if this is inside an async task later
                current_ws = ws
                error_response = Dict(
                    "type" => "ERROR_MESSAGE",
                    "payload" => Dict("message" => "An error occurred while processing your request.", "details" => sprint(showerror, e))
                )
                try HTTP.WebSockets.send(current_ws, JSON3.write(error_response)) catch end
            end
        end
    end
    println("WebSocket client disconnected: ", ws)
end

# Define a dedicated handler for the /ws route
function ws_route_handler(req::HTTP.Request)
    if HTTP.WebSockets.is_upgrade(req)
        # The `HTTP.WebSockets.upgrade(req) do client_ws ... end` syntax
        # directly takes the HTTP.Request and handles the upgrade.
        # The `client_ws` is the resulting WebSocket.
        HTTP.WebSockets.upgrade(req; binary=true) do client_ws
            try
                websocket_handler(client_ws)
            catch ex
                println("Error in websocket_handler during WebSocket session: $ex")
                # Ensure WebSocket is closed if an error escapes the handler
                if isopen(client_ws)
                    try close(client_ws) catch close_ex println("Error closing WebSocket: $close_ex") end
                end
            end
        end
        # HTTP.WebSockets.upgrade handles returning the 101 Switching Protocols response
        # for a successful upgrade, or an error response if the upgrade fails.
    else
        # If the request to /ws is not a WebSocket upgrade request
        return HTTP.Response(400, "This endpoint is for WebSocket connections only.")
    end
end

function run_makie_server(host="127.0.0.1", port=8080)
    println("Starting Makie Plot Editor WebSocket server on ws://$host:$port/ws")
    router = HTTP.Router()

    # Register the dedicated WebSocket route handler
    HTTP.register!(router, "GET", "/ws", ws_route_handler)

    HTTP.register!(router, "GET", "/hello", req -> HTTP.Response(200, "Hello from Makie Server!"))

    try
        println("Server listening on $host:$port...")
        HTTP.serve(router, host, port; verbose=false)
    catch e
        if e isa Base.IOError && occursin("EADDRINUSE", string(e))
            println("ERROR: Port $port is already in use. Please choose a different port or stop the existing service.")
        elseif e isa InterruptException
            println("\nServer interrupted by user. Exiting.")
        else
            println("ERROR starting server: $e")
            showerror(stdout, e, catch_backtrace())
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Running server.jl directly for testing (prefer using run_server.jl)...")
    try
        using Pkg
        project_dir = dirname(dirname(@__DIR__))
        Pkg.activate(project_dir)
        Pkg.instantiate()
    catch e
        println("Could not activate/instantiate local Pkg environment: $e")
    end
    run_makie_server()
end
