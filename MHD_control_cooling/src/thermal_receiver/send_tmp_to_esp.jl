# thermal_sender/gui_send_max_temp_to_esp32.jl

using Sockets
using Printf
using GLMakie # For the GUI
using Dates
using LibSerialPort # For ESP32 communication
using Statistics  # For maximum

# Configure Makie to run in a non-blocking way for interactive use
GLMakie.activate!(inline=false)


# These files are expected to be in the same directory or accessible via Julia's load path.
try
    include("config.jl")      # Defines ServerConfig with TCP server IP, port, frame dimensions etc.
    include("image_processor.jl") # Defines ImageProcessor.process_tcp_frame
catch e
    println("Error including config.jl or image_processor.jl: $e")
    println("Please ensure these files exist and are correctly formatted.")
    println("Using placeholder definitions for ServerConfig and ImageProcessor for script structure.")

    function process_tcp_frame(tcp_frame_bytes::Vector{UInt8}, config)::ProcessedImageData
        expected_bytes = config.FRAME_WIDTH * config.FRAME_HEIGHT * 2
        if length(tcp_frame_bytes) < expected_bytes
            @warn "Incomplete frame received. Expected $expected_bytes, got $(length(tcp_frame_bytes))"
            return ProcessedImageData(false, nothing)
        end
        try
            raw_pixel_data = reinterpret(UInt16, tcp_frame_bytes[1:expected_bytes])
            image_matrix = reshape(raw_pixel_data, (config.FRAME_HEIGHT, config.FRAME_WIDTH))
            return ProcessedImageData(true, image_matrix)
        catch e
            @error "Error processing TCP frame: $e"
            return ProcessedImageData(false, nothing)
        end
    end
end # End of try-catch for includes


# --- Temperature Conversion Functions ---
function raw_to_celsius(raw_value::UInt16)::Float32
    return Float32(Float64(raw_value) * 0.0984 - 265.82)
end

function raw_matrix_to_celsius_matrix(raw_matrix::Matrix{UInt16})::Matrix{Float32}
    return Float32.((Float64.(raw_matrix) .* 0.0984) .- 265.82)
end
# --- End Temperature Conversion Functions ---


# --- Serial Port Functions ---
function list_serial_ports_lib()
    ports = LibSerialPort.get_port_list()
    if isempty(ports)
        println("No serial ports found by LibSerialPort.")
        return nothing
    end
    println("Available serial ports (via LibSerialPort):")
    for (i, port) in enumerate(ports)
        println("[$i] $port")
    end
    return ports
end

function select_serial_port(ports::Union{Vector{String},Nothing})::Union{String,Nothing}
    if ports === nothing || isempty(ports)
        return nothing
    end

    print("Enter the number of the serial port for ESP32: ")
    idx_str = readline()
    idx = tryparse(Int, idx_str)

    if idx === nothing || !(1 <= idx <= length(ports))
        println("Invalid selection.")
        return nothing
    end
    return ports[idx]
end
# --- End Serial Port Functions ---


function main_gui_with_makie()
    println("Thermal Camera GUI and ESP32 Max Temp Sender (with Duty Cycle Plot)")

    available_ports = list_serial_ports_lib()
    esp32_port_name = select_serial_port(available_ports)

    if esp32_port_name === nothing
        println("No serial port selected for ESP32. Exiting.")
        return
    end
    esp32_baud_rate = 115200

    # --- Makie Figure and Layout Setup ---
    fig = Figure(size=(1200, 900))

    # Main layout: Heatmap | Charts (Temp + Duty)
    # Sliders below everything
    main_grid = fig[1, 1] = GridLayout()
    slider_grid = fig[2, 1] = GridLayout()

    # Left panel for heatmap
    left_panel = main_grid[1, 1] = GridLayout()
    ax_heatmap = Axis(left_panel[1, 1], title="Live Thermal Image (°C)", aspect=DataAspect())

    frame_height = try
        ServerConfig.FRAME_HEIGHT
    catch
        60
    end
    frame_width = try
        ServerConfig.FRAME_WIDTH
    catch
        80
    end
    obs_thermal_data = Observable(zeros(Float32, frame_height, frame_width))

    slider_min_val_initial = 15.0
    slider_max_val_initial = 40.0
    slider_overall_min = -20.0
    slider_overall_max = 120.0
    slider_step = 0.5

    obs_color_min = Observable(slider_min_val_initial)
    obs_color_max = Observable(slider_max_val_initial)
    color_range_obs = lift((cmin, cmax) -> (min(cmin, cmax - slider_step), max(cmax, cmin + slider_step)), obs_color_min, obs_color_max)

    heatmap!(ax_heatmap, obs_thermal_data, colormap=:hot, colorrange=color_range_obs)
    Colorbar(left_panel[1, 2], colormap=:hot, colorrange=color_range_obs, label="Temperature (°C)", width=30)
    colsize!(left_panel, 2, Auto())


    # Right panel for charts
    right_panel = main_grid[1, 2] = GridLayout()
    ax_temp_chart = Axis(right_panel[1, 1], title="Max Temp Sent to ESP32 (°C)", xlabel="Frame Count / Update", ylabel="Max Temp (°C)")
    ax_duty_chart = Axis(right_panel[2, 1], title="Actuator Duty Cycles (%)", xlabel="Update Count", ylabel="Duty Cycle (0-100%)")

    max_temp_history_x = Observable(Int[])
    max_temp_history_y = Observable(Float32[])
    lines!(ax_temp_chart, max_temp_history_x, max_temp_history_y, color=:blue, label="Max Temp")
    axislegend(ax_temp_chart, position=:lt)


    # Observables for duty cycle charts
    duty_hot_history_x = Observable(Int[])
    duty_hot_history_y = Observable(Float32[])
    duty_cold_history_x = Observable(Int[])
    duty_cold_history_y = Observable(Float32[])

    lines!(ax_duty_chart, duty_hot_history_x, duty_hot_history_y, color=:red, label="Hot Duty (%)")
    lines!(ax_duty_chart, duty_cold_history_x, duty_cold_history_y, color=:cyan, label="Cold Duty (%)")
    axislegend(ax_duty_chart, position=:lt)

    # Sliders
    label_obs_color_min = lift(x -> @sprintf("%.1f°C", x), obs_color_min)
    label_obs_color_max = lift(x -> @sprintf("%.1f°C", x), obs_color_max)

    Label(slider_grid[1, 1], "Min Temp (Heatmap):", halign=:right)
    sl_min = Slider(slider_grid[1, 2], range=slider_overall_min:slider_step:slider_overall_max, startvalue=slider_min_val_initial)
    Label(slider_grid[1, 3], label_obs_color_min, halign=:left, width=100)

    Label(slider_grid[2, 1], "Max Temp (Heatmap):", halign=:right)
    sl_max = Slider(slider_grid[2, 2], range=slider_overall_min:slider_step:slider_overall_max, startvalue=slider_max_val_initial)
    Label(slider_grid[2, 3], label_obs_color_max, halign=:left, width=100)

    on(sl_min.value) do val
        obs_color_min[] = val
    end
    on(sl_max.value) do val
        obs_color_max[] = val
    end

    colsize!(slider_grid, 1, Auto())
    colsize!(slider_grid, 2, Relative(0.5))
    colsize!(slider_grid, 3, Auto())
    rowsize!(fig.layout, 2, Auto()) # Make slider row height auto

    display(fig)
    println("Makie GUI initialized. Waiting for data...")

    # --- Main Communication Logic ---
    try
        LibSerialPort.open(esp32_port_name, esp32_baud_rate) do esp32_serial_port
            println("ESP32 serial port $esp32_port_name opened successfully.")

            println("Attempting to connect to thermal camera data source at $(ServerConfig.SERVER_IP):$(ServerConfig.SERVER_PORT)...")
            thermal_data_socket = nothing
            receive_buffer = UInt8[]

            try
                thermal_data_socket = connect(ServerConfig.SERVER_IP, ServerConfig.SERVER_PORT)
                println("Successfully connected to thermal data source.")

                packet_count = 0      # For max temp x-axis
                duty_update_count = 0 # Separate counter for duty cycle x-axis
                running = true
                max_history_length = 200

                while running && isopen(thermal_data_socket) && events(fig.scene).window_open[]
                    try
                        data_chunk = readavailable(thermal_data_socket)
                        if !isempty(data_chunk)
                            append!(receive_buffer, data_chunk)
                        elseif !isopen(thermal_data_socket)
                            println("Thermal data socket closed unexpectedly.")
                            running = false
                        end

                        while length(receive_buffer) >= ServerConfig.TCP_TOTAL_FRAME_SIZE && running
                            packet_count += 1
                            current_tcp_frame = receive_buffer[1:ServerConfig.TCP_TOTAL_FRAME_SIZE]
                            receive_buffer = receive_buffer[(ServerConfig.TCP_TOTAL_FRAME_SIZE+1):end]
                            parsed_image_data = ImageProcessor.process_tcp_frame(current_tcp_frame, ServerConfig)

                            if parsed_image_data.is_valid && parsed_image_data.thermal_image_matrix !== nothing
                                raw_img_matrix = parsed_image_data.thermal_image_matrix
                                temp_matrix = raw_matrix_to_celsius_matrix(raw_img_matrix)
                                obs_thermal_data[] = temp_matrix

                                if !isempty(temp_matrix)
                                    max_temp = maximum(temp_matrix)

                                    push!(max_temp_history_x[], packet_count)
                                    push!(max_temp_history_y[], max_temp)
                                    if length(max_temp_history_x[]) > max_history_length
                                        deleteat!(max_temp_history_x[], 1)
                                        deleteat!(max_temp_history_y[], 1)
                                    end
                                    notify(max_temp_history_x)
                                    notify(max_temp_history_y)
                                    autolimits!(ax_temp_chart)

                                    esp_message = @sprintf("SET_HOT_TEMP:%.2f\n", max_temp)
                                    try
                                        LibSerialPort.write(esp32_serial_port, esp_message)
                                    catch e_serial_write
                                        println("Error writing to ESP32: $e_serial_write")
                                        running = false
                                    end
                                end
                            end
                        end

                        try # Read from ESP32
                            if LibSerialPort.bytesavailable(esp32_serial_port) > 0
                                esp_response_bytes = LibSerialPort.read(esp32_serial_port)
                                esp_response_str = String(copy(esp_response_bytes))
                                print("ESP32: $(strip(esp_response_str))\n")

                                # Parse duty cycle messages
                                for line in split(strip(esp_response_str), '\n')
                                    if occursin("ESP32_HOT_DUTY:", line)
                                        try
                                            val_str = split(line, ':')[2]
                                            duty_val = parse(Float32, val_str) * 100.0f0 # Convert to percentage
                                            duty_update_count += 1 # Increment for each duty update received
                                            push!(duty_hot_history_x[], duty_update_count)
                                            push!(duty_hot_history_y[], duty_val)
                                            if length(duty_hot_history_x[]) > max_history_length
                                                deleteat!(duty_hot_history_x[], 1)
                                                deleteat!(duty_hot_history_y[], 1)
                                            end
                                            notify(duty_hot_history_x)
                                            notify(duty_hot_history_y)
                                        catch e_parse
                                            println("Error parsing HOT_DUTY value: $e_parse from line: $line")
                                        end
                                    elseif occursin("ESP32_COLD_DUTY:", line)
                                        try
                                            val_str = split(line, ':')[2]
                                            duty_val = parse(Float32, val_str) * 100.0f0 # Convert to percentage
                                            # Using same duty_update_count for simplicity, or use separate if needed
                                            push!(duty_cold_history_x[], duty_update_count)
                                            push!(duty_cold_history_y[], duty_val)
                                            if length(duty_cold_history_x[]) > max_history_length
                                                deleteat!(duty_cold_history_x[], 1)
                                                deleteat!(duty_cold_history_y[], 1)
                                            end
                                            notify(duty_cold_history_x)
                                            notify(duty_cold_history_y)
                                        catch e_parse
                                            println("Error parsing COLD_DUTY value: $e_parse from line: $line")
                                        end
                                    end
                                end
                                autolimits!(ax_duty_chart) # Re-adjust duty chart limits
                            end
                        catch e_serial_read
                            # println("Error reading from ESP32: $e_serial_read")
                        end
                        sleep(0.01)
                    catch e
                        if e isa InterruptException
                            println("\nInterrupted.")
                            running = false
                        elseif e isa EOFError
                            println("Thermal source EOF.")
                            running = false
                        elseif e isa Base.IOError && !isopen(thermal_data_socket)
                            println("Thermal source IOError & closed: $e.")
                            running = false
                        else
                            println("Error in main loop: $e")
                        end
                    end
                end
            catch e_tcp
                if !(e_tcp isa InterruptException)
                    println("TCP connection error: $e_tcp")
                end
            finally
                if thermal_data_socket !== nothing && isopen(thermal_data_socket)
                    close(thermal_data_socket)
                    println("Thermal data source socket closed.")
                end
            end
        end
    catch e_serial_open
        println("Error opening ESP32 serial port $esp32_port_name: $e_serial_open")
    end
    println("GUI and temperature sending script stopped.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_gui_with_makie()
end
