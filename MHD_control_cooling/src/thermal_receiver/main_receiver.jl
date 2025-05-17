# thermal_receiver/main_receiver_stream.jl

using Sockets
using Printf
using GLMakie
using Dates

include("config.jl")
include("image_processor.jl")

# --- Temperature Conversion Functions ---
"""
Converts a single raw UInt16 pixel value to Celsius temperature.
"""
function raw_to_celsius(raw_value::UInt16)::Float32
    return Float32(Float64(raw_value) * 0.0984 - 265.82)
end

"""
Converts a matrix of raw UInt16 pixel values to a matrix of Celsius temperatures (Float32).
"""
function raw_matrix_to_celsius_matrix(raw_matrix::Matrix{UInt16})::Matrix{Float32}
    # Broadcasting the conversion formula.
    # Intermediate calculations use Float64 for precision before casting to Float32.
    return Float32.((Float64.(raw_matrix) .* 0.0984) .- 265.82)
end
# --- End Temperature Conversion Functions ---

function main_stream_with_makie()
    println("Thermal Data Receiver Client (Temperatures with Makie Visualization)")
    println("Attempting to connect to ESP32 server at $(ServerConfig.SERVER_IP):$(ServerConfig.SERVER_PORT)...")

    socket = nothing
    receive_buffer = UInt8[]
    data_file_stream = nothing
    output_filename = ""

    # --- Makie Setup for Temperatures ---
    # Observable now holds Float32 temperature data
    initial_temp_data = Observable(zeros(Float32, ServerConfig.FRAME_HEIGHT, ServerConfig.FRAME_WIDTH))

    # Adjusted slider ranges for Celsius temperatures
    # Raw 2900 -> ~19.5°C; Raw 3500 -> ~78.6°C
    # Raw 2000 -> ~-69°C; Raw 6000 -> ~324°C
    # Let's set a more common initial and overall range for typical scenes.
    slider_min_val_initial = 15.0  # °C
    slider_max_val_initial = 40.0  # °C
    slider_overall_min = -20.0 # °C (Adjust if you expect very low temperatures)
    slider_overall_max = 120.0 # °C (Adjust if you expect very high temperatures)
    slider_step = 0.5          # °C

    obs_color_min = Observable(slider_min_val_initial)
    obs_color_max = Observable(slider_max_val_initial)

    color_range_obs = lift(obs_color_min, obs_color_max) do cmin, cmax
        (min(cmin, cmax - slider_step), max(cmax, cmin + slider_step))
    end

    fig = Figure(size=(ServerConfig.FRAME_WIDTH * 7, ServerConfig.FRAME_HEIGHT * 7 + 150))
    ax = Axis(fig[1, 1], title="Live Thermal Image (°C)", aspect=DataAspect())
    # Heatmap now displays temperature data
    hm = heatmap!(ax, initial_temp_data, colormap=:hot, colorrange=color_range_obs)
    Colorbar(fig[1, 2], hm, label="Temperature (°C)") # Updated label

    slider_layout = GridLayout(fig[2, 1])
    # Labels for sliders now reflect temperature values
    label_obs_color_min = lift(x -> @sprintf("%.1f°C", x), obs_color_min)
    label_obs_color_max = lift(x -> @sprintf("%.1f°C", x), obs_color_max)

    Label(slider_layout[1, 1], "Min Temp:", halign=:right)
    sl_min = Slider(slider_layout[1, 2], range=slider_overall_min:slider_step:slider_overall_max, startvalue=slider_min_val_initial)
    Label(slider_layout[1, 3], label_obs_color_min, halign=:left, width=100)

    Label(slider_layout[2, 1], "Max Temp:", halign=:right)
    sl_max = Slider(slider_layout[2, 2], range=slider_overall_min:slider_step:slider_overall_max, startvalue=slider_max_val_initial)
    Label(slider_layout[2, 3], label_obs_color_max, halign=:left, width=100)

    on(sl_min.value) do val
        obs_color_min[] = val
    end
    on(sl_max.value) do val
        obs_color_max[] = val
    end

    colsize!(slider_layout, 1, Auto())
    colsize!(slider_layout, 2, Relative(0.6))
    colsize!(slider_layout, 3, Auto())

    display(fig)
    println("Makie plot initialized for temperatures. Waiting for data...")
    # --- End Makie Setup ---

    try
        data_dir = "thermal_recordings_celsius" # New directory for temperature data
        if !isdir(data_dir)
            try
                mkpath(data_dir)
                println("Created directory for recordings: $data_dir")
            catch e
                println("Warning: Could not create directory $data_dir: $e. Recordings will be saved in current directory.")
                data_dir = "."
            end
        end

        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        # Indicate in filename that it stores temperatures, and uses Float32
        output_filename = joinpath(data_dir, "thermal_temps_f32_$(timestamp).bin")

        try
            data_file_stream = open(output_filename, "w")
            println("Attempting to record thermal data (Temperatures, Float32) to: $output_filename")

            if !(0 < ServerConfig.FRAME_WIDTH <= typemax(UInt16)) || !(0 < ServerConfig.FRAME_HEIGHT <= typemax(UInt16))
                @error "Frame dimensions (W:$(ServerConfig.FRAME_WIDTH), H:$(ServerConfig.FRAME_HEIGHT)) exceed UInt16 limits for header. Adjust header saving logic or check ServerConfig."
                close(data_file_stream)
                data_file_stream = nothing
            else
                write(data_file_stream, UInt16(ServerConfig.FRAME_WIDTH))
                write(data_file_stream, UInt16(ServerConfig.FRAME_HEIGHT))
                # Optionally: write a version or data type flag here if format might change often
                # e.g., write(data_file_stream, UInt8(1)) # Version 1: Float32 data
                flush(data_file_stream)
                println("Saved data header to $output_filename: Width=$(ServerConfig.FRAME_WIDTH), Height=$(ServerConfig.FRAME_HEIGHT). Data will be Float32 temperatures.")
            end
        catch e
            println("Error opening file $output_filename for recording or writing header: $e")
            if data_file_stream !== nothing && isopen(data_file_stream)
                close(data_file_stream)
            end
            data_file_stream = nothing
        end

        socket = connect(ServerConfig.SERVER_IP, ServerConfig.SERVER_PORT)
        println("Successfully connected to server on port $(ServerConfig.SERVER_PORT).")
        println("Waiting for raw data stream (will be converted to Celsius)...")

        packet_count = 0
        running = true
        while running && isopen(socket) && isopen(fig.scene)
            try
                data_chunk = readavailable(socket)

                if isempty(data_chunk)
                    if !isopen(socket)
                        println("Socket closed while trying to read data.")
                        running = false
                        continue
                    end
                    sleep(0.001)
                    continue
                end

                append!(receive_buffer, data_chunk)

                while length(receive_buffer) >= ServerConfig.TCP_TOTAL_FRAME_SIZE
                    packet_count += 1
                    current_tcp_frame = receive_buffer[1:ServerConfig.TCP_TOTAL_FRAME_SIZE]
                    receive_buffer = receive_buffer[(ServerConfig.TCP_TOTAL_FRAME_SIZE+1):end]

                    # ImageProcessor still returns raw UInt16 data
                    parsed_image_data = ImageProcessor.process_tcp_frame(current_tcp_frame, ServerConfig)

                    if parsed_image_data.is_valid && parsed_image_data.thermal_image_matrix !== nothing
                        raw_img_matrix = parsed_image_data.thermal_image_matrix # This is Matrix{UInt16}

                        # Convert raw data to temperatures (Matrix{Float32})
                        temp_matrix = raw_matrix_to_celsius_matrix(raw_img_matrix)

                        # Update Makie plot with temperature data
                        initial_temp_data[] = temp_matrix

                        # Save the temperature matrix (Float32) to the file
                        if data_file_stream !== nothing && isopen(data_file_stream)
                            try
                                write(data_file_stream, temp_matrix) # Write Float32 matrix
                                flush(data_file_stream)
                            catch e_write
                                println("Error writing temperature frame $packet_count to file: $e_write")
                                close(data_file_stream)
                                data_file_stream = nothing
                            end
                        end

                        if packet_count % 30 == 0 # Print stats less frequently
                            rows, cols = size(temp_matrix)
                            min_val_temp = minimum(temp_matrix)
                            max_val_temp = maximum(temp_matrix)
                            # Stats now refer to temperatures
                            stats_line = @sprintf "Frame #%d: %dx%d, Temp Min: %.2f°C, Temp Max: %.2f°C, Buf: %d, CRange: (%.1f, %.1f)°C" packet_count rows cols min_val_temp max_val_temp length(receive_buffer) obs_color_min[] obs_color_max[]
                            println(stats_line)
                        end
                    else
                        # @warn "Failed to process extracted TCP frame #$packet_count."
                    end
                end
            catch e
                if !(e isa InterruptException)
                    # ... (existing error handling: EOFError, IOError, etc.)
                    if e isa EOFError
                        println("Connection closed by server (EOF). Stopping receiver.")
                    elseif e isa IOError && !isopen(socket)
                        println("IOError and socket is now closed. Likely connection lost: $e. Stopping receiver.")
                    elseif e isa IOError
                        println("IOError while reading from socket: $e. Connection might be unstable. Stopping receiver.")
                    else
                        println("Unhandled error ($typeof(e)) encountered in packet processing loop: $e")
                    end
                end
                running = false
            end
        end
    catch e
        if !(e isa InterruptException)
            println("Failed to connect to server or other critical network error: $e")
        else
            println("Connection or operation interrupted by user (Ctrl-C).")
        end
    finally
        if socket !== nothing && isopen(socket)
            close(socket)
            println("Socket closed.")
        elseif socket !== nothing
            println("Socket was already closed or never opened successfully.")
        end

        if data_file_stream !== nothing && isopen(data_file_stream)
            try
                close(data_file_stream)
                println("Temperature data recording saved to: $output_filename")
            catch e_close
                println("Error closing data file $output_filename: $e_close")
            end
            # ... (existing logic for data_file_stream already closed or not opened) ...
        elseif data_file_stream !== nothing && !isopen(data_file_stream) && output_filename != ""
            println("Temperature data recording to $output_filename was stopped prematurely. Check file for partial data.")
        elseif output_filename != "" && data_file_stream === nothing
            println("No temperature data was recorded, or file $output_filename could not be properly opened/written to.")
        end

        println("Makie window might remain open; close it manually if script terminates.")
    end
    println("Thermal data receiver client (temperatures) stopped.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_stream_with_makie()
end
