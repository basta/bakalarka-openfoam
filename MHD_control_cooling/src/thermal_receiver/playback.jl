# playback_thermal_data.jl

using Printf # For formatted printing

"""
    read_thermal_frames(filename::String)

Reads thermal image frames from a binary file created by the thermal_receiver script.
The file format is expected to be:
1. Frame Width (UInt16)
2. Frame Height (UInt16)
3. Sequentially stored frames, where each frame consists of (width * height) UInt16 pixel values.
"""
function read_thermal_frames(filename::String)
    if !isfile(filename)
        println("Error: File not found: $filename")
        return
    end

    println("Attempting to read thermal data from: $filename")
    file_stream = nothing
    frames_read = 0
    frame_width::UInt16 = 0
    frame_height::UInt16 = 0

    try
        file_stream = open(filename, "r") # Open in binary read mode

        # Read header: Frame Width and Frame Height
        try
            frame_width = read(file_stream, UInt16)
            frame_height = read(file_stream, UInt16)
        catch e
            if e isa EOFError
                println("Error: Could not read header. File is too short or corrupted.")
                return
            else
                rethrow()
            end
        end

        println("Successfully read header:")
        println("  Frame Width: $frame_width pixels")
        println("  Frame Height: $frame_height pixels")

        if frame_width == 0 || frame_height == 0
            println("Error: Invalid frame dimensions in header (width or height is zero). Cannot process frames.")
            return
        end

        # Calculate the number of pixels per frame (convert dimensions to Int for calculations)
        num_pixels_per_frame = Int(frame_width) * Int(frame_height)
        println("Each frame should contain $num_pixels_per_frame pixels (UInt16 values).")

        # Calculate expected bytes per frame for more robust EOF checking if needed
        bytes_per_frame = num_pixels_per_frame * sizeof(UInt16)

        println("\nReading frames...")

        while !eof(file_stream)
            # Allocate a flat buffer for one frame's data
            frame_data_flat = Vector{UInt16}(undef, num_pixels_per_frame)

            try
                # Read the raw pixel data for one frame into the buffer
                # read! will fill the array. It throws EOFError if the stream ends prematurely.
                read!(file_stream, frame_data_flat)
            catch e
                if e isa EOFError
                    # This occurs if the file ends before a full frame can be read.
                    # Check how many bytes were actually available if needed for diagnostics.
                    # bytes_available_at_eof = bytesavailable(file_stream) # This would be 0 if EOFError was thrown by read!
                    # For simplicity, we just note that the last frame might be partial.
                    println("\nEOF reached while attempting to read frame $(frames_read + 1).")
                    println("The file might be truncated, or this is the end of the stream after the last complete frame.")
                    break # Exit the loop, as a full frame could not be read
                else
                    rethrow() # Propagate other types of errors
                end
            end

            frames_read += 1

            # Reshape the flat 1D data into a 2D matrix (height x width)
            # Julia's arrays are column-major, and `reshape` fills them column by column by default.
            # This matches the order in which `write(io, matrix)` serializes matrix data.
            thermal_image_matrix = reshape(frame_data_flat, (Int(frame_height), Int(frame_width)))

            # Process the frame (example: print statistics)
            min_val = minimum(thermal_image_matrix)
            max_val = maximum(thermal_image_matrix)
            # sum promotes to a larger integer type (e.g., Int64) before division
            mean_val = Float64(sum(thermal_image_matrix)) / num_pixels_per_frame

            @printf "Frame #%04d: Min=%5d, Max=%5d, Mean=%.2f\n" frames_read min_val max_val mean_val

            # --- Optional: Add visualization here if desired ---
            # Example: Display the first few frames using GLMakie (uncomment to use)
            # if frames_read <= 3
            #   try
            #     # This would require `using GLMakie` at the top of the file
            #     # Ensure GLMakie is in your Project.toml or environment
            #     makie_fig = Figure()
            #     makie_ax = Axis(makie_fig[1,1], aspect=DataAspect(), title="Frame $frames_read: $frame_height x $frame_width")
            #     heatmap!(makie_ax, thermal_image_matrix, colormap=:hot) # Use desired colormap
            #     display(makie_fig)
            #     println("Displaying frame $frames_read. Close window or press Enter in console to continue...")
            #     # readline() # Uncomment if you want to pause after each displayed frame
            #   catch vis_error
            #     println("Could not display frame $frames_read with GLMakie: $vis_error. (GLMakie might not be loaded or working).")
            #   end
            # end
            # --- End Optional Visualization ---
        end

        println("\nFinished reading.")
        if frames_read > 0
            println("Total complete frames processed: $frames_read")
        else
            println("No complete frames were found or processed after the header.")
        end

    catch e
        # Catch errors not handled by more specific blocks (e.g., permission issues)
        println("An unexpected error occurred: $e")
        # For debugging, you might want to see the stack trace:
        # Base.showerror(stderr, e, catch_backtrace())
    finally
        if file_stream !== nothing && isopen(file_stream)
            close(file_stream)
            println("File stream closed.")
        end
    end
end

# This part allows running the script from the command line with a filename argument
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) == 1
        data_filename = ARGS[1]
        read_thermal_frames(data_filename)
    else
        println("Usage: julia playback_thermal_data.jl <path_to_thermal_data.bin>")

        # Helper: List .bin files in a common recording directory
        default_recording_dir = "thermal_recordings"
        println("\nLooking for recordings in './$default_recording_dir/'...")
        if isdir(default_recording_dir)
            bin_files = filter(f -> endswith(f, ".bin"), readdir(default_recording_dir))
            if !isempty(bin_files)
                println("Available recordings found:")
                for (i, f_name) in enumerate(bin_files)
                    println("  $i. $f_name")
                end
                println("\nExample: julia playback_thermal_data.jl $default_recording_dir/$(bin_files[1])")

            else
                println("  (No .bin files found in '$default_recording_dir')")
            end
        else
            println("  (Directory '$default_recording_dir' not found)")
        end
    end
end
