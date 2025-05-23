
"""
Module for processing raw thermal image data frames based on client.js logic.
"""
module ImageProcessor

# This module will rely on ServerConfig being included and available
# in the main script's global scope.

export ThermalImage, process_tcp_frame

"""
Structure to hold the processed raw thermal image.
"""
struct ThermalImage
    is_valid::Bool
    thermal_image_matrix::Union{Matrix{UInt16}, Nothing} # FRAME_HEIGHT x FRAME_WIDTH
    # Optionally, keep a view of the raw bytes for debugging
    # raw_bytes_view::Union{SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int64}}, true}, Nothing}
end

"""
    process_tcp_frame(frame_bytes::Vector{UInt8}, cfg_raw::Module)::ThermalImage

Extracts the raw thermal image data from a full TCP frame.
The input `frame_bytes` is expected to be exactly `cfg_raw.TCP_TOTAL_FRAME_SIZE`.
"""
function process_tcp_frame(frame_bytes::Vector{UInt8}, cfg_raw)::ThermalImage
    # This function assumes frame_bytes is exactly one complete TCP_TOTAL_FRAME_SIZE
    # The calling function (main_receiver) will handle buffering and providing such frames.

    # Calculate start and end indices for the raw image data within the total frame
    # Julia is 1-indexed.
    image_data_start_idx = cfg_raw.STRIP_HEAD_BYTES + 1
    image_data_end_idx = cfg_raw.STRIP_HEAD_BYTES + cfg_raw.RAW_IMAGE_SIZE_BYTES

    if image_data_end_idx > length(frame_bytes)
        @error "process_tcp_frame: Not enough bytes in provided frame_bytes to extract image. Expected at least $image_data_end_idx, got $(length(frame_bytes))."
        return ThermalImage(false, nothing)
    end

    # Extract the raw image data part
    raw_image_data_bytes_segment = view(frame_bytes, image_data_start_idx:image_data_end_idx)

    if length(raw_image_data_bytes_segment) != cfg_raw.RAW_IMAGE_SIZE_BYTES
        @warn "process_tcp_frame: Extracted raw image data segment size mismatch. Expected $(cfg_raw.RAW_IMAGE_SIZE_BYTES), got $(length(raw_image_data_bytes_segment))"
        return ThermalImage(false, nothing)
    end

    local thermal_matrix::Union{Matrix{UInt16}, Nothing} = nothing
    try
        # ESP32 is typically little-endian.
        # Reinterpret the UInt8 raw vector view as a vector of UInt16.
        # Need to ensure raw_image_data_bytes_segment is a concrete array for reinterpret if it's a view from a larger buffer
        # that might change. However, here it's a view of `frame_bytes` which is a full, distinct frame.
        temp_u16_vector_view = reinterpret(UInt16, raw_image_data_bytes_segment)

        # Convert to host endianness (ltoh: little-endian to host) and make a concrete copy.
        thermal_data_uint16_flat = ltoh.(temp_u16_vector_view)

        # Reshape the flat vector into a matrix.
        # Data is assumed to be in row-major order in the flat vector from the sensor.
        # Julia's reshape fills column-major, so we reshape to (width, height)
        # then permute dimensions to get (height, width).
        thermal_matrix = permutedims(reshape(thermal_data_uint16_flat, cfg_raw.FRAME_WIDTH, cfg_raw.FRAME_HEIGHT), (2,1))
    catch e
        @error "process_tcp_frame: Could not convert or reshape thermal image data. Error: $e"
        # Consider logging stacktrace(catch_backtrace()) for more detail
        return ThermalImage(false, nothing)
    end

    return ThermalImage(true, thermal_matrix)
end

end # module ImageProcessor
