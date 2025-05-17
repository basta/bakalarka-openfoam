# thermal_receiver/packet_parser.jl

"""
Module for parsing thermal data packets from the ESP32 server.
"""
module PacketParser

# Assuming ServerConfig is in a file config.jl in the same directory
# and will be included by the main script before this one.
# If running this module standalone, you might need to adjust paths or include directly.
# For the structure main_receiver.jl -> include("config.jl") -> include("packet_parser.jl"),
# ServerConfig will be in the Main module scope.
# A more robust way for larger projects would be proper module management.
# For this script structure, we rely on ServerConfig being available.

export ParsedPacket, parse_thermal_packet

"""
Structure to hold the parsed data from a thermal packet.
"""
struct ParsedPacket
    is_valid::Bool
    header_prefix_ok::Bool
    payload_len_str::String
    payload_len_val::Union{UInt32, Nothing} # Expected: 0x2808 (10248)
    frame_type::String
    metadata::Vector{UInt8} # 160 bytes
    thermal_image_raw::Vector{UInt8} # 10080 bytes
    thermal_image::Union{Matrix{UInt16}, Nothing} # IMAGE_HEIGHT x IMAGE_WIDTH
    crc_string::String
    crc_value::Union{UInt16, Nothing} # Typically 16-bit for 4 hex chars
end

"""
    parse_thermal_packet(packet_bytes::Vector{UInt8}, cfg::Module)::ParsedPacket

Parses a raw byte vector representing one thermal data packet.
`cfg` should be the ServerConfig module.
"""
function parse_thermal_packet(packet_bytes::Vector{UInt8}, cfg)::ParsedPacket
    local is_overall_valid = true
    local is_prefix_ok = false
    local parsed_payload_len_str = ""
    local parsed_payload_len_val::Union{UInt32, Nothing} = nothing
    local parsed_frame_type = ""
    local parsed_metadata = UInt8[]
    local parsed_thermal_raw = UInt8[]
    local parsed_thermal_image::Union{Matrix{UInt16}, Nothing} = nothing
    local parsed_crc_str = ""
    local parsed_crc_val::Union{UInt16, Nothing} = nothing

    if length(packet_bytes) != cfg.EXPECTED_PACKET_SIZE
        @error "Invalid packet size: $(length(packet_bytes)). Expected $(cfg.EXPECTED_PACKET_SIZE)."
        return ParsedPacket(false, false, "", nothing, "", UInt8[], UInt8[], nothing, "", nothing)
    end

    # Bytes are 0-indexed in documentation, Julia is 1-indexed.
    # 1. Bytes 0-3: "   #"
    expected_prefix = UInt8[' ', ' ', ' ', '#']
    actual_prefix = packet_bytes[1:4]
    is_prefix_ok = (actual_prefix == expected_prefix)
    if !is_prefix_ok
        @warn "Packet prefix mismatch. Expected $expected_prefix, got $actual_prefix"
        is_overall_valid = false
    end

    # 2. Bytes 4-7: "2808" (ASCII hex string for payload length 0x2808 = 10248)
    parsed_payload_len_str = String(packet_bytes[5:8])
    try
        parsed_payload_len_val = parse(UInt32, parsed_payload_len_str, base=16)
        if parsed_payload_len_val != 0x2808
            @warn "Payload length field value mismatch. Expected 0x2808, got 0x$(string(parsed_payload_len_val, base=16, pad=4)) (from \"$parsed_payload_len_str\")"
            is_overall_valid = false
        end
    catch e
        @warn "Could not parse payload length string: \"$parsed_payload_len_str\". Error: $e"
        is_overall_valid = false
        parsed_payload_len_val = nothing
    end

    # 3. Bytes 8-11: "GFRA" (Frame type)
    parsed_frame_type = String(packet_bytes[9:12])
    if parsed_frame_type != "GFRA"
        @warn "Unexpected frame type. Expected \"GFRA\", got \"$parsed_frame_type\""
        # Not necessarily invalidating the whole packet structure, but it's a warning.
    end

    # 4. Bytes 12-171 (0-indexed): 160 bytes of metadata/padding -> Julia indices 13:172
    parsed_metadata = packet_bytes[13:172]

    # 5. Bytes 172-10251 (0-indexed): 10080 bytes of thermal image data -> Julia indices 173:10252
    parsed_thermal_raw = packet_bytes[173:10252]

    expected_raw_image_size = cfg.IMAGE_WIDTH * cfg.IMAGE_HEIGHT * sizeof(UInt16)
    if length(parsed_thermal_raw) == expected_raw_image_size
        try
            # ESP32 is little-endian.
            # Reinterpret the UInt8 raw vector as a vector of UInt16 (creates a view).
            temp_u16_vector_view = reinterpret(UInt16, parsed_thermal_raw)
            # Convert to host endianness (ltoh: little-endian to host) and make a concrete copy.
            thermal_data_uint16_flat = ltoh.(temp_u16_vector_view)

            # Reshape the flat vector into a matrix.
            # Data is assumed to be in row-major order in the flat vector.
            # reshape(vector, W, H) creates a matrix with W rows, H columns (column-major fill).
            # permutedims(..., (2,1)) transposes it to H rows, W columns.
            parsed_thermal_image = permutedims(reshape(thermal_data_uint16_flat, cfg.IMAGE_WIDTH, cfg.IMAGE_HEIGHT), (2,1))
        catch e
            @error "Could not convert or reshape thermal image data. Error: $e"
            parsed_thermal_image = nothing
            is_overall_valid = false
        end
    else
        @warn "Thermal raw data size mismatch. Expected $expected_raw_image_size, got $(length(parsed_thermal_raw))"
        is_overall_valid = false
    end

    # 6. Bytes 10252-10255 (0-indexed): 4-char ASCII CRC -> Julia indices 10253:10256
    parsed_crc_str = String(packet_bytes[10253:10256])
    try
        parsed_crc_val = parse(UInt16, parsed_crc_str, base=16)
    catch e
        @warn "Could not parse CRC string: \"$parsed_crc_str\". Error: $e"
        # Not invalidating packet based on CRC parse failure alone for now,
        # as CRC validation is not yet implemented.
        parsed_crc_val = nothing
    end

    # Final validity check (can be more stringent if needed)
    is_overall_valid = is_overall_valid && (parsed_thermal_image !== nothing)

    return ParsedPacket(
        is_overall_valid,
        is_prefix_ok,
        parsed_payload_len_str,
        parsed_payload_len_val,
        parsed_frame_type,
        parsed_metadata,
        parsed_thermal_raw,
        parsed_thermal_image,
        parsed_crc_str,
        parsed_crc_val
    )
end

end # module PacketParser
