# thermal_receiver/config.jl

"""
Configuration for the ESP32 Thermal Server Client (Raw Stream Mode).
Based on analysis of client.js.
"""
module ServerConfig

# !!! IMPORTANT: Verify ESP32_HOST is correct for your setup if not 192.168.4.1 !!!
const SERVER_IP = "192.168.4.1" # Default from client.js
const SERVER_PORT = 3333        # From client.js

const FRAME_WIDTH = 80
const FRAME_HEIGHT = 62
const RAW_IMAGE_SIZE_BYTES = FRAME_WIDTH * FRAME_HEIGHT * 2 # 9920 bytes (80 * 62 * 2)

const STRIP_HEAD_BYTES = 160
# STRIP_TAIL_BYTES is implied by TCP_TOTAL_FRAME_SIZE - STRIP_HEAD_BYTES - RAW_IMAGE_SIZE_BYTES
const STRIP_TAIL_BYTES = 160
const TCP_TOTAL_FRAME_SIZE = RAW_IMAGE_SIZE_BYTES + STRIP_HEAD_BYTES + STRIP_TAIL_BYTES # 10240 bytes

end
