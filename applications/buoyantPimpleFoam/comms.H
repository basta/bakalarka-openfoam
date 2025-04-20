// #############################################################################
// # File: comms.H - Modifications for sending temperature data              #
// #############################################################################
#ifndef COMMS_H // Include guard
#define COMMS_H

#include "fvCFD.H" // Include OpenFOAM headers for volScalarField etc.
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <memory> // For std::unique_ptr
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <errno.h> // For errno
#include <cstring> // For strerror

class SocketClient {
private:
    int sock;
    struct sockaddr_in server;
    bool connected;
    std::string serverIP;
    int port;

public:
    // Constructor: Stores connection details but doesn't connect yet
    SocketClient(const std::string& ip, int p) :
        sock(-1), connected(false), serverIP(ip), port(p)
    {
        // Prepare the sockaddr_in structure
        server.sin_addr.s_addr = inet_addr(serverIP.c_str());
        server.sin_family = AF_INET;
        server.sin_port = htons(port);
    }

    // Destructor: Ensures disconnection and socket closure
    ~SocketClient() {
        disconnect();
    }

    // Attempts to establish a connection
    bool connect() {
        if (connected) {
            // Optional: Add a check to see if the connection is still alive
            // For now, assume it is if connected is true.
            return true;
        }

        // Create socket
        sock = socket(AF_INET, SOCK_STREAM, 0);
        if (sock == -1) {
            std::cerr << "Could not create socket: " << strerror(errno) << std::endl;
            return false;
        }

        // Connect to remote server
        if (::connect(sock, (struct sockaddr *)&server, sizeof(server)) < 0) {
            Foam::Warning << "Connect failed to " << serverIP << ":" << port
                          << ". Error: " << strerror(errno) << " (errno: " << errno << ")" << Foam::endl;
            close(sock);
            sock = -1;
            connected = false;
            return false;
        }

        Foam::Info << "Connected to server " << serverIP << ":" << port << Foam::endl;
        connected = true;
        return true;
    }

    // Disconnects and closes the socket
    void disconnect() {
        if (sock != -1) {
            Foam::Info << "Closing socket" << Foam::endl;
            close(sock);
            sock = -1;
        }
        connected = false;
    }

    // Sends a generic request string
    bool sendRequest(const std::string& request) {
        if (!connected && !connect()) { // Try to reconnect if not connected
             Foam::Warning << "Send failed: Not connected and cannot reconnect." << Foam::endl;
             return false;
        }

        ssize_t sentBytes = send(sock, request.c_str(), request.length(), 0);
        if (sentBytes < 0) {
            Foam::Warning << "Send failed with error: " << strerror(errno) << " (errno: " << errno << ")" << Foam::endl;
            disconnect(); // Assume connection is broken
            return false;
        }
         if (static_cast<size_t>(sentBytes) != request.length()) {
             Foam::Warning << "Send failed: Incomplete send." << Foam::endl;
             disconnect(); // Assume connection is broken
             return false;
         }
        return true;
    }

    // Sends a header string, followed by size, then float vector data
    bool sendData(const char* header, const std::vector<float>& data) {
        if (!connected && !connect()) { // Try to reconnect
            Foam::Warning << "Send data failed: Not connected and cannot reconnect." << Foam::endl;
            return false;
        }

        // 1. Send Header
        std::string fullHeader = header; // Assumes header includes newline if needed
        if (!sendRequest(fullHeader)) {
             Foam::Warning << "Failed to send data header: " << header << Foam::endl;
             // sendRequest handles disconnect on error
             return false;
        }

        // 2. Send Size (as int)
        int dataSize = data.size();
        ssize_t sentBytes = send(sock, &dataSize, sizeof(int), 0);
        if (sentBytes != sizeof(int)) {
             Foam::Warning << "Failed to send data size. Sent " << sentBytes << " bytes, expected " << sizeof(int) << ". Error: " << strerror(errno) << Foam::endl;
             disconnect();
             return false;
        }

        // 3. Send Data (vector of floats)
        size_t totalBytesToSend = data.size() * sizeof(float);
        size_t totalBytesSent = 0;
        const char* dataPtr = reinterpret_cast<const char*>(data.data());

        while (totalBytesSent < totalBytesToSend) {
            sentBytes = send(sock, dataPtr + totalBytesSent, totalBytesToSend - totalBytesSent, 0);
            if (sentBytes < 0) {
                Foam::Warning << "Send data failed during transmission. Error: " << strerror(errno) << Foam::endl;
                disconnect();
                return false;
            }
             if (sentBytes == 0) {
                 Foam::Warning << "Send data failed: Socket closed during transmission." << Foam::endl;
                 disconnect();
                 return false;
            }
            totalBytesSent += sentBytes;
        }

        Foam::Info << "Sent " << dataSize << " floats (" << totalBytesSent << " bytes) with header: " << header << Foam::endl;
        return true;
    }


    // Receives a dynamic array of floats (size first, then data)
    std::vector<float> receiveArray() {
        std::vector<float> result;
         if (!connected && !connect()) { // Try to reconnect
             Foam::Warning << "Receive array failed: Not connected and cannot reconnect." << Foam::endl;
             return result; // Return empty vector
         }

        // 1. Receive the array size (number of floats)
        int arraySize = 0;
        ssize_t bytesReceived = recv(sock, &arraySize, sizeof(int), MSG_WAITALL); // Wait for all size bytes
        if (bytesReceived <= 0) {
            Foam::Warning << "Failed to receive array size. Return: " << bytesReceived << ", Error: " << strerror(errno) << Foam::endl;
            disconnect();
            return result; // Return empty vector
        }
         if (bytesReceived != sizeof(int)) {
             Foam::Warning << "Failed to receive complete array size. Received " << bytesReceived << " bytes." << Foam::endl;
             disconnect();
             return result; // Return empty vector
         }

        if (arraySize <= 0) {
             Foam::Warning << "Received invalid array size: " << arraySize << Foam::endl;
             // Might still be connected, but data is bad. Decide whether to disconnect.
             // For now, just return empty.
             return result;
        }


        // 2. Receive the array data
        size_t totalBytesToReceive = static_cast<size_t>(arraySize) * sizeof(float);
        result.resize(arraySize); // Pre-allocate vector
        size_t totalBytesReceived = 0;
        char* bufferPtr = reinterpret_cast<char*>(result.data());

        while (totalBytesReceived < totalBytesToReceive) {
            bytesReceived = recv(sock, bufferPtr + totalBytesReceived, totalBytesToReceive - totalBytesReceived, 0); // Don't necessarily need MSG_WAITALL here, loop handles partial reads

            if (bytesReceived < 0) {
                Foam::Warning << "Error receiving array data: " << strerror(errno) << Foam::endl;
                disconnect();
                result.clear(); // Clear partial data
                return result;
            }
            if (bytesReceived == 0) {
                 Foam::Warning << "Error receiving array data: Socket closed prematurely." << Foam::endl;
                 disconnect();
                 result.clear(); // Clear partial data
                 return result;
            }
            totalBytesReceived += bytesReceived;
        }

         Foam::Info << "Received " << arraySize << " floats (" << totalBytesReceived << " bytes)." << Foam::endl;
        return result;
    }

    // Check if currently connected
    bool isConnected() const {
        return connected;
    }
};

// --- Global Socket Client Pointer ---
// Using unique_ptr for automatic memory management
std::unique_ptr<SocketClient> sock_client;

// --- Initialization Function ---
// Creates and connects the socket client. Call this once in main().
bool init_socket(const std::string& ip = "127.0.0.1", int port = 8080) {
    if (!sock_client) { // Only initialize if not already done
        try {
            sock_client = std::make_unique<SocketClient>(ip, port);
            if (!sock_client->connect()) {
                 Foam::Warning << "Initial socket connection failed." << Foam::endl;
                 sock_client.reset(); // Release the pointer if connection failed
                 return false;
            }
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Exception during socket initialization: " << e.what() << Foam::endl;
            sock_client.reset();
            return false;
        }
    }
    // If already initialized, try connecting if needed
    return sock_client->connect();
}

// --- Function to Request Force Field ---
// Reuses the global sock_client connection
std::vector<std::array<float, 3>> request_field(float t) {
    std::vector<std::array<float, 3>> result;
    if (!sock_client || !sock_client->isConnected()) {
        Foam::Warning << "Cannot request field: Socket not initialized or connected." << Foam::endl;
        // Attempt to re-initialize/reconnect
        if (!init_socket()) {
             return result; // Return empty if reconnect fails
        }
    }

    // Send request
    std::string req_string = "REQUEST_FORCES_" + std::to_string(t) + "\n"; // Changed header for clarity
    if (!sock_client->sendRequest(req_string)) {
        Foam::Warning << "Failed to send force field request." << Foam::endl;
        // sendRequest handles disconnect on error
        return result; // Return empty
    }

    // Receive data
    std::vector<float> dataArray = sock_client->receiveArray();
    if (dataArray.empty() || dataArray.size() % 3 != 0) {
         Foam::Warning << "Received invalid or empty force field data." << Foam::endl;
         // receiveArray handles disconnect on error
         return result; // Return empty
    }

    // Convert flat array to vector of 3D vectors
    result.reserve(dataArray.size() / 3);
    for (size_t i = 0; i < dataArray.size(); i += 3) {
        result.push_back({dataArray[i], dataArray[i + 1], dataArray[i + 2]});
    }
    Foam::Info << "Successfully received and processed force field for t=" << t << Foam::endl;
    return result;
}


// --- Function to Send Temperature Field ---
// Extracts data from volScalarField and sends using the global sock_client
bool send_temperature(const Foam::volScalarField& T, float currentTime) {
     if (!sock_client || !sock_client->isConnected()) {
        Foam::Warning << "Cannot send temperature: Socket not initialized or connected." << Foam::endl;
         // Attempt to re-initialize/reconnect
         if (!init_socket()) {
             return false; // Return false if reconnect fails
         }
    }

    // 1. Extract Temperature Data
    // Access the internal field data which is a Foam::Field<scalar>
    const Foam::Field<Foam::scalar>& internalT = T.internalField();

    // Convert Foam::scalar (usually double) to float for sending
    // Note: This involves potential precision loss if scalar is double.
    std::vector<float> tempData(internalT.size());
    for (Foam::label i = 0; i < internalT.size(); ++i) {
        tempData[i] = static_cast<float>(internalT[i]);
    }

    // 2. Send Data
    std::string header = "SEND_TEMP_" + std::to_string(currentTime) + "\n";
    if (!sock_client->sendData(header.c_str(), tempData)) {
        Foam::Warning << "Failed to send temperature field." << Foam::endl;
        // sendData handles disconnect on error
        return false;
    }

    Foam::Info << "Successfully sent temperature field (" << tempData.size() << " cells) for t=" << currentTime << Foam::endl;
    return true;
}

// --- Function to Send Wall Heat Flux Data ---
// Sends pre-calculated heat flux values using the global sock_client
bool send_wall_heat_flux(const std::vector<float>& flux_data, float currentTime) {
     if (!sock_client || !sock_client->isConnected()) {
        Foam::Warning << "Cannot send wall heat flux: Socket not initialized or connected." << Foam::endl;
         // Attempt to re-initialize/reconnect
         if (!init_socket()) {
             return false; // Return false if reconnect fails
         }
    }

    // Send Data
    std::string header = "SEND_WALLHEATFLUX_" + std::to_string(currentTime) + "\n";
    if (!sock_client->sendData(header.c_str(), flux_data)) {
        Foam::Warning << "Failed to send wall heat flux data." << Foam::endl;
        // sendData handles disconnect on error
        return false;
    }

    // Foam::Info << "Successfully sent wall heat flux data (" << flux_data.size() << " values) for t=" << currentTime << Foam::endl; // Redundant with sendData info
    return true;
}


// --- Example usage function ---
int example() {
    if (!init_socket()) { // Initialize connection
        return 1;
    }

    // Example: Request forces at t=0.1
    auto forces = request_field(0.1f);
    if (!forces.empty()) {
         Foam::Info << "Example received " << forces.size() << " force vectors." << Foam::endl;
    }

    // Example: Send a dummy temperature field
    // In real code, you'd pass the actual thermo.T() field
    // Create a dummy field for demonstration
    // Note: This requires a mesh object, which isn't available here.
    // This part is just illustrative.
    /*
    Foam::fvMesh& mesh = ...; // Need a mesh object
    Foam::volScalarField dummyT
    (
        Foam::IOobject
        (
            "dummyT",
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE
        ),
        mesh,
        Foam::dimensionedScalar("dummyT", Foam::dimTemperature, 0.0)
    );
    // Populate dummyT if needed...
    send_temperature(dummyT, 0.1f);
    */

    return 0;
}

#endif // COMMS_H
