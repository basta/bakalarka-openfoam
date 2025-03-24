#include <iostream>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>

class SocketClient {
private:
    int sock;
    struct sockaddr_in server;

public:
    SocketClient(const char* serverIP, int port) {
        // Create socket
        sock = socket(AF_INET, SOCK_STREAM, 0);
        if (sock == -1) {
            throw std::runtime_error("Could not create socket");
        }
        // Prepare the sockaddr_in structure
        server.sin_addr.s_addr = inet_addr(serverIP);
        server.sin_family = AF_INET;
        server.sin_port = htons(port);
    }

    ~SocketClient() {
        std::cout << "Closing socket" << std::endl;
        close(sock);
    }

    bool connect() {
        // Connect to remote server
        if (::connect(sock, (struct sockaddr *)&server, sizeof(server)) < 0) {
            std::cerr << "Connect failed. Error: " << strerror(errno) << " (errno: " << errno << ")" << std::endl;
            return false;
        }
    
        std::cout << "Connected to server" << std::endl;
        return true;
    }
    bool sendRequest(const std::string& request) {
        // Send the request message
        auto err = send(sock, request.c_str(), request.length(), 0);
        if (err < 0) {
            std::cerr << "Send failed with error: " << strerror(errno) << " (errno: " << errno << ")" << std::endl;
            return false;
        }

        return true;
    } 

    std::vector<float> receiveArray() {
        std::vector<float> result;

        // First receive the array size
        int arraySize = 0;
        if (recv(sock, &arraySize, sizeof(int), 0) <= 0) {
            std::cerr << "Failed to receive array size" << std::endl;
            return result;
        }

        auto floatBuffer = std::make_unique<float[]>(512);

        int ctr = 0;
        while (true) {
            // int bytesReceived = recv(sock, buffer + totalBytesReceived, bytesRemaining, 0);
            int bytesReceived = recv(sock, floatBuffer.get(), 512*4, 0);
            for (int i = 0; i < bytesReceived / sizeof(float); i++) {
                result.push_back(floatBuffer.get()[i]);
            }


            if (bytesReceived < 0) {
                std::cerr << "Error receiving array data" << std::endl;
                result.clear();
                return result;
            }
            if (bytesReceived == 0) {
                break;
            }

            // totalBytesReceived += bytesReceived;
            // bytesRemaining -= bytesReceived;
        }

        return result;
    }
};

std::unique_ptr<SocketClient> sock_client;

int init() {

}

std::vector<std::array<float, 3>> request_field(float t) {
    try {
        // Create socket client connecting to localhost on port 8080
        sock_client = std::make_unique<SocketClient>("127.0.0.1", 8080);

        // Connect to server
        if (!sock_client->connect()) {
            return {};
        }

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return {};
    }

    // TODO
    // if (!sock_client->sendRequest("REQUEST_ARRAY " + std::to_string(t))) {
    std::string req_string = "REQUEST_ARRAY_" + std::to_string(t) + "\n";
    if (!sock_client->sendRequest(req_string)) {
        Warning << "Failed to send request" << endl;
        return {};
    }
    std::vector<float> dataArray = sock_client->receiveArray();
    std::vector<std::array<float, 3>> result;
    for (size_t i = 0; i < dataArray.size(); i += 3) {
        result.push_back({dataArray[i], dataArray[i + 1], dataArray[i + 2]});
    }
    return result;
}

int example() {
    try {
        // Create socket client connecting to localhost on port 8080
        SocketClient client("127.0.0.1", 8080);

        // Connect to server
        if (!client.connect()) {
            return 1;
        }

        // Send request for data
        if (!client.sendRequest("REQUEST_ARRAY\n")) {
            return 1;
        }

        // Receive array of integers
        std::vector<float> dataArray = client.receiveArray();

        // Print received data
        std::cout << "Received array data:" << std::endl;
        for (size_t i = 0; i < dataArray.size(); i++) {
            std::cout << dataArray[i] << " ";
            if ((i + 1) % 10 == 0) std::cout << std::endl;
        }
        std::cout << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
