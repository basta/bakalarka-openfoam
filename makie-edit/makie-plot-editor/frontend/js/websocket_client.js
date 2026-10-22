// frontend/js/websocket_client.js

let socket = null;

function connect(url) {
    if (socket && (socket.readyState === WebSocket.OPEN || socket.readyState === WebSocket.CONNECTING)) {
        console.log("WebSocket is already connected or connecting.");
        return;
    }

    socket = new WebSocket(url);

    socket.onopen = function(event) {
        console.log("WebSocket connection established.");
        // You could send an initial message if needed
        // sendMessage({ type: "CLIENT_HELLO", payload: { message: "Frontend connected!" } });
    };

    socket.onmessage = function(event) {
        console.log("Message from server: ", event.data);
        try {
            const message = JSON.parse(event.data);
            handleServerMessage(message);
        } catch (e) {
            console.error("Failed to parse server message or handle it:", e);
            // Handle non-JSON messages or parsing errors if necessary
        }
    };

    socket.onerror = function(error) {
        console.error("WebSocket Error: ", error);
    };

    socket.onclose = function(event) {
        console.log("WebSocket connection closed: ", event.code, event.reason);
        socket = null;
        // Optionally, try to reconnect
        // setTimeout(() => connect(url), 5000);
    };
}

function sendMessage(messageObject) {
    if (socket && socket.readyState === WebSocket.OPEN) {
        try {
            const messageString = JSON.stringify(messageObject);
            socket.send(messageString);
            console.log("Message sent to server: ", messageObject);
        } catch (e) {
            console.error("Failed to stringify message or send:", e);
        }
    } else {
        console.error("WebSocket is not connected. Cannot send message.");
    }
}

function handleServerMessage(message) {
    // Dispatch messages based on their type
    // This is where you'll integrate with ui.js or plot_display.js
    console.log("Handling server message:", message);
    switch (message.type) {
        case "PLOT_UPDATE":
            if (window.PlotDisplayModule && window.PlotDisplayModule.updatePlot) {
                window.PlotDisplayModule.updatePlot(message.payload);
            }
            break;
        case "INITIAL_PLOT_STATE":
            // Could update UI controls and plot display
            if (window.PlotDisplayModule && window.PlotDisplayModule.updatePlot) {
                window.PlotDisplayModule.updatePlot(message.payload.plotUpdate);
            }
            if (window.UiModule && window.UiModule.updateInspectorControls) {
                // You'd need a function in ui.js to take attributes and set control values
                // window.UiModule.updateInspectorControls(message.payload.attributes);
            }
            console.log("Received initial plot state.");
            break;
        case "STATUS_MESSAGE":
            // Display status to user (e.g., in the status bar)
            console.log(`Status: [${message.payload.level}] ${message.payload.message}`);
            const statusBar = document.querySelector('.status-bar'); // Assuming you have one
            if (statusBar) {
                statusBar.textContent = `Status: ${message.payload.message}`;
            }
            break;
        case "ERROR_MESSAGE":
            console.error(`Server Error: ${message.payload.message}`, message.payload.details || '');
            // Display error to user
            alert(`Error from server: ${message.payload.message}`);
            break;
        // Add more cases as per your protocol
        default:
            console.warn("Received unknown message type from server:", message.type);
    }
}

function disconnect() {
    if (socket) {
        socket.close();
    }
}

// Expose functions
window.WebSocketClient = {
    connect,
    sendMessage,
    disconnect
};
