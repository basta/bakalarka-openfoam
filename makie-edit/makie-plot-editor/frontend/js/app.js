// frontend/js/app.js

document.addEventListener('DOMContentLoaded', () => {
    console.log("DOM fully loaded and parsed");

    // Build the inspector tree view
    if (window.UiModule && window.UiModule.buildTree) {
        window.UiModule.buildTree(window.UiModule.TREE_DATA, 'inspectorTreeView');
    } else {
        console.error("UiModule or UiModule.buildTree is not defined.");
    }

    // Setup sidebar toggle functionality
    if (window.UiModule && window.UiModule.setupSidebarToggle) {
        window.UiModule.setupSidebarToggle();
    } else {
        console.error("UiModule or UiModule.setupSidebarToggle is not defined.");
    }

    // Setup menu dropdown functionality
    if (window.UiModule && window.UiModule.setupMenuDropdowns) {
        window.UiModule.setupMenuDropdowns();
    } else {
        console.error("UiModule or UiModule.setupMenuDropdowns is not defined.");
    }

    // Initialize WebSocket connection (example)
    if (window.WebSocketClient && window.WebSocketClient.connect) {
        window.WebSocketClient.connect('ws://localhost:8080/ws'); // Replace with your WebSocket server URL
        console.log("WebSocketClient connect function is available. Call it when backend is ready.");
    } else {
        console.error("WebSocketClient or WebSocketClient.connect is not defined.");
    }

    // Example: Add event listeners to dynamically created form elements within the tree
    // This needs to be more robust, possibly using event delegation on the tree container
    // or re-attaching listeners after tree rebuild/content update.
    const inspectorPanel = document.getElementById('controlPanel');
    if (inspectorPanel) {
        inspectorPanel.addEventListener('change', function(event) {
            const target = event.target;
            if (target.classList.contains('form-input') || target.classList.contains('form-select') || target.classList.contains('form-color-picker') || target.classList.contains('form-checkbox')) {
                console.log(`Inspector control changed: ID=${target.id}, Value=${target.type === 'checkbox' ? target.checked : target.value}`);
                // Here you would typically get the associated tree node ID/attribute
                // and send a SET_PLOT_ATTRIBUTE message via WebSocket
                // Example:
                // if (window.WebSocketClient && window.WebSocketClient.sendMessage) {
                //     const nodeId = target.closest('.tree-node')?.dataset.nodeId; // Simplistic way to get node context
                //     const attributeName = target.id; // This needs better mapping
                //     window.WebSocketClient.sendMessage({
                //         type: "SET_PLOT_ATTRIBUTE",
                //         payload: {
                //             targetId: nodeId || "unknown_target", // Refine this
                //             attribute: attributeName, // Refine this
                //             value: target.type === 'checkbox' ? target.checked : target.value
                //         }
                //     });
                // }
            }
        });
    }


    console.log("Frontend application initialized.");
});
