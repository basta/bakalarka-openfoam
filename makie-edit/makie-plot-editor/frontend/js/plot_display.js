// frontend/js/plot_display.js

function updatePlot(payload) {
    const plotDisplayArea = document.getElementById('plot-display');
    if (!plotDisplayArea) {
        console.error("Plot display area not found.");
        return;
    }

    if (payload.updateType === "IMAGE") {
        // Assuming payload.data is a base64 encoded image string
        // e.g., "data:image/png;base64,iVBORw0KGgo..."
        plotDisplayArea.innerHTML = `<img src="${payload.data}" alt="Makie Plot" style="max-width: 100%; max-height: 100%; object-fit: contain;">`;
        console.log("Plot display updated with image.");
    } else if (payload.updateType === "WGLMAKIE_SCENE") {
        // This part is highly dependent on how JSServe/WGLMakie sends scene updates.
        // For now, it's a placeholder. You'd likely interact with a JSServe-managed canvas.
        plotDisplayArea.innerHTML = `<p>WGLMakie scene update received (actual rendering TBD).</p><pre>${JSON.stringify(payload.sceneModifications, null, 2)}</pre>`;
        console.log("WGLMakie scene update received.");
        // Example: If JSServe creates a canvas with a specific ID, you might interact with it here
        // or JSServe's own client-side JS would handle it.
    } else {
        console.warn("Unknown plot update type:", payload.updateType);
        plotDisplayArea.textContent = `Received plot update of unknown type: ${payload.updateType}`;
    }
}

// Expose functions
window.PlotDisplayModule = {
    updatePlot
};
