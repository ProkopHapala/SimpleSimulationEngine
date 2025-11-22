
// Global instances
window.engine = null;
window.renderer = null;

document.addEventListener('DOMContentLoaded', () => {
    console.log("Initializing Spacecraft Editor...");

    // 1. Initialize Engine
    window.engine = new SpaceCraftEngine();

    // 2. Initialize Renderer
    window.renderer = new SpaceCraftRenderer(window.engine);
    window.renderer.init(document.getElementById('canvas-container'));

    // 3. Initialize GUI
    window.gui = new GUI(window.engine, window.renderer);

    // 4. Start Loop
    // --- UI Elements ---
    const logElement = document.getElementById('log');

    // Init Logger
    if (window.logger) {
        window.logger.setContainer(logElement);
        window.logger.info("Logger initialized.");
    }

    // --- Animation Loop ---
    function animate() {
        requestAnimationFrame(animate);

        // Update Engine (Physics/Logic)
        // window.engine.update(); // TODO: Physics

        // Update Renderer
        window.renderer.update();
        window.renderer.render();
    }
    animate();
});

// Helper for logging to UI
// Deprecated: Use window.logger instead
window.logToUI = function (msg) {
    window.logger.info(msg);
};
