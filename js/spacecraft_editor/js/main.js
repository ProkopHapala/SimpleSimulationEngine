import { SpaceCraftEngine } from './SpaceCraftEngine.js';
import { SpaceCraftRenderer } from './SpaceCraftRenderer.js';
import { GUI } from './GUI.js';

// Global instances (kept on window for debugging)
window.engine = null;
window.renderer = null;

document.addEventListener('DOMContentLoaded', () => {
    console.log("Initializing Spacecraft Editor...");

    // NOTE: Global logger is initialized in js/common_js/Logger.js which should be loaded before this script.
    // We ensure it is properly set up with the UI container below.


    // Shader Loading Helper
    const loadShader = async (url) => {
        const response = await fetch(url);
        return await response.text();
    };

    // Load all shaders
    Promise.all([
        loadShader('../common_resources/shaders/atom.glslv'),
        loadShader('../common_resources/shaders/atom.glslf'),
        loadShader('../common_resources/shaders/bond.glslv'),
        loadShader('../common_resources/shaders/bond_color.glslf'), // bond fragment uses palette-based colors
        loadShader('../common_resources/shaders/label.glslv'),
        loadShader('../common_resources/shaders/label.glslf')
    ]).then(([atomVert, atomFrag, bondVert, bondFrag, labelVert, labelFrag]) => {
        console.log("Shaders loaded.");

        const shaders = {
            atom: { vertex: atomVert, fragment: atomFrag },
            bond: { vertex: bondVert, fragment: bondFrag },
            label: { vertex: labelVert, fragment: labelFrag },
            // Map for SpaceCraftRenderer (it expects specific names if we didn't change it fully, 
            // but we updated SpaceCraftRenderer to use MeshRenderer which uses atom/bond/label keys if we pass them correctly)
            // Wait, MeshRenderer expects shaders.atom, shaders.bond, shaders.label.
            // SpaceCraftRenderer init(container, shaders) passes this object.
            // So we should structure it as MeshRenderer expects.

            // But SpaceCraftRenderer.initMeshes (old) used nodeVertex etc.
            // We updated SpaceCraftRenderer to call MeshRenderer.init which uses shaders.atom etc.
            // So this structure is correct for MeshRenderer.

            // However, SpaceCraftRenderer.initMeshes also used to look for nodeVertex/Fragment.
            // But we replaced initMeshes with MeshRenderer.init.
            // So we are good.
        };

        // 1. Initialize Engine
        window.engine = new SpaceCraftEngine();

        // 2. Initialize Renderer
        window.renderer = new SpaceCraftRenderer(window.engine);
        window.renderer.init(document.getElementById('canvas-container'), shaders);

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
            window.renderer.update();
            window.renderer.render();
        }
        animate();

    }).catch(err => {
        console.error("Failed to load shaders:", err);
    });
});

// Helper for logging to UI
// Deprecated: Use window.logger instead
window.logToUI = function (msg) {
    window.logger.info(msg);
};
