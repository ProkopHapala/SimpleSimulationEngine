import { initSimulationState, stepSimulation, initTwoStates } from './physics.js';
import { DemoRenderer } from './render.js';
import { initUI } from './ui.js';

window.mhdSim = null;
window.mhdRenderer = null;
window.mhdUI = null;

function init() {
    const container = document.getElementById('canvas-container');
    if (!container) {
        console.error('Missing #canvas-container in DOM');
        return;
    }

    // Initialize simulation state
    window.mhdSim = initSimulationState();

    // Create renderer
    window.mhdRenderer = new DemoRenderer(container);

    // Initialize UI (this will also call initTwoStates with default config)
    window.mhdUI = initUI(window.mhdSim, window.mhdRenderer);

    let lastT = performance.now();
    function animate(now) {
        requestAnimationFrame(animate);
        const dtVis = (now - lastT) * 0.001;
        lastT = now;

        // Run solver if requested
        if (window.mhdSim.runSolver) {
            stepSimulation(window.mhdSim);
            window.mhdSim.runSolver = false;
        }

        // Update controls
        if (window.mhdRenderer && window.mhdRenderer.controls) {
            window.mhdRenderer.controls.update();
        }

        // Render
        window.mhdRenderer.render(window.mhdSim, dtVis);
    }
    requestAnimationFrame(animate);
}

document.addEventListener('DOMContentLoaded', init);
