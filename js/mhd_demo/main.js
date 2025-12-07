import { initSimulationState, stepSimulation } from './physics.js';
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

    window.mhdSim = initSimulationState();
    window.mhdRenderer = new DemoRenderer(container);
    window.mhdUI = initUI(window.mhdSim, window.mhdRenderer);

    let lastT = performance.now();
    function animate(now) {
        requestAnimationFrame(animate);
        const dtVis = (now - lastT) * 0.001;
        lastT = now;
        if (!window.mhdSim.paused && window.mhdSim.autoStep) {
            const sub = 2;
            for (let i = 0; i < sub; i++) stepSimulation(window.mhdSim);
        }
        if (window.mhdRenderer && window.mhdRenderer.controls) {
            window.mhdRenderer.controls.update();
        }
        window.mhdRenderer.render(window.mhdSim, dtVis);
    }
    requestAnimationFrame(animate);
}

document.addEventListener('DOMContentLoaded', init);
