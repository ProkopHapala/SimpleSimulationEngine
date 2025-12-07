// Minimal UI binding
import { initSimulationState, applyCoils, rebuildPlasma } from './physics.js';
import { Logger } from '../common_js/Logger.js';

export function initUI(sim, renderer) {
    const coilSlider = document.getElementById('coil-slider');
    const coilVal = document.getElementById('coil-val');
    const pSlider = document.getElementById('p-slider');
    const pVal = document.getElementById('p-val');
    const dtInput = document.getElementById('dt');
    const resetBtn = document.getElementById('reset-btn');
    const pauseBtn = document.getElementById('pause-btn');
    const applyBtn = document.getElementById('apply-btn');
    const solveBtn = document.getElementById('solve-btn');
    const verbSel = document.getElementById('verb');
    const plasmaRadius = document.getElementById('plasma-radius');
    const plasmaNodes = document.getElementById('plasma-nodes');
    const scText = document.getElementById('sc-coils');
    const cageText = document.getElementById('cage-coils');
    const logEl = document.getElementById('log');

    if (typeof window !== 'undefined') {
        window.logger = window.logger || new Logger();
        window.logger.setContainer(logEl);
    }

    const updateLabels = () => {
        // Display in kA
        coilVal.textContent = (sim.params.scCurrent / 1000.0).toFixed(1) + ' kA';
        pVal.textContent = sim.params.plasmaPressure0.toFixed(0);
    };

    coilSlider.addEventListener('input', () => {
        // Slider value 0..100 maps to 0..100 kA
        sim.params.scCurrent = parseFloat(coilSlider.value) * 1000.0;
        updateLabels();
    });
    pSlider.addEventListener('input', () => {
        sim.params.plasmaPressure0 = parseFloat(pSlider.value) * 100.0;
        sim.P0 = sim.params.plasmaPressure0;
        updateLabels();
    });
    dtInput.addEventListener('change', () => {
        sim.params.dt = parseFloat(dtInput.value) || sim.params.dt;
    });
    const parseCoils = (text) => {
        const lines = text.split('\n').map(l => l.trim()).filter(l => l.length > 0);
        return lines.map(l => {
            const parts = l.split(/\s+/).map(Number).filter(v => !Number.isNaN(v));
            if (parts.length >= 3) return { r: parts[0], z: parts[1], I: parts[2] };
            if (parts.length >= 2) return { r: parts[0], z: 0, I: parts[1] };
            return { r: 1, z: 0, I: 0 };
        });
    };

    applyBtn.addEventListener('click', () => {
        const radius = parseFloat(plasmaRadius.value) || 0.6;
        const nNodes = parseInt(plasmaNodes.value) || 24;
        rebuildPlasma(sim, radius, nNodes);
        const scList = parseCoils(scText.value);
        const cageList = parseCoils(cageText.value);
        applyCoils(sim, scList, cageList);
        sim.runSolver = true;
        log('Applied config and queued solve.');
    });

    solveBtn.addEventListener('click', () => {
        sim.runSolver = true;
        log('Solve requested.');
    });

    resetBtn.addEventListener('click', () => {
        const fresh = initSimulationState();
        Object.keys(sim).forEach((k) => delete sim[k]);
        Object.assign(sim, fresh);
        updateLabels();
        log('Reset simulation.');
    });
    pauseBtn.addEventListener('click', () => {
        sim.paused = !sim.paused;
        pauseBtn.textContent = sim.paused ? 'Resume' : 'Pause';
    });
    verbSel.addEventListener('change', () => {
        sim.params.solverVerb = parseInt(verbSel.value) || 0;
        log(`Verbosity set to ${sim.params.solverVerb}`);
    });

    function log(msg) {
        if (!logEl) return;
        const line = document.createElement('div');
        line.textContent = msg;
        logEl.prepend(line);
        while (logEl.childElementCount > 50) logEl.removeChild(logEl.lastChild);
        if (window.logger) window.logger.info(msg);
    }

    updateLabels();
    return { log };
}
