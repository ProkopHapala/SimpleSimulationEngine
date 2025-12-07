// UI bindings for MHD demo with side panel
import { initSimulationState, initTwoStates, switchState, sampleFieldAtPoint, stepSimulation, solveFluxConservation } from './physics.js';
import { Logger } from '../common_js/Logger.js';

export function initUI(sim, renderer) {
    // State buttons
    const btnState0 = document.getElementById('btn-state0');
    const btnState1 = document.getElementById('btn-state1');

    // Plasma parameters
    const plasmaR0 = document.getElementById('plasma-r0');
    const plasmaR1 = document.getElementById('plasma-r1');
    const plasmaZ0 = document.getElementById('plasma-z0');

    // SC parameters
    const scCurrent = document.getElementById('sc-current');
    const scRadius = document.getElementById('sc-radius');
    const scZ = document.getElementById('sc-z');

    // Cage parameters
    const cageF = document.getElementById('cage-f');
    const cageZMin = document.getElementById('cage-zmin');
    const cageZMax = document.getElementById('cage-zmax');
    const cageCount = document.getElementById('cage-count');

    // Solver
    const solverMethod = document.getElementById('solver-method');
    const btnSolve = document.getElementById('btn-solve');
    const btnReset = document.getElementById('btn-reset');

    // Visualization
    const chkField = document.getElementById('chk-field');
    const chkControlPts = document.getElementById('chk-control-pts');
    const chkRings = document.getElementById('chk-rings');
    const chkProfile = document.getElementById('chk-profile');
    const fieldScale = document.getElementById('field-scale');

    const logEl = document.getElementById('log');

    // Setup logger
    if (typeof window !== 'undefined') {
        window.logger = window.logger || new Logger();
        window.logger.setContainer(logEl);
    }

    function log(msg) {
        if (!logEl) return;
        const line = document.createElement('div');
        line.textContent = `[${new Date().toLocaleTimeString()}] ${msg}`;
        logEl.prepend(line);
        while (logEl.childElementCount > 50) logEl.removeChild(logEl.lastChild);
        if (window.logger) window.logger.info(msg);
    }

    function getConfig() {
        return {
            plasmaR0: parseFloat(plasmaR0.value) || 0.1,
            plasmaR1: parseFloat(plasmaR1.value) || 0.5,
            plasmaZ0: parseFloat(plasmaZ0.value) || 0.25,
            scCurrent: (parseFloat(scCurrent.value) || 50) * 1000, // kA to A
            scRadius: parseFloat(scRadius.value) || 1.2,
            scZ: parseFloat(scZ.value) || -0.5,
            cageF: parseFloat(cageF.value) || 0.25,
            cageZMin: parseFloat(cageZMin.value) || 0.1,
            cageZMax: parseFloat(cageZMax.value) || 3.0,
            cageCount: parseInt(cageCount.value) || 8,
            plasmaCount: 8,
        };
    }

    function applyConfig() {
        const config = getConfig();
        initTwoStates(sim, config);
        updateStateButtons();
        log(`Applied config: plasma r0=${config.plasmaR0}, r1=${config.plasmaR1}, SC=${config.scCurrent / 1000}kA`);
    }

    function updateStateButtons() {
        btnState0.classList.toggle('active', sim.currentStateIndex === 0);
        btnState1.classList.toggle('active', sim.currentStateIndex === 1);
    }

    function updateRenderOptions() {
        if (renderer) {
            renderer.showField = chkField.checked;
            renderer.showControlPoints = chkControlPts.checked;
            renderer.showRings = chkRings.checked;
            renderer.showProfile = chkProfile.checked;
            renderer.fieldScale = parseFloat(fieldScale.value) || 1.0;
        }
    }

    const verbSel = document.createElement('select');
    verbSel.innerHTML = `
        <option value="0">0 (None)</option>
        <option value="1" selected>1 (Basic)</option>
        <option value="2">2 (Detailed)</option>
        <option value="3">3 (Matrix)</option>
    `;
    // Add verbosity control
    const verbRow = document.createElement('div');
    verbRow.className = 'control-row';
    const verbLabel = document.createElement('label');
    verbLabel.textContent = 'Verbosity';
    verbRow.appendChild(verbLabel);
    verbRow.appendChild(verbSel);

    // Append after Visualization section? Or to Log section.
    // Let's put it in Solver section for now or Log.
    logEl.parentElement.insertBefore(verbRow, logEl);

    // ... existing initialization ...

    verbSel.addEventListener('change', () => {
        sim.params.solverVerb = parseInt(verbSel.value) || 0;
        log(`Verbosity set to ${sim.params.solverVerb}`);
    });

    // Solve button
    btnSolve.addEventListener('click', () => {
        const method = solverMethod.value;
        log(`Solving currents using ${method} method...`);

        let solved = false;
        if (method === 'flux') {
            solveFluxConservation(sim);
            log('Flux conservation solver completed');
            solved = true;
        } else {
            // Control points method - use stepSimulation
            sim.runSolver = true;
            stepSimulation(sim);
            log('Control point solver completed');
            solved = true;
        }

        if (solved && renderer) {
            // Force label update if renderer supports it
            if (renderer.updateLabels) renderer.updateLabels(sim);
        }
    });

    // Update labels when switching state if solution exists
    function onStateSwitch() {
        if (renderer && renderer.updateLabels) renderer.updateLabels(sim);
        // Maybe log currents?
    }

    // State toggle buttons
    btnState0.addEventListener('click', () => {
        switchState(sim, 0);
        updateStateButtons();
        onStateSwitch();
        log('Switched to t=0 (Initial state)');
    });

    btnState1.addEventListener('click', () => {
        switchState(sim, 1);
        updateStateButtons();
        onStateSwitch();
        log('Switched to t=1 (Expanded state)');
    });

    // Reset button
    btnReset.addEventListener('click', () => {
        applyConfig();
        if (renderer && renderer.updateLabels) renderer.updateLabels(sim);
        log('Reset simulation');
    });

    // Visualization checkboxes
    chkField.addEventListener('change', updateRenderOptions);
    chkControlPts.addEventListener('change', updateRenderOptions);
    chkRings.addEventListener('change', updateRenderOptions);
    chkProfile.addEventListener('change', updateRenderOptions);
    fieldScale.addEventListener('change', updateRenderOptions);

    // Parameter inputs - apply on change
    const paramInputs = [plasmaR0, plasmaR1, plasmaZ0, scCurrent, scRadius, scZ,
        cageF, cageZMin, cageZMax, cageCount];
    paramInputs.forEach(input => {
        input.addEventListener('change', () => {
            applyConfig();
            if (renderer && renderer.updateLabels) renderer.updateLabels(sim);
        });
    });

    // Initial setup
    applyConfig();
    updateRenderOptions();
    if (renderer && renderer.updateLabels) renderer.updateLabels(sim);
    log('MHD Demo initialized with two-state visualization');

    return { log, applyConfig, updateRenderOptions };
}
