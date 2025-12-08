// UI bindings for MHD demo with side panel
import { initTwoStates, switchState, stepSimulation, solveFluxConservation, initFromCoilList } from './physics.js';
import { Logger } from '../common_js/Logger.js';

export function initUI(sim, renderer) {
    // Tabs
    const tabSettings = document.getElementById('tab-settings');
    const tabEditor = document.getElementById('tab-editor');
    const panelSettings = document.getElementById('panel-settings');
    const panelEditor = document.getElementById('panel-editor');

    // Editor
    const configTextT0 = document.getElementById('config-text-t0');
    const configTextT1 = document.getElementById('config-text-t1');
    const btnGenConfig = document.getElementById('btn-gen-config');
    const btnLoadConfig = document.getElementById('btn-load-config');

    // Settings Inputs
    const btnState0 = document.getElementById('btn-state0');
    const btnState1 = document.getElementById('btn-state1');

    const scCurrent = document.getElementById('sc-current');
    const scRadius = document.getElementById('sc-radius');
    const scZ = document.getElementById('sc-z');

    // Cage
    const cageCount = document.getElementById('cage-count');
    const cageRMin = document.getElementById('cage-r-min');
    const cageRMax = document.getElementById('cage-r-max');
    const cageZStart = document.getElementById('cage-z-start');
    const cageZEnd = document.getElementById('cage-z-end');

    // Plasma
    const plasmaCount = document.getElementById('plasma-count');
    const plasmaZ0 = document.getElementById('plasma-z0');
    const plasmaR0 = document.getElementById('plasma-r0');
    const plasmaR1 = document.getElementById('plasma-r1');

    // Trajectory
    const tSlider = document.getElementById('t-slider');
    const tValue = document.getElementById('t-value');

    // Solver & Viz
    const solverMethod = document.getElementById('solver-method');
    const verbosity = document.getElementById('verbosity');
    const btnSolve = document.getElementById('btn-solve');
    const btnReset = document.getElementById('btn-reset');
    const chkAutoSolve = document.getElementById('chk-auto-solve');

    const chkField = document.getElementById('chk-field');
    const chkControlPts = document.getElementById('chk-control-pts');
    const chkRings = document.getElementById('chk-rings');
    const chkShader = document.getElementById('chk-shader');
    const fieldScale = document.getElementById('field-scale');
    const bMax = document.getElementById('b-max');

    const energyPlot = document.getElementById('energy-plot');
    const logEl = document.getElementById('log');

    // Trajectory data (pre-computed)
    let trajectory = [];
    let currentTIndex = 0;

    // Logger Setup
    if (typeof window !== 'undefined') {
        window.logger = window.logger || new Logger();
        if (logEl) window.logger.setContainer(logEl);
        sim.params = sim.params || {};
        sim.params.solverVerb = 1;
    }

    function log(msg) {
        if (!logEl) return;
        const line = document.createElement('div');
        line.textContent = `[${new Date().toLocaleTimeString()}] ${msg}`;
        logEl.prepend(line);
        while (logEl.childElementCount > 50) logEl.removeChild(logEl.lastChild);
    }

    // Generate config from settings - EQUIDISTANT R for cage!
    function generateConfigText() {
        const lines = ["# TYPE  R0     Z0     R1     Z1     I0"];

        const sc_r = parseFloat(scRadius?.value) || 0.3;
        const sc_z = parseFloat(scZ?.value) || 0.3;
        const sc_I = parseFloat(scCurrent?.value) || 1e6;
        lines.push(`SC      ${sc_r.toFixed(3)}  ${sc_z.toFixed(3)}   ${sc_r.toFixed(3)}  ${sc_z.toFixed(3)}   ${sc_I.toExponential(1)}`);

        const nCage = parseInt(cageCount?.value) || 5;
        const rMin = parseFloat(cageRMin?.value) || 0.2;
        const rMax = parseFloat(cageRMax?.value) || 2.0;
        const zStart = parseFloat(cageZStart?.value) || 0.5;
        const zEnd = parseFloat(cageZEnd?.value) || 1.8;

        // EQUIDISTANT R: z = z0 + A * (r^2 - rMin^2)
        // At rMin: z = zStart
        // At rMax: z = zEnd
        // A = (zEnd - zStart) / (rMax^2 - rMin^2)
        const A = (zEnd - zStart) / Math.max(0.001, rMax * rMax - rMin * rMin);

        for (let i = 0; i < nCage; i++) {
            const t = nCage > 1 ? i / (nCage - 1) : 0;
            const r = rMin + (rMax - rMin) * t;  // Equidistant in R
            const z = zStart + A * (r * r - rMin * rMin);
            lines.push(`CAGE    ${r.toFixed(3)}  ${z.toFixed(3)}   ${r.toFixed(3)}  ${z.toFixed(3)}   0.0`);
        }

        const nPlasma = parseInt(plasmaCount?.value) || 1;
        const pZ = parseFloat(plasmaZ0?.value) || 1.0;
        const pR0 = parseFloat(plasmaR0?.value) || 0.01;
        const pR1 = parseFloat(plasmaR1?.value) || 1.0;

        for (let i = 0; i < nPlasma; i++) {
            lines.push(`PLASMA  ${pR0.toFixed(3)}  ${pZ.toFixed(3)}   ${pR1.toFixed(3)}  ${pZ.toFixed(3)}   0.0`);
        }

        return lines.join('\n');
    }

    function parseConfigText(text) {
        const coils = [];
        for (let line of text.split('\n')) {
            line = line.trim();
            if (!line || line.startsWith('#')) continue;
            const parts = line.split(/\s+/);
            if (parts.length < 6) continue;
            coils.push({
                type: parts[0],
                R0: parseFloat(parts[1]), Z0: parseFloat(parts[2]),
                R1: parseFloat(parts[3]), Z1: parseFloat(parts[4]),
                I0: parseFloat(parts[5])
            });
        }
        return coils;
    }

    function loadConfigFromTextarea() {
        const text = configTextT0?.value || '';
        const coils = parseConfigText(text);
        if (coils.length > 0) {
            initFromCoilList(sim, coils);
            trajectory = [];
            if (renderer?.updateLabels) renderer.updateLabels(sim);
        }
    }

    // Interpolate plasma between t=0 and t=1
    function setStateAtT(t) {
        if (!sim.state0 || !sim.state1) return;

        const p0 = sim.state0.plasmaNodes;
        const p1 = sim.state1.plasmaNodes;

        sim.plasmaNodes = p0.map((n0, i) => {
            const n1 = p1[i] || n0;
            return {
                r: n0.r * (1 - t) + n1.r * t,
                z: n0.z * (1 - t) + n1.z * t,
                vr: 0, vz: 0, m: 1.0
            };
        });

        if (trajectory.length > 0) {
            const idx = Math.min(trajectory.length - 1, Math.round(t * (trajectory.length - 1)));
            const step = trajectory[idx];
            sim.plasmaCurrents = step.plasmaCurrents.slice();
            sim.cageCurrents = step.cageCurrents.slice();
            currentTIndex = idx;
        } else {
            const c0 = sim.state0.plasmaCurrents || [];
            const c1 = sim.state1.plasmaCurrents || [];
            sim.plasmaCurrents = c0.map((i0, idx) => {
                const i1 = c1[idx] || 0;
                return i0 * (1 - t) + i1 * t;
            });
        }

        sim.currentT = t;
    }

    function updateStateButtons() {
        const t = sim.currentT || 0;
        if (btnState0) btnState0.classList.toggle('active', t < 0.5);
        if (btnState1) btnState1.classList.toggle('active', t >= 0.5);
    }

    function updateRenderOptions() {
        if (renderer) {
            renderer.showField = chkField?.checked ?? true;
            renderer.showControlPoints = chkControlPts?.checked ?? false;
            renderer.showRings = chkRings?.checked ?? true;
            renderer.showShader = chkShader?.checked ?? false;
            renderer.fieldScale = parseFloat(fieldScale?.value) || 1.0;
            renderer.bMax = parseFloat(bMax?.value) || 0.1;
        }
    }

    function drawEnergyPlot() {
        if (!energyPlot) return;
        const ctx = energyPlot.getContext('2d');
        const w = energyPlot.width;
        const h = energyPlot.height;

        ctx.fillStyle = '#111';
        ctx.fillRect(0, 0, w, h);

        if (trajectory.length < 2) {
            ctx.fillStyle = '#555';
            ctx.font = '12px sans-serif';
            ctx.fillText('Click "Solve" to compute', 10, h / 2);
            return;
        }

        const maxE = Math.max(...trajectory.map(e => Math.abs(e.energy))) || 1;
        const maxI = Math.max(...trajectory.map(e => Math.abs(e.totalCurrent))) || 1;

        ctx.strokeStyle = '#0f0';
        ctx.lineWidth = 2;
        ctx.beginPath();
        trajectory.forEach((step, i) => {
            const x = (i / (trajectory.length - 1)) * w;
            const y = h - (Math.abs(step.energy) / maxE) * h * 0.85 - 5;
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        });
        ctx.stroke();

        ctx.strokeStyle = '#f80';
        ctx.beginPath();
        trajectory.forEach((step, i) => {
            const x = (i / (trajectory.length - 1)) * w;
            const y = h - (Math.abs(step.totalCurrent) / maxI) * h * 0.85 - 5;
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        });
        ctx.stroke();

        if (currentTIndex >= 0 && currentTIndex < trajectory.length) {
            const x = (currentTIndex / (trajectory.length - 1)) * w;
            ctx.strokeStyle = '#fff';
            ctx.lineWidth = 1;
            ctx.beginPath();
            ctx.moveTo(x, 0);
            ctx.lineTo(x, h);
            ctx.stroke();
        }

        ctx.font = '10px sans-serif';
        ctx.fillStyle = '#0f0';
        ctx.fillText('Energy', 5, 12);
        ctx.fillStyle = '#f80';
        ctx.fillText('Current', 50, 12);
    }

    function computeTrajectory(nSteps = 20) {
        trajectory = [];

        for (let i = 0; i <= nSteps; i++) {
            const t = i / nSteps;
            setStateAtT(t);
            solveFluxConservation(sim);

            const plasmaCurrents = sim.plasmaCurrents.slice();
            const cageCurrents = sim.cageCurrents.slice();
            const plasmaNodes = sim.plasmaNodes.map(n => ({ ...n }));

            let totalCurrent = 0;
            let energy = 0;
            plasmaCurrents.forEach(I => {
                totalCurrent += I;
                energy += 0.5 * I * I * 1e-9;
            });
            cageCurrents.forEach(I => {
                energy += 0.5 * I * I * 1e-9;
            });

            trajectory.push({ t, plasmaNodes, plasmaCurrents, cageCurrents, totalCurrent, energy });
        }

        drawEnergyPlot();
        log(`Trajectory computed (${nSteps + 1} steps)`);
    }

    // Auto-solve if checkbox is checked
    function autoSolveIfEnabled() {
        if (chkAutoSolve?.checked) {
            computeTrajectory(20);
            const t = parseFloat(tSlider?.value) || 0;
            setStateAtT(t);
            drawEnergyPlot();
        }
    }

    function addWheelHandler(el) {
        if (!el) return;
        el.addEventListener('wheel', (e) => {
            e.preventDefault();
            const step = parseFloat(el.step) || 0.1;
            const delta = e.deltaY < 0 ? step : -step;
            el.value = (parseFloat(el.value) + delta).toFixed(4);
            el.dispatchEvent(new Event('change'));
        }, { passive: false });
    }

    // Tabs
    if (tabSettings) tabSettings.addEventListener('click', () => {
        tabSettings.classList.add('active');
        tabEditor?.classList.remove('active');
        if (panelSettings) panelSettings.style.display = 'block';
        if (panelEditor) panelEditor.style.display = 'none';
    });

    if (tabEditor) tabEditor.addEventListener('click', () => {
        tabSettings?.classList.remove('active');
        tabEditor.classList.add('active');
        if (panelSettings) panelSettings.style.display = 'none';
        if (panelEditor) panelEditor.style.display = 'block';
    });

    // Editor Buttons
    if (btnGenConfig) btnGenConfig.addEventListener('click', () => {
        if (configTextT0) configTextT0.value = generateConfigText();
        log("Generated config");
    });

    if (btnLoadConfig) btnLoadConfig.addEventListener('click', () => {
        loadConfigFromTextarea();
        autoSolveIfEnabled();
        log("Loaded config");
    });

    // State Buttons
    if (btnState0) btnState0.addEventListener('click', () => {
        if (tSlider) tSlider.value = 0;
        if (tValue) tValue.textContent = '0.00';
        setStateAtT(0);
        updateStateButtons();
        drawEnergyPlot();
        if (renderer?.updateLabels) renderer.updateLabels(sim);
    });

    if (btnState1) btnState1.addEventListener('click', () => {
        if (tSlider) tSlider.value = 1;
        if (tValue) tValue.textContent = '1.00';
        setStateAtT(1);
        updateStateButtons();
        drawEnergyPlot();
        if (renderer?.updateLabels) renderer.updateLabels(sim);
    });

    // T-slider
    if (tSlider) tSlider.addEventListener('input', () => {
        const t = parseFloat(tSlider.value);
        if (tValue) tValue.textContent = t.toFixed(2);
        setStateAtT(t);
        updateStateButtons();
        drawEnergyPlot();
        if (renderer?.updateLabels) renderer.updateLabels(sim);
    });

    // Solver
    if (btnSolve) btnSolve.addEventListener('click', () => {
        computeTrajectory(20);
        const t = parseFloat(tSlider?.value) || 0;
        setStateAtT(t);
        drawEnergyPlot();
        if (renderer?.updateLabels) renderer.updateLabels(sim);
    });

    if (btnReset) btnReset.addEventListener('click', () => {
        if (configTextT0) configTextT0.value = generateConfigText();
        loadConfigFromTextarea();
        trajectory = [];
        currentTIndex = 0;
        if (tSlider) tSlider.value = 0;
        if (tValue) tValue.textContent = '0.00';
        autoSolveIfEnabled();
        drawEnergyPlot();
        log("Reset");
    });

    if (verbosity) verbosity.addEventListener('change', () => {
        sim.params.solverVerb = parseInt(verbosity.value) || 0;
    });

    // Viz checkboxes
    [chkField, chkControlPts, chkRings, chkShader, fieldScale, bMax].forEach(el => {
        if (el) el.addEventListener('change', updateRenderOptions);
    });

    // Param inputs - auto-solve on change
    const allInputs = [plasmaR0, plasmaR1, plasmaZ0, plasmaCount, scCurrent, scRadius, scZ,
        cageRMin, cageRMax, cageZStart, cageZEnd, cageCount, fieldScale, bMax];

    allInputs.forEach(el => {
        addWheelHandler(el);
        if (el) el.addEventListener('change', () => {
            if (configTextT0) configTextT0.value = generateConfigText();
            loadConfigFromTextarea();
            autoSolveIfEnabled();
            drawEnergyPlot();
        });
    });

    // INIT: Generate config, load, and auto-solve on startup
    if (configTextT0) configTextT0.value = generateConfigText();
    loadConfigFromTextarea();
    updateRenderOptions();
    autoSolveIfEnabled();  // Run at startup!
    drawEnergyPlot();
    log('MHD Demo ready');

    return { log, computeTrajectory, setStateAtT };
}
