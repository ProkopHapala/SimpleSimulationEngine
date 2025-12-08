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
    const trajModeSel = document.getElementById('traj-mode');
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
    let trajMode = 'interp'; // 'interp' or 'dynamic'

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

    function cloneNodes(nodes) {
        return nodes.map(n => ({ ...n }));
    }

    // Set state by interpolating or sampling precomputed trajectory
    function setStateAtT(t) {
        if (trajMode === 'dynamic' && trajectory.length > 0) {
            const idx = Math.min(trajectory.length - 1, Math.round(t * (trajectory.length - 1)));
            const step = trajectory[idx];
            sim.plasmaNodes = cloneNodes(step.plasmaNodes);
            sim.plasmaCurrents = step.plasmaCurrents.slice();
            sim.cageCurrents = step.cageCurrents.slice();
            currentTIndex = idx;
            sim.currentT = step.t;
            return;
        }

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
        trajMode = trajModeSel?.value || 'interp';
        sim.trajMode = trajMode;

        if (trajMode === 'dynamic') {
            // ============================================================
            // Dynamic simulation following Python demo_plasma_dynamics_flux.py
            // Single plasma coil as point mass, flux-conserving currents,
            // Lorentz + gas pressure forces, Euler integration
            // ============================================================
            
            // Parameters matching Python defaults
            const dt = 5e-6;           // time step [s]
            const mass = 0.1;          // plasma coil mass [kg]
            const P0 = 1e6;            // initial gas pressure [Pa]
            const gamma_gas = 5.0/3.0; // adiabatic index
            const L_eff = 1.0;         // effective axial length for pressure model
            const r_min = 0.01;        // minimum radius guard
            const nDynSteps = 400;     // integration steps (like Python default)
            
            // Get initial plasma geometry from state0 (same as interpolated t=0)
            const p0 = sim.state0?.plasmaNodes || sim.plasmaNodes;
            if (!p0 || p0.length === 0) {
                log('No plasma nodes for dynamic simulation');
                return;
            }
            
            // Use middle node as the single "plasma coil" (like Python's single plasma loop)
            const midIdx = Math.floor(p0.length / 2);
            let r_p = p0[midIdx].r;
            let z_p = p0[midIdx].z;
            
            // Initial velocity (outward radial kick, like Python v0 = [1000, 0])
            let vr = sim.params?.initVr ?? 1000.0;
            let vz = sim.params?.initVz ?? 0.0;
            
            // Initial volume for gas pressure model: V = pi * r^2 * L_eff
            const r0 = r_p;
            const V0_gas = Math.PI * r0 * r0 * L_eff;
            
            // Build initial flux Phi0 from SC coil through all loops at t=0
            // For simplicity, we track flux through the plasma coil only
            const MU0 = 4e-7 * Math.PI;
            
            // Helper: mutual inductance between two loops
            const mutualInd = (r1, z1, r2, z2) => {
                const eps = 1e-9;
                const rp1 = Math.max(eps, r1), rp2 = Math.max(eps, r2);
                const dz = z1 - z2;
                const num = 4 * rp1 * rp2;
                const den = (rp1 + rp2) * (rp1 + rp2) + dz * dz + eps;
                let k2 = Math.min(1 - 1e-9, Math.max(0, num / den));
                if (k2 < 1e-12) return 0;
                const k = Math.sqrt(k2);
                // Approximate elliptic integrals
                let a = 1.0, b = Math.sqrt(1.0 - k2);
                for (let i = 0; i < 15; i++) { const an = 0.5*(a+b); b = Math.sqrt(a*b); a = an; }
                const K = Math.PI / (2.0 * a);
                a = 1.0; b = Math.sqrt(1.0 - k2); let sum = 0, twoPow = 1;
                for (let i = 0; i < 15; i++) { const an = 0.5*(a+b), cn = 0.5*(a-b); sum += twoPow*cn*cn; twoPow *= 2; b = Math.sqrt(a*b); a = an; }
                const E = Math.PI * 0.25 * (2.0*a*a - sum) / a;
                const M = MU0 * Math.sqrt(rp1*rp2) * ((2.0/k - k)*K - (2.0/k)*E);
                return isFinite(M) ? M : 0;
            };
            
            // Self inductance
            const selfInd = (r) => {
                const rp = Math.max(1e-9, r), rw = 0.01;
                return MU0 * rp * (Math.log(8*rp/rw) - 2);
            };
            
            // B-field from a loop at (a, zc) with current I, evaluated at (rp, zp)
            const fieldFromLoop = (a, zc, I, rp, zp) => {
                const eps = 1e-9;
                const rSafe = Math.max(eps, rp);
                const dz = zp - zc;
                const num = 4 * a * rSafe;
                const den = (a + rSafe)*(a + rSafe) + dz*dz + eps;
                let k2 = Math.min(1 - 1e-9, Math.max(0, num / den));
                let a_ = 1.0, b_ = Math.sqrt(1.0 - k2);
                for (let i = 0; i < 15; i++) { const an = 0.5*(a_+b_); b_ = Math.sqrt(a_*b_); a_ = an; }
                const K = Math.PI / (2.0 * a_);
                a_ = 1.0; b_ = Math.sqrt(1.0 - k2); let sum = 0, twoPow = 1;
                for (let i = 0; i < 15; i++) { const an = 0.5*(a_+b_), cn = 0.5*(a_-b_); sum += twoPow*cn*cn; twoPow *= 2; b_ = Math.sqrt(a_*b_); a_ = an; }
                const E = Math.PI * 0.25 * (2.0*a_*a_ - sum) / a_;
                const denom = Math.sqrt(den);
                const common = MU0 * I / (2 * Math.PI * (denom + eps));
                const denom2 = (a - rSafe)*(a - rSafe) + dz*dz + eps;
                let Br = common * (dz / rSafe) * (-K + (a*a + rSafe*rSafe + dz*dz) / denom2 * E);
                let Bz = common * (K + (a*a - rSafe*rSafe - dz*dz) / denom2 * E);
                if (rp < eps) Br = 0;
                return { Br: isFinite(Br) ? Br : 0, Bz: isFinite(Bz) ? Bz : 0 };
            };
            
            // SC coil info
            const sc = sim.scCoils[0] || { r: 1.0, z: 0.0, I: 1e6 };
            
            // Cage coils (fixed)
            const cage = sim.cageCoils || [];
            const nCage = cage.length;
            
            // Total loops = cage + 1 plasma
            const nLoops = nCage + 1;
            
            // Initial flux Phi0 through each loop from SC only (at t=0 positions)
            const Phi0 = [];
            for (let i = 0; i < nCage; i++) {
                Phi0.push(mutualInd(cage[i].r, cage[i].z, sc.r, sc.z) * sc.I);
            }
            Phi0.push(mutualInd(r_p, z_p, sc.r, sc.z) * sc.I); // plasma at initial position
            
            // Gaussian solve helper
            const gaussSolve = (A, b) => {
                const n = b.length;
                const aug = A.map((row, i) => [...row, b[i]]);
                for (let col = 0; col < n; col++) {
                    let maxRow = col, maxVal = Math.abs(aug[col][col]);
                    for (let row = col+1; row < n; row++) if (Math.abs(aug[row][col]) > maxVal) { maxVal = Math.abs(aug[row][col]); maxRow = row; }
                    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
                    const pivot = aug[col][col];
                    if (Math.abs(pivot) < 1e-12) continue;
                    for (let row = col+1; row < n; row++) { const f = aug[row][col]/pivot; for (let j = col; j <= n; j++) aug[row][j] -= f*aug[col][j]; }
                }
                const x = new Array(n).fill(0);
                for (let i = n-1; i >= 0; i--) { let s = aug[i][n]; for (let j = i+1; j < n; j++) s -= aug[i][j]*x[j]; if (Math.abs(aug[i][i]) > 1e-12) x[i] = s/aug[i][i]; }
                return x;
            };
            
            // Snapshot helper
            const pushSnapshot = (tVal, r_cur, z_cur, I_cage, I_plasma) => {
                // Reconstruct plasmaNodes with current plasma position
                const plasmaNodes = cloneNodes(p0);
                plasmaNodes[midIdx].r = r_cur;
                plasmaNodes[midIdx].z = z_cur;
                const plasmaCurrents = sim.plasmaCurrents.map(() => 0);
                plasmaCurrents[midIdx] = I_plasma;
                const cageCurrents = I_cage.slice();
                let totalCurrent = I_plasma;
                let energy = 0.5 * I_plasma * I_plasma * 1e-9;
                cageCurrents.forEach(I => { energy += 0.5 * I * I * 1e-9; });
                
                if (typeof window !== 'undefined' && window.logger) {
                    window.logger.info(`dyn t=${tVal.toFixed(3)} r=${r_cur.toFixed(4)} z=${z_cur.toFixed(4)} I_p=${(I_plasma/1e6).toFixed(3)}MA`, window.logger.INFO);
                }
                
                trajectory.push({ t: tVal, plasmaNodes, plasmaCurrents, cageCurrents, totalCurrent, energy });
            };
            
            // Time integration loop (like Python run_plasma_dynamics)
            for (let step = 0; step <= nDynSteps; step++) {
                // 1) Build inductance matrix K at current positions
                const K = [];
                for (let i = 0; i < nLoops; i++) {
                    const row = [];
                    for (let j = 0; j < nLoops; j++) {
                        const ri = i < nCage ? cage[i].r : r_p;
                        const zi = i < nCage ? cage[i].z : z_p;
                        const rj = j < nCage ? cage[j].r : r_p;
                        const zj = j < nCage ? cage[j].z : z_p;
                        row.push(i === j ? selfInd(ri) : mutualInd(ri, zi, rj, zj));
                    }
                    K.push(row);
                }
                
                // 2) Solve K @ I = Phi0 for currents
                const I_all = gaussSolve(K, Phi0.slice());
                const I_cage = I_all.slice(0, nCage);
                const I_plasma = I_all[nCage] || 0;
                
                // 3) Store snapshot (map step to t in [0,1])
                const tVal = step / nDynSteps;
                pushSnapshot(tVal, r_p, z_p, I_cage, I_plasma);
                
                if (step === nDynSteps) break; // last step, don't integrate further
                
                // 4) Compute Lorentz force on plasma from SC and cage (not self)
                let Br_tot = 0, Bz_tot = 0;
                // From SC
                const Bsc = fieldFromLoop(sc.r, sc.z, sc.I, r_p, z_p);
                Br_tot += Bsc.Br; Bz_tot += Bsc.Bz;
                // From cage coils
                for (let i = 0; i < nCage; i++) {
                    const Bc = fieldFromLoop(cage[i].r, cage[i].z, I_cage[i], r_p, z_p);
                    Br_tot += Bc.Br; Bz_tot += Bc.Bz;
                }
                const L_len = 2 * Math.PI * r_p;
                const Fr_lorentz = I_plasma * L_len * Bz_tot;
                const Fz_lorentz = I_plasma * L_len * (-Br_tot);
                
                // 5) Gas pressure force (adiabatic P = P0 * (V0/V)^gamma)
                const V_cur = Math.PI * Math.max(r_p, r_min) * Math.max(r_p, r_min) * L_eff;
                const P_cur = P0 * Math.pow(V0_gas / V_cur, gamma_gas);
                const Fr_gas = P_cur * (2 * Math.PI * r_p * L_eff);
                
                // Total force
                const Fr = Fr_lorentz + Fr_gas;
                const Fz = Fz_lorentz;
                
                // 6) Euler integration
                const ar = Fr / mass;
                const az = Fz / mass;
                vr += ar * dt;
                vz += az * dt;
                r_p += vr * dt;
                z_p += vz * dt;
                
                // Guard against collapse through axis
                if (r_p < r_min) { r_p = r_min; vr = -0.5 * vr; }
            }
            
            log(`Dynamic trajectory computed (${nDynSteps + 1} steps)`);
        } else {
            for (let i = 0; i <= nSteps; i++) {
                const t = i / nSteps;
                setStateAtT(t);
                solveFluxConservation(sim);

                const plasmaCurrents = sim.plasmaCurrents.slice();
                const cageCurrents = sim.cageCurrents.slice();
                const plasmaNodes = cloneNodes(sim.plasmaNodes);

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
            log(`Interpolated trajectory computed (${nSteps + 1} steps)`);
        }

        drawEnergyPlot();
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

    // Trajectory mode
    if (trajModeSel) {
        trajModeSel.addEventListener('change', () => {
            trajMode = trajModeSel.value || 'interp';
            computeTrajectory(20);
            const t = parseFloat(tSlider?.value) || 0;
            setStateAtT(t);
            drawEnergyPlot();
            if (renderer?.updateLabels) renderer.updateLabels(sim);
        });
    }

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
