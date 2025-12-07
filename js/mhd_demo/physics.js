// Core physics and state for the MHD / potential-flow demo
// Minimal axisymmetric ring model with springs + gas pressure + iterative current solve

import { Vec3 } from '../common_js/Vec3.js';
import { Logger } from '../common_js/Logger.js';

const logger = (typeof window !== 'undefined' && window.logger) ? window.logger : new Logger();

const MU0 = 4e-7 * Math.PI;

function ellipticK(m) {
    // Carlson symmetric form via AGM approximation (sufficient for demo)
    if (m < 0) return NaN;
    if (m > 1) m = 1;
    const eps = 1e-12;
    let a = 1.0;
    let b = Math.sqrt(1.0 - m);
    let c = Math.sqrt(m);
    let sum = 0.0;
    let twoPow = 1.0;
    for (let i = 0; i < 20; i++) {
        const an = 0.5 * (a + b);
        const bn = Math.sqrt(a * b);
        const cn = 0.5 * (a - b);
        sum += twoPow * cn * cn;
        twoPow *= 2.0;
        a = an; b = bn; c = cn;
        if (cn < eps) break;
    }
    return Math.PI / (2.0 * a);
}

function ellipticE(m) {
    if (m < 0) return NaN;
    if (m > 1) m = 1;
    const eps = 1e-12;
    let a = 1.0;
    let b = Math.sqrt(1.0 - m);
    let c = Math.sqrt(m);
    let sum = 0.0;
    let twoPow = 1.0;
    for (let i = 0; i < 20; i++) {
        const an = 0.5 * (a + b);
        const bn = Math.sqrt(a * b);
        const cn = 0.5 * (a - b);
        sum += twoPow * cn * cn;
        twoPow *= 2.0;
        a = an; b = bn; c = cn;
        if (cn < eps) break;
    }
    return Math.PI * 0.25 * (2.0 * a * a - sum) / a;
}

export function coilField(r, z, I, a) {
    // Returns {Br, Bz} for a circular loop of radius a at origin in r-z plane
    const rp = Math.max(0, r);
    const zp = z;
    const eps = 1e-9;
    const num = 4 * a * rp;
    const den = (a + rp) * (a + rp) + zp * zp + eps;
    let k2 = num / den;
    if (!isFinite(k2)) k2 = 0;
    k2 = Math.min(1 - 1e-9, Math.max(0, k2));
    const K = ellipticK(k2);
    const E = ellipticE(k2);
    const denom = Math.sqrt(den);
    const common = MU0 * I / (2 * Math.PI * (denom + eps));
    const denom2 = ((a - rp) * (a - rp) + zp * zp + eps);
    const invR = 1.0 / (rp + eps);
    const Br = common * (zp * invR) * (-K + (a * a + rp * rp + zp * zp) / denom2 * E);
    const Bz = common * (K + (a * a - rp * rp - zp * zp) / denom2 * E);
    return {
        Br: isFinite(Br) ? Br : 0,
        Bz: isFinite(Bz) ? Bz : 0
    };
}

function revolveVolume(points) {
    // points: array of {r,z} closed loop
    let V = 0;
    const n = points.length;
    for (let i = 0; i < n; i++) {
        const p0 = points[i];
        const p1 = points[(i + 1) % n];
        const dz = p1.z - p0.z;
        const r0 = p0.r;
        const r1 = p1.r;
        V += -Math.PI * dz * (r1 * r1 + r0 * r1 + r0 * r0) / 3;
    }
    return Math.max(V, 1e-6);
}

// Geometry Generators
function makeParabolicCage(f, zMin, zMax, count) {
    const coils = [];
    // User wants "deeper cone", "length bigger than radius".
    // Parabola z = a * r^2.
    // If we want r to be smaller for given z, 'a' should be large.
    // Let's use simple z = k * r^2  => r = sqrt(z/k).
    // Or just parameterized: linear in Z, calculate R.
    // User mentioned: "set of point with increasing diameter (like parabola...)"
    // Let's stick to z = r^2 / (4f) + vertexZ.
    // If f is small, it's narrow. If f is large, it's wide.
    // User said "too open", so make f smaller?
    // r = 2*sqrt(f*z). If f=0.25, r = 2*sqrt(0.25*z) = sqrt(z). At z=4, r=2. Deep enough.
    // Let's try f=0.25.

    for (let i = 0; i < count; i++) {
        const t = i / (count - 1);
        const z = zMin + (zMax - zMin) * t;
        const r = 2 * Math.sqrt(f * z);
        coils.push({ r, z, I: 0 });
    }
    return coils;
}

function makeSphericalPlasma(centerZ, radius, count) {
    const nodes = [];
    const currents = [];
    const restLen = [];
    // User suggestion: r=sin(theta) z=cos(theta) with constant step dTheta.
    // Range: We want "ball". In R-Z, this is a semi-circle from pole to pole.
    // Theta is angle from Z-axis.
    // Range: 0 (North Pole) to PI (South Pole).
    // But r must be >= 0. sin(0)=0, sin(PI)=0. sin(theta) > 0 for 0..PI.
    // So theta goes from 0 to PI.
    // Z = centerZ + R * cos(theta).
    // R_cyl = R * sin(theta). (Cylindrical radius).

    // NOTE: Avoid theta=0 and theta=PI exactly to prevent r=0 singularity/NaNs in Solver.
    const eps = 0.1; // modest padding from axis
    const startTheta = eps;
    const endTheta = Math.PI - eps;

    for (let i = 0; i < count; i++) {
        const t = i / (count - 1);
        const theta = startTheta + (endTheta - startTheta) * t;

        const r = radius * Math.sin(theta);
        const z = centerZ + radius * Math.cos(theta);

        nodes.push({ r, z, vr: 0, vz: 0, m: 1.0 });
        currents.push(0);
        restLen.push(0);
    }

    // Rest lengths (simple chain)
    for (let i = 0; i < count - 1; i++) {
        const n0 = nodes[i];
        const n1 = nodes[i + 1];
        const dist = Math.sqrt(Math.pow(n1.r - n0.r, 2) + Math.pow(n1.z - n0.z, 2));
        restLen[i] = dist;
    }

    return { nodes, currents, restLen };
}

export function initParabolicNozzle(sim) {
    // Reduced counts for testing
    const countCage = 8;
    const countPlasma = 8;

    // 1. SC Coil: At z=-0.5 (near vertex).
    // Vertex of parabola at z=0.
    // User wants "outside ... shielded by the metallic cage".
    // Cage starts at z=0.2 (below).
    // SC at z=-0.5 is safe.
    sim.scCoils = [{ r: 1.2, z: -0.5, I: sim.params.scCurrent }];

    // 2. Cage: Parabola. Deep. f=0.25 (Focus at 0.25).
    // Range z=0.1 to z=3.0.
    // At z=3.0, r = 2*sqrt(0.25*3) = 2*0.866 = 1.73. 
    // Length (3.0) > Radius (1.73). Matches criteria.
    sim.cageCoils = makeParabolicCage(0.25, 0.1, 3.0, countCage);
    sim.cageCurrents = sim.cageCoils.map(() => 0);

    // 3. Plasma: Sphere at focus? 
    // Parabola focus f=0.25.
    // Let's put plasma slightly further out for visibility or exactly at focus.
    // User: "initially we have plasma ball very small in the focus".
    // Dynamic logic usually expands it. Preset starts it somewhere?
    // Let's start it at z=0.8 (inside the cone) with radius 0.4.
    const centerZ = 0.8;
    const radius = 0.4;
    const plasma = makeSphericalPlasma(centerZ, radius, countPlasma);
    sim.plasmaNodes = plasma.nodes;
    sim.plasmaCurrents = plasma.currents;
    sim.restLen = plasma.restLen;

    sim.V0 = revolveVolume(sim.plasmaNodes);
    sim.P0 = sim.params.plasmaPressure0;

    // Set view to static solver mode by default?
    // User requested "Static Solver Logic".
    // We update params/flags if needed.
}

export function initSimulationState() {
    const state = {
        params: {
            scCurrent: 50000.0, // 50 kA
            coilRadius: 1.0,
            plasmaPressure0: 1e4,
            gamma: 5 / 3,
            alphaI: 1e8,
            lambdaI: 1e-3,
            springK: 50,
            damping: 0.995,
            dt: 0.01,
            solverVerb: 1,
        },
        paused: false,
        autoStep: false, // do not step unless user presses Solve (or toggles auto)
        runSolver: false,

        staticSolver: true, // NEW Flag for simplified mode

        scCoils: [{ r: 1.0, z: 0.0, I: 1.0 }],
        cageCoils: [{ r: 1.2, z: -0.4, I: 0 }, { r: 1.2, z: 0.4, I: 0 }],
        plasmaNodes: [],
        plasmaCurrents: [],
        cageCurrents: [],
        restLen: [],
        V0: 0,
        P0: 1e4,
    };

    // Default load (can be overwritten by UI calling initParabolicNozzle)
    initParabolicNozzle(state);

    return state;
}

export function rebuildPlasma(state, radius, n) {
    state.plasmaNodes = [];
    state.plasmaCurrents = [];
    state.restLen = [];
    const N = Math.max(4, Math.floor(n));
    for (let i = 0; i < N; i++) {
        const t = (i / N) * Math.PI * 2;
        const r = radius;
        const z = radius * 0.0 + 0.5 * Math.sin(t);
        state.plasmaNodes.push({ r, z, vr: 0, vz: 0, m: 1.0 });
        state.plasmaCurrents.push(0);
        state.restLen.push(0);
    }
    for (let i = 0; i < N; i++) {
        const n0 = state.plasmaNodes[i];
        const n1 = state.plasmaNodes[(i + 1) % N];
        const dr = n1.r - n0.r;
        const dz = n1.z - n0.z;
        state.restLen[i] = Math.sqrt(dr * dr + dz * dz);
    }
    state.V0 = revolveVolume(state.plasmaNodes);
    state.P0 = state.params.plasmaPressure0;
}

export function applyCoils(state, scList, cageList) {
    state.scCoils = scList.map(c => ({ r: c.r, z: c.z ?? 0, I: c.I ?? 0 }));
    state.cageCoils = cageList.map(c => ({ r: c.r, z: c.z ?? 0, I: c.I ?? 0 }));
    state.cageCurrents = state.cageCoils.map(() => 0);
}

export function sampleFieldAtPoint(r, z, state) {
    // Sample the total magnetic field at point (r, z) in the r-z plane
    // coilField expects z to be the distance from the coil plane (coil at z=0)
    // So we must pass (z - c.z) to get the field relative to each coil's position
    let Br = 0, Bz = 0;
    for (const c of state.scCoils) {
        const f = coilField(r, z - c.z, c.I, c.r);
        Br += f.Br; Bz += f.Bz;
    }
    state.cageCoils.forEach((c, i) => {
        const f = coilField(r, z - c.z, state.cageCurrents[i], c.r);
        Br += f.Br; Bz += f.Bz;
    });
    state.plasmaNodes.forEach((p, i) => {
        const f = coilField(r, z - p.z, state.plasmaCurrents[i], p.r);
        Br += f.Br; Bz += f.Bz;
    });
    return { Br, Bz };
}

function precomputeUnitFields(points, coils) {
    // returns array fields[point][coil] = {Br,Bz}
    // coilField expects z to be distance from coil plane, so we use (p.z - c.z)
    return points.map((p) => coils.map((c) => coilField(p.r, p.z - c.z, 1.0, c.r)));
}

export function stepSimulation(state) {
    if (state.paused) return;
    if (!state.runSolver && !state.autoStep) return;
    state.runSolver = false;
    const { alphaI, lambdaI, springK, damping, dt: dtSim, gamma } = state.params;

    // 1) set SC currents
    state.scCoils.forEach(c => c.I = state.params.scCurrent);

    // Control points: slightly offset to avoid singularity at the wire itself
    // Plasma: just inside (-r). Cage: just outside (+r? or inside?). 
    // Let's bias them away from the source loops.
    const eps = 0.02;
    const controlPlasma = state.plasmaNodes.map(p => ({ r: Math.max(0.01, p.r - eps), z: p.z, targetBr: 0, targetBz: 0 }));
    const controlCage = state.cageCoils.map(c => ({ r: c.r - eps, z: c.z, targetBr: 0, targetBz: 0 }));
    const controls = controlPlasma.concat(controlCage);

    const coilsDynamic = [
        ...state.plasmaNodes.map((p, i) => ({ r: p.r, z: p.z, idx: i, type: 'plasma' })),
        ...state.cageCoils.map((c, i) => ({ r: c.r, z: c.z, idx: i, type: 'cage' })),
    ];

    // Pass coils with z-coordinates so precomputeUnitFields can compute relative distances
    const unitFields = precomputeUnitFields(controls, coilsDynamic);

    // Static Solver Targets:
    // 1. Cage Control Points: Preserve initial flux (from SC only).
    //    Target B = B_sc_only.
    // 2. Plasma Control Points: Expel field (Diamagnetic).
    //    Target B = 0.

    // First, compute B_sc at all control points
    const B_sc_only = controls.map(cp => {
        let Br = 0, Bz = 0;
        state.scCoils.forEach(c => {
            // coilField expects z to be distance from coil plane
            const f = coilField(cp.r, cp.z - c.z, c.I, c.r);
            Br += f.Br; Bz += f.Bz;
        });
        return { Br, Bz };
    });

    // Set targets based on type
    controls.forEach((cp, i) => {
        const isCage = i >= state.plasmaNodes.length; // Cage points appended after plasma
        if (isCage) {
            // Cage: Target = B_sc_only (Flux conservation: Total B stays what it was initially)
            cp.targetBr = B_sc_only[i].Br;
            cp.targetBz = B_sc_only[i].Bz;
        } else {
            // Plasma: Target = 0 (Diamagnetism)
            cp.targetBr = 0;
            cp.targetBz = 0;
        }
    });

    // Precompute constant part of B (fixed fields from SC) needed for residual calc
    // Residual = (A*x + B_fixed) - Target
    // Here B_fixed IS B_sc_only (since cage/plasma are the variables 'x')
    // So Residual = A*x + B_sc_only - Target
    // For Cage: Target = B_sc_only => Residual = A*x
    // For Plasma: Target = 0 => Residual = A*x + B_sc_only

    const fixedB = B_sc_only; // The field from "fixed" sources (SC)

    // Iterative gradient descent for currents
    const nIter = 8;
    for (let iter = 0; iter < nIter; iter++) {
        const errors = controls.map((cp, ci) => {
            const base = fixedB[ci];
            let Br = base.Br, Bz = base.Bz;
            coilsDynamic.forEach((c, j) => {
                const Ij = c.type === 'plasma' ? state.plasmaCurrents[c.idx] : state.cageCurrents[c.idx];
                const uf = unitFields[ci][j];
                Br += Ij * uf.Br;
                Bz += Ij * uf.Bz;
            });
            return { Br: Br - cp.targetBr, Bz: Bz - cp.targetBz };
        });

        // accumulate gradient
        const grad = new Array(coilsDynamic.length).fill(0);
        controls.forEach((cp, ci) => {
            const err = errors[ci];
            coilsDynamic.forEach((c, j) => {
                const uf = unitFields[ci][j];
                grad[j] += err.Br * uf.Br + err.Bz * uf.Bz;
            });
        });

        coilsDynamic.forEach((c, j) => {
            const cur = c.type === 'plasma' ? state.plasmaCurrents[c.idx] : state.cageCurrents[c.idx];
            const delta = -alphaI * grad[j] - lambdaI * cur;
            const next = cur + delta;
            if (c.type === 'plasma') state.plasmaCurrents[c.idx] = isFinite(next) ? next : 0;
            else state.cageCurrents[c.idx] = isFinite(next) ? next : 0;
        });
    }

    // Logging
    if (logger && logger.verb(state.params.solverVerb)) {
        const cpN = controls.length;
        const coilN = coilsDynamic.length;
        // Level 1: Basic
        logger.info(`Solve: coils=${coilN} controlPoints=${cpN}`, Logger.INFO);

        // Level 2: RHS and Solution
        if (state.params.solverVerb >= 2) {
            // Calculate initial errors (RHS) before solve (approximate, using unit fields and zero currents)
            // Actually, the "RHS" in A*x = b is the 'b' vector, which is -B_fixed. 
            // The residuals are A*x - b. 
            // We can log the initial target B vs static B.
            const initErrors = controls.map((cp, ci) => {
                const base = fixedB[ci];
                return Math.sqrt(base.Br * base.Br + base.Bz * base.Bz);
            });
            logger.info(`RHS (Fixed Field Magnitude): ${logger.formatVector(initErrors, 4)}`, Logger.WARN);

            const Ivec = coilsDynamic.map(c => c.type === 'plasma' ? state.plasmaCurrents[c.idx] : state.cageCurrents[c.idx]);
            logger.info(`Solution I: ${logger.formatVector(Ivec, 4)}`, Logger.WARN);
        }

        // Level 3: Matrix
        if (state.params.solverVerb >= 3) {
            const rows = controls.length;
            const cols = coilsDynamic.length;
            const mat = [];
            controls.forEach((cp, ci) => {
                coilsDynamic.forEach((c, j) => {
                    const uf = unitFields[ci][j];
                    // For visualization, we can log magnitude of influence or just one component
                    // Let's log magnitude
                    mat.push(Math.sqrt(uf.Br * uf.Br + uf.Bz * uf.Bz));
                });
            });
            logger.info(`Interaction Matrix |A|:\n${logger.formatMatrix(mat, rows, cols, 3)}`, Logger.INFO);
        }
    }

    if (state.staticSolver) return;

    // Pressure
    const V = revolveVolume(state.plasmaNodes);
    const P = state.P0 * Math.pow(state.V0 / V, gamma);

    // Forces per node
    const N = state.plasmaNodes.length;
    const Fr = new Array(N).fill(0);
    const Fz = new Array(N).fill(0);

    // Springs + pressure normals
    for (let i = 0; i < N; i++) {
        const i1 = (i + 1) % N;
        const p0 = state.plasmaNodes[i];
        const p1 = state.plasmaNodes[i1];
        const dr = p1.r - p0.r;
        const dz = p1.z - p0.z;
        const len = Math.sqrt(dr * dr + dz * dz) + 1e-9;
        const n_r = dz / len; // outward normal (approx)
        const n_z = -dr / len;
        const avgR = 0.5 * (p0.r + p1.r);
        const area = 2 * Math.PI * avgR * len;
        const fp = P * area;
        Fr[i] += fp * n_r * 0.5;
        Fz[i] += fp * n_z * 0.5;
        Fr[i1] += fp * n_r * 0.5;
        Fz[i1] += fp * n_z * 0.5;

        // springs
        const k = springK;
        const fSpring = k * (len - state.restLen[i]);
        const t_r = dr / len;
        const t_z = dz / len;
        Fr[i] += fSpring * t_r;
        Fz[i] += fSpring * t_z;
        Fr[i1] -= fSpring * t_r;
        Fz[i1] -= fSpring * t_z;
    }

    // Lorentz per node
    for (let i = 0; i < N; i++) {
        const p = state.plasmaNodes[i];
        const B = sampleFieldAtPoint(p.r, p.z, state);
        const I = state.plasmaCurrents[i];
        Fr[i] += I * 2 * Math.PI * p.r * B.Bz;
        Fz[i] += I * 2 * Math.PI * p.r * (-B.Br);
    }

    // Integrate
    for (let i = 0; i < N; i++) {
        const p = state.plasmaNodes[i];
        const invm = 1.0 / p.m;
        p.vr = (p.vr + dtSim * Fr[i] * invm) * damping;
        p.vz = (p.vz + dtSim * Fz[i] * invm) * damping;
        p.r = Math.max(0.05, p.r + dtSim * p.vr);
        p.z = p.z + dtSim * p.vz;
        if (!isFinite(p.r)) p.r = 0.1;
        if (!isFinite(p.z)) p.z = 0;
    }
}
