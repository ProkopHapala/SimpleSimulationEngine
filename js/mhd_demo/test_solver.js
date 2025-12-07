
import { initSimulationState, stepSimulation, applyCoils } from './physics.js';
import { Logger } from '../common_js/Logger.js';

// Setup global logger for Node environment
global.logger = new Logger();
global.logger.setConsoleVerbosity(Logger.INFO);

console.log("=== MHD Solver Offline Test ===");

// 1. Initialize State
const sim = initSimulationState();
console.log("State initialized.");

// 2. Configure Coils (Use default Parabolic Nozzle preset)
// const scCoils = [{ r: 1.0, z: -0.5, I: 1.0 }, { r: 1.0, z: 0.5, I: 1.0 }];
// const cageCoils = [{ r: 1.2, z: 0.0, I: 0 }];
// applyCoils(sim, scCoils, cageCoils);
console.log(`Coils configured: SC=${sim.scCoils.length}, Cage=${sim.cageCoils.length}, Plasma=${sim.plasmaNodes.length}`);

// 3. Force Solver to Run
sim.runSolver = true;
sim.params.solverVerb = 2; // Verify logging

// 4. Step Simulation
console.log("Running stepSimulation...");
stepSimulation(sim);

// 5. Check Results
console.log("\n--- Results ---");
console.log("Cage Currents:", sim.cageCurrents);
console.log("Plasma Currents (First 5):", sim.plasmaCurrents.slice(0, 5));

const p0 = sim.plasmaNodes[0];
console.log(`Plasma Node 0 Pos: r=${p0.r.toFixed(4)}, z=${p0.z.toFixed(4)}`);

if (sim.cageCurrents.some(i => !isFinite(i)) || sim.plasmaCurrents.some(i => !isFinite(i))) {
    console.error("FAIL: NaN detected in currents!");
    process.exit(1);
}

console.log("SUCCESS: Simulation stepped without errors.");
