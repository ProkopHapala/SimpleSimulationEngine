
import {
    SimVec3Obj_Manual, SimVec3Obj_Methods,
    SimVec3Arr_Manual, SimVec3Arr_Methods,
    SimVec3Lst_Manual, SimVec3Lst_Methods,
    SimTHREE_Manual, SimTHREE_Methods,
    SimGLM_Manual, SimGLM_Methods
} from './NBody_bench.js';

const N = 1000;
const STEPS = 100;

// ============================================================================
// Runner
// ============================================================================

function runBenchmark(SimClass, name) {
    const sim = new SimClass(N);
    sim.init();

    // Warmup
    sim.run(1, 0.01);

    const start = performance.now();
    sim.run(STEPS, 0.01);
    const duration = performance.now() - start;

    const stepsPerSec = (STEPS / duration) * 1000;
    const interactionsPerSec = (N * N * STEPS) / (duration / 1000);
    console.log(`${name.padEnd(25)} | Time: ${duration.toFixed(2).padStart(8)} ms | Steps/s: ${stepsPerSec.toFixed(2).padStart(6)} | Int/s: ${interactionsPerSec.toExponential(2)}`);
}

console.log("Running Benchmarks (N=" + N + ", Steps=" + STEPS + ")...");
console.log("--------------------------------------------------------------------------------");

runBenchmark(SimVec3Obj_Manual, "Vec3Obj (Manual)");
runBenchmark(SimVec3Obj_Methods, "Vec3Obj (Methods)");
console.log("-");
runBenchmark(SimVec3Lst_Manual, "Vec3Lst (Manual)");
runBenchmark(SimVec3Lst_Methods, "Vec3Lst (Methods)");
console.log("-");
runBenchmark(SimVec3Arr_Manual, "Vec3Arr (Manual)");
runBenchmark(SimVec3Arr_Methods, "Vec3Arr (Methods)");
console.log("-");
runBenchmark(SimTHREE_Manual, "THREE (Manual)");
runBenchmark(SimTHREE_Methods, "THREE (Methods)");
console.log("-");
runBenchmark(SimGLM_Manual, "GLM (Manual)");
runBenchmark(SimGLM_Methods, "GLM (Methods)");
