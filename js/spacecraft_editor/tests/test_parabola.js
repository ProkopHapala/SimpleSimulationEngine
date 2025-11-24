const { saveObj } = require('../test_utils.js');

// Initialize
const mesh = new MeshBuilder();

// Implement Parabola Logic
function makeParabola(mb) {
    mb.logOperation("makeParabola", { test: "Parabola_Wire_new" });

    // C++ call example: 
    // Parabola_Wire_new( truss, {3,6}, Vec2f{0.2,0.0}, Vec2f{1.0,M_PI*2.0}, 10.0, 10.0, 0.0, WireFlags::STAR  );

    const n = { x: 3, y: 6 };
    const UVmin = { x: 0.2, y: 0.0 };
    const UVmax = { x: 1.0, y: Math.PI * 2.0 };
    const R = 10.0;
    const L = 10.0;
    const voff = 0.0;
    const flags = WireFlags.STAR;

    mb.Parabola_Wire_new(n, UVmin, UVmax, R, L, voff, flags);

    mb.logOperation("makeParabola", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

// Run
try {
    makeParabola(mesh);
    // Export
    saveObj('output_parabola.obj', mesh);
} catch (e) {
    logger.error("Error running makeParabola: " + e.message);
    logger.error(e.stack);
}
