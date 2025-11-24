const { saveObj } = require('../test_utils.js');

// Initialize
const mesh = new MeshBuilder();

function makePanel(mb) {
    mb.logOperation("makePanel", { test: "QuadPanel" });

    // C++ call: QuadPanel(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {80.0,0.0,0.0}, {0.0,100.0,0.0}, {120.0,120.0,0.0}, 10.0, Quat4i{0,0,0,0} );

    const n = { x: 10, y: 10 };
    const UVmin = { x: 0.0, y: 0.0 };
    const UVmax = { x: 1.0, y: 1.0 };
    const p00 = new Vec3(0.0, 0.0, 0.0);
    const p01 = new Vec3(80.0, 0.0, 0.0);
    const p10 = new Vec3(0.0, 100.0, 0.0);
    const p11 = new Vec3(120.0, 120.0, 0.0);
    const width = 10.0;
    const stickTypes = { x: 0, y: 0, z: 0, w: 0 };

    mb.QuadPanel(n, UVmin, UVmax, p00, p01, p10, p11, width, stickTypes);

    mb.logOperation("makePanel", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

function makeSlab(mb) {
    mb.logOperation("makeSlab", { test: "QuadSlab" });

    // C++ call: QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, {0.333333,0.333333,7.0}, 0b101010111, Quat4i{0,0,0,0} );

    const n = { x: 10, y: 10 };
    const UVmin = { x: 0.0, y: 0.0 };
    const UVmax = { x: 1.0, y: 1.0 };
    const p00 = new Vec3(0.0, 0.0, 0.0);
    const p01 = new Vec3(86.6, 50.0, 0.0);
    const p10 = new Vec3(0.0, 100.0, 0.0);
    const p11 = new Vec3(86.6, 150.0, 0.0);
    const up = new Vec3(0.33, 0.33, 7.0);
    const dirMask = 0b101010111;
    const stickTypes = { x: 0, y: 0, z: 0, w: 0 };

    mb.QuadSlab(n, UVmin, UVmax, p00, p01, p10, p11, up, dirMask, stickTypes);

    mb.logOperation("makeSlab", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

try {
    makePanel(mesh);
    saveObj('output_panel.obj', mesh);

    mesh.clear();
    makeSlab(mesh);
    saveObj('output_slab.obj', mesh);
} catch (e) {
    logger.error("Error running test_panel: " + e.message);
    logger.error(e.stack);
}
