const { saveObj } = require('../test_utils.js');

// Initialize
const mesh = new MeshBuilder();

function makeQuadSheet(mb) {
    mb.logOperation("makeQuadSheet", { test: "QuadSheet" });

    const n = { x: 10, y: 10 };
    const UVmin = { x: 0.0, y: 0.0 };
    const UVmax = { x: 1.0, y: 1.0 };
    const p00 = new Vec3(0.0, 0.0, 0.0);
    const p01 = new Vec3(80.0, 0.0, 0.0);
    const p10 = new Vec3(0.0, 100.0, 0.0);
    const p11 = new Vec3(120.0, 120.0, 0.0);
    const dirMask = 0b1111; // All directions
    const stickTypes = { x: 0, y: 0, z: 0, w: 0 };

    mb.QuadSheet(n, UVmin, UVmax, p00, p01, p10, p11, dirMask, stickTypes);

    mb.logOperation("makeQuadSheet", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

function makeTubeSheet(mb) {
    mb.logOperation("makeTubeSheet", { test: "TubeSheet" });

    const n = { x: 20, y: 10 };
    const UVmin = { x: 0.0, y: 0.0 };
    const UVmax = { x: 1.0, y: 1.0 };
    const Rs = { x: 10.0, y: 10.0 }; // R1, R2
    const L = 50.0;
    const dirMask = 0b1011;
    const twist = 0.5;
    const stickTypes = { x: 0, y: 0, z: 0, w: 0 };

    mb.TubeSheet(n, UVmin, UVmax, Rs, L, dirMask, twist, stickTypes);

    mb.logOperation("makeTubeSheet", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

function makeTorusSheet(mb) {
    mb.logOperation("makeTorusSheet", { test: "TorusSheet" });

    const n = { x: 20, y: 20 };
    const UVmin = { x: 0.0, y: 0.0 };
    const UVmax = { x: 1.0, y: 1.0 };
    const Rs = { x: 5.0, y: 20.0 }; // r, R
    const dirMask = 0b1011;
    const twist = 0.5;
    const stickTypes = { x: 0, y: 0, z: 0, w: 0 };

    mb.TorusSheet(n, UVmin, UVmax, Rs, dirMask, twist, stickTypes);

    mb.logOperation("makeTorusSheet", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

function makeSlabTube(mb) {
    mb.logOperation("makeSlabTube", { test: "SlabTube" });

    const n = { x: 20, y: 10 };
    const UVmin = { x: 0.0, y: 0.0 };
    const UVmax = { x: 1.0, y: 1.0 };
    const Rs = { x: 10.0, y: 10.0 }; // R1, R2
    const L = 50.0;
    const up = new Vec3(0.0, 0.0, 2.0); // Thickness/Offset
    const dirMask = 0b101010111;
    const stickTypes = { x: 0, y: 0, z: 0, w: 0 };

    mb.SlabTube(n, UVmin, UVmax, Rs, L, up, dirMask, stickTypes);

    mb.logOperation("makeSlabTube", { status: "complete", verts: mb.verts.length, edges: mb.edges.length });
}

try {
    makeQuadSheet(mesh);
    saveObj('output_quadsheet.obj', mesh);
    mesh.clear();

    makeTubeSheet(mesh);
    saveObj('output_tubesheet.obj', mesh);
    mesh.clear();

    makeTorusSheet(mesh);
    saveObj('output_torussheet.obj', mesh);
    mesh.clear();

    makeSlabTube(mesh);
    saveObj('output_slabtube.obj', mesh);
    mesh.clear();

} catch (e) {
    logger.error("Error running test_sheets: " + e.message);
    logger.error(e.stack);
}
