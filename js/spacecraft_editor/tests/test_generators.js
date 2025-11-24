const { saveObj } = require('../test_utils.js');

// Initialize
const mesh = new MeshBuilder();

function test_skeleton(mb) {
    mb.logOperation("test_skeleton", { test: "Skeleton girder network" });

    // Define nodes (positions and sizes)
    const nodes = [
        { pos: new Vec3(0, 0, 0), size: 1.0 },
        { pos: new Vec3(10, 0, 0), size: 1.0 },
        { pos: new Vec3(10, 10, 0), size: 1.0 },
        { pos: new Vec3(0, 10, 0), size: 1.0 },
        { pos: new Vec3(5, 5, 10), size: 1.0 }
    ];

    // Define connections (girders between nodes)
    const connections = [
        { x: 0, y: 1 },  // Bottom square
        { x: 1, y: 2 },
        { x: 2, y: 3 },
        { x: 3, y: 0 },
        { x: 0, y: 4 },  // Pyramid from bottom to apex
        { x: 1, y: 4 },
        { x: 2, y: 4 },
        { x: 3, y: 4 }
    ];

    mb.skeleton(nodes, connections, {
        nseg: 3,
        width: 0.3,
        stickTypes: { x: 0, y: 0, z: 0, w: 0 }
    });

    mb.logOperation("test_skeleton", {
        status: "complete",
        verts: mb.verts.length,
        edges: mb.edges.length
    });
}

function test_block(mb) {
    mb.logOperation("test_block", { test: "Subdivided block" });

    mb.block(
        new Vec3(0, 0, 0),
        new Vec3(10, 8, 6),
        { x: 3, y: 3, z: 2 },
        true  // wireframe
    );

    mb.logOperation("test_block", {
        status: "complete",
        verts: mb.verts.length,
        edges: mb.edges.length
    });
}

try {
    test_skeleton(mesh);
    saveObj('output_skeleton.obj', mesh);
    mesh.clear();

    test_block(mesh);
    saveObj('output_block.obj', mesh);
    mesh.clear();

} catch (e) {
    logger.error("Error running test_generators: " + e.message);
    logger.error(e.stack);
}
