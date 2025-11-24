const { saveObj } = require('../test_utils.js');

// Initialize
const mesh = new MeshBuilder();

function test_bevel_cube(mb) {
    mb.logOperation("test_bevel_cube", { test: "Bevel vertices of a cube" });

    // Create a cube
    mb.addCube(new Vec3(0, 0, 0), 4, false);

    // Build edge adjacency
    mb.build_edgesOfVerts();

    // Bevel all 8 vertices
    const L = 0.3;  // Bevel distance
    const h = 0.0;  // No offset along normal
    const bPoly = true;  // Create polygon faces
    const bEdgeWedge = false;  // Don't connect to original vertex

    logger.info(`Beveling ${mb.verts.length} vertices...`);

    // Bevel each vertex with its normal pointing outward from center
    for (let iv = 0; iv < 8; iv++) {
        const pos = mb.verts[iv].pos;
        const nor = pos.clone();
        nor.normalize();  // Normal points away from origin

        mb.bevel_vert(iv, L, h, bPoly, bEdgeWedge, nor);
    }

    mb.logOperation("test_bevel_cube", {
        status: "complete",
        verts: mb.verts.length,
        edges: mb.edges.length,
        chunks: mb.chunks.length
    });
}

function test_bevel_octahedron(mb) {
    mb.logOperation("test_bevel_octahedron", { test: "Bevel octahedron" });

    // Create an octahedron (vertices at unit axes)
    const s = 2.0;
    mb.vert(new Vec3(s, 0, 0));   // 0
    mb.vert(new Vec3(-s, 0, 0));  // 1
    mb.vert(new Vec3(0, s, 0));   // 2
    mb.vert(new Vec3(0, -s, 0));  // 3
    mb.vert(new Vec3(0, 0, s));   // 4
    mb.vert(new Vec3(0, 0, -s));  // 5

    // Edges (12 edges forming 8 triangular faces)
    // Top pyramid (apex at +Z)
    mb.edge(0, 2); mb.edge(2, 1); mb.edge(1, 3); mb.edge(3, 0);
    mb.edge(4, 0); mb.edge(4, 2); mb.edge(4, 1); mb.edge(4, 3);
    // Bottom pyramid (apex at -Z)
    mb.edge(5, 0); mb.edge(5, 2); mb.edge(5, 1); mb.edge(5, 3);

    // Build edge adjacency
    mb.build_edgesOfVerts();

    // Bevel the 6 vertices
    const L = 0.4;
    const h = 0.0;
    const bPoly = true;
    const bEdgeWedge = false;

    for (let iv = 0; iv < 6; iv++) {
        const pos = mb.verts[iv].pos;
        const nor = pos.clone();
        nor.normalize();

        mb.bevel_vert(iv, L, h, bPoly, bEdgeWedge, nor);
    }

    mb.logOperation("test_bevel_octahedron", {
        status: "complete",
        verts: mb.verts.length,
        edges: mb.edges.length,
        chunks: mb.chunks.length
    });
}

try {
    test_bevel_cube(mesh);
    saveObj('output_bevel_cube.obj', mesh);
    mesh.clear();

    test_bevel_octahedron(mesh);
    saveObj('output_bevel_octahedron.obj', mesh);
    mesh.clear();

} catch (e) {
    logger.error("Error running test_bevel: " + e.message);
    logger.error(e.stack);
}
