
// constructionBlockTests.js

// constructionBlockTests.js

window.ConstructionBlockTests = {

    // Helper to clear and setup
    setup: (engine) => {
        if (document.getElementById('chkAutoClear').checked) {
            engine.reset();
        }
        return engine.mesh;
    },

    // Tests
    tests: {
        "Girder Simple": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            window.logger.info("Running Girder Simple Test...");

            const p0 = [0, 0, 0];
            const p1 = [0, 10, 0];
            const p2 = [10, 10, 0];

            // Use MeshBuilder's girder1
            mesh.girder1(p0, p1, [1, 0, 0], 5, 0.5);
            mesh.girder1(p1, p2, [0, 1, 0], 5, 0.5);

            // Update Renderer directly (bypass Engine.processCommands which rebuilds from Craft)
            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length / 3}, Edges: ${mesh.edges.length / 4}`);
        },

        "Cube Nodes (Manual)": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            window.logger.info("Running Cube Nodes Test...");

            // Manually build what -cube_nodes does roughly
            // 1. Create nodes (cubes)
            const nodes = [
                [0, 0, 0],
                [0, 10, 0],
                [10, 0, 0]
            ];

            // 2. For each node, add a cube
            nodes.forEach(p => {
                window.logger.info(`Adding Cube at [${p}]`);
                mesh.addCube(p, 1.0);
            });

            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length / 3}, Edges: ${mesh.edges.length / 4}`);
        }
    }
};
