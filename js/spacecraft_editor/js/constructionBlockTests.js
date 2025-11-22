
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

            const p0 = new Vec3(0, 0, 0);
            const p1 = new Vec3(0, 10, 0);
            const p2 = new Vec3(10, 10, 0);

            mesh.girder1(p0, p1, new Vec3(1, 0, 0), 5, 0.5);
            mesh.girder1(p1, p2, new Vec3(0, 1, 0), 5, 0.5);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Cube Nodes (Manual)": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            window.logger.info("Running Cube Nodes Test...");

            const nodes = [
                new Vec3(0, 0, 0),
                new Vec3(0, 10, 0),
                new Vec3(10, 0, 0)
            ];

            nodes.forEach(p => {
                window.logger.info(`Adding Cube at ${p.toString()}`);
                mesh.addCube(p, 1.0);
            });

            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Bridge Quads Simple": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            window.logger.info("Running Bridge Quads Simple Test...");

            // Create two quads facing each other
            const s = 2.0; // half-size

            // Quad 1 at z=0
            const v0 = mesh.vert(new Vec3(-s, -s, 0));
            const v1 = mesh.vert(new Vec3(s, -s, 0));
            const v2 = mesh.vert(new Vec3(s, s, 0));
            const v3 = mesh.vert(new Vec3(-s, s, 0));

            // Quad 2 at z=10
            const v4 = mesh.vert(new Vec3(-s, -s, 10));
            const v5 = mesh.vert(new Vec3(s, -s, 10));
            const v6 = mesh.vert(new Vec3(s, s, 10));
            const v7 = mesh.vert(new Vec3(-s, s, 10));

            // Bridge them
            const q1 = { x: v0, y: v1, z: v2, w: v3 };
            const q2 = { x: v4, y: v5, z: v6, w: v7 };

            mesh.bridge_quads(q1, q2, 4, { x: 0, y: 1, z: 2, w: 3 }, { x: 1, y: 1, z: 1, w: 1 }, true);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Bridge Quads Twisted": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            window.logger.info("Running Bridge Quads Twisted Test...");

            const s = 2.0;

            // Quad 1 at z=0 (normal orientation)
            const v0 = mesh.vert(new Vec3(-s, -s, 0));
            const v1 = mesh.vert(new Vec3(s, -s, 0));
            const v2 = mesh.vert(new Vec3(s, s, 0));
            const v3 = mesh.vert(new Vec3(-s, s, 0));

            // Quad 2 at z=10 (rotated 45 degrees)
            const angle = Math.PI / 4;
            const c = Math.cos(angle);
            const sn = Math.sin(angle);
            const v4 = mesh.vert(new Vec3(-s * c - (-s) * sn, -s * sn + (-s) * c, 10));
            const v5 = mesh.vert(new Vec3(s * c - (-s) * sn, s * sn + (-s) * c, 10));
            const v6 = mesh.vert(new Vec3(s * c - s * sn, s * sn + s * c, 10));
            const v7 = mesh.vert(new Vec3(-s * c - s * sn, -s * sn + s * c, 10));

            const q1 = { x: v0, y: v1, z: v2, w: v3 };
            const q2 = { x: v4, y: v5, z: v6, w: v7 };

            // Test alignment
            mesh.bridge_quads(q1, q2, 6, { x: 0, y: 1, z: 2, w: 3 }, { x: 1, y: 1, z: 0, w: 0 }, true);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Cube Nodes with Bridges": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            window.logger.info("Running Cube Nodes with Bridges Test...");

            // Create two cubes
            const p1 = new Vec3(0, 0, 0);
            const p2 = new Vec3(0, 0, 15);
            const size = 2.0;
            const s = size * 0.5;

            // Cube 1
            const c1v0 = mesh.vert(new Vec3(p1.x - s, p1.y - s, p1.z - s));
            const c1v1 = mesh.vert(new Vec3(p1.x + s, p1.y - s, p1.z - s));
            const c1v2 = mesh.vert(new Vec3(p1.x + s, p1.y + s, p1.z - s));
            const c1v3 = mesh.vert(new Vec3(p1.x - s, p1.y + s, p1.z - s));
            const c1v4 = mesh.vert(new Vec3(p1.x - s, p1.y - s, p1.z + s));
            const c1v5 = mesh.vert(new Vec3(p1.x + s, p1.y - s, p1.z + s));
            const c1v6 = mesh.vert(new Vec3(p1.x + s, p1.y + s, p1.z + s));
            const c1v7 = mesh.vert(new Vec3(p1.x - s, p1.y + s, p1.z + s));

            // Cube 1 edges
            mesh.edge(c1v0, c1v1); mesh.edge(c1v1, c1v2); mesh.edge(c1v2, c1v3); mesh.edge(c1v3, c1v0);
            mesh.edge(c1v4, c1v5); mesh.edge(c1v5, c1v6); mesh.edge(c1v6, c1v7); mesh.edge(c1v7, c1v4);
            mesh.edge(c1v0, c1v4); mesh.edge(c1v1, c1v5); mesh.edge(c1v2, c1v6); mesh.edge(c1v3, c1v7);

            // Cube 2
            const c2v0 = mesh.vert(new Vec3(p2.x - s, p2.y - s, p2.z - s));
            const c2v1 = mesh.vert(new Vec3(p2.x + s, p2.y - s, p2.z - s));
            const c2v2 = mesh.vert(new Vec3(p2.x + s, p2.y + s, p2.z - s));
            const c2v3 = mesh.vert(new Vec3(p2.x - s, p2.y + s, p2.z - s));
            const c2v4 = mesh.vert(new Vec3(p2.x - s, p2.y - s, p2.z + s));
            const c2v5 = mesh.vert(new Vec3(p2.x + s, p2.y - s, p2.z + s));
            const c2v6 = mesh.vert(new Vec3(p2.x + s, p2.y + s, p2.z + s));
            const c2v7 = mesh.vert(new Vec3(p2.x - s, p2.y + s, p2.z + s));

            // Cube 2 edges
            mesh.edge(c2v0, c2v1); mesh.edge(c2v1, c2v2); mesh.edge(c2v2, c2v3); mesh.edge(c2v3, c2v0);
            mesh.edge(c2v4, c2v5); mesh.edge(c2v5, c2v6); mesh.edge(c2v6, c2v7); mesh.edge(c2v7, c2v4);
            mesh.edge(c2v0, c2v4); mesh.edge(c2v1, c2v5); mesh.edge(c2v2, c2v6); mesh.edge(c2v3, c2v7);

            // Bridge top face of cube 1 to bottom face of cube 2
            const q1 = { x: c1v4, y: c1v5, z: c1v6, w: c1v7 }; // top face of cube 1
            const q2 = { x: c2v0, y: c2v1, z: c2v2, w: c2v3 }; // bottom face of cube 2

            mesh.bridge_quads(q1, q2, 5, { x: 0, y: 1, z: 2, w: 3 }, { x: 1, y: 1, z: 1, w: 1 }, true);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            window.logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        }
    }
};
