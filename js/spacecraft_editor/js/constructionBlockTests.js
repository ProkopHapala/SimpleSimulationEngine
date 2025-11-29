import { Vec3 } from '../../common_js/Vec3.js';
import { logger } from '../../common_js/Logger.js';

// constructionBlockTests.js

export const ConstructionBlockTests = {

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
            logger.info("Running Girder Simple Test...");

            const p0 = new Vec3(0, 0, 0);
            const p1 = new Vec3(0, 10, 0);
            const p2 = new Vec3(10, 10, 0);

            mesh.girder1(p0, p1, new Vec3(1, 0, 0), 5, 0.5, { x: 1, y: 1, z: 1, w: 1 }, true);
            mesh.girder1(p1, p2, new Vec3(0, 1, 0), 5, 0.5, { x: 1, y: 1, z: 1, w: 1 }, true);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Cube Nodes (Manual)": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Cube Nodes Test...");

            const nodes = [
                new Vec3(0, 0, 0),
                new Vec3(0, 10, 0),
                new Vec3(10, 0, 0)
            ];

            nodes.forEach(p => {
                logger.info(`Adding Cube at ${p.toString()}`);
                mesh.addCube(p, 1.0);
            });

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Bridge Quads Simple": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Bridge Quads Simple Test...");

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
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Bridge Quads Twisted": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Bridge Quads Twisted Test...");

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
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },



        "Bridge Cubes": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Bridge Cubes Test...");

            // Create two cubes
            const p1 = new Vec3(0, 0, 0);
            const p2 = new Vec3(0, 0, 15);
            const size = 2.0;

            // Cube 1
            const chRange1 = mesh.addCube(p1, size, true);
            // Cube 2
            const chRange2 = mesh.addCube(p2, size, true);

            // Bridge them using automatic face detection
            // stickTypes: {x: longitudinal, y: ring, z: spiral, w: internal}
            mesh.bridgeFacingPolygons(p1, p2, chRange1, chRange2, 5, { x: 0, y: 1, z: 2, w: 3 }, { x: 1, y: 1, z: 1, w: 1 });

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        }
    }
};
