import { Vec3 } from '../../common_js/Vec3.js';
import { logger } from '../../common_js/Logger.js';

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

        "Ropes V-Shape": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Ropes V-Shape Test (shared endpoint)...");

            const nseg = 10;

            // Shared root
            const A = new Vec3(0, 0, 0);
            const B = new Vec3(0, 10, 0);
            const C = new Vec3(10, 10, 0);

            // First polyline A->B
            const vertsAB = [];
            for (let i = 0; i <= nseg; i++) {
                const t = i / nseg;
                const p = new Vec3(
                    A.x * (1 - t) + B.x * t,
                    A.y * (1 - t) + B.y * t,
                    A.z * (1 - t) + B.z * t
                );
                vertsAB.push(mesh.vert(p));
            }
            for (let i = 0; i < nseg; i++) {
                mesh.edge(vertsAB[i], vertsAB[i + 1]);
            }

            // Second polyline A->C (shares vertex 0 with first polyline)
            const vertsAC = [vertsAB[0]];
            for (let i = 1; i <= nseg; i++) {
                const t = i / nseg;
                const p = new Vec3(
                    A.x * (1 - t) + C.x * t,
                    A.y * (1 - t) + C.y * t,
                    A.z * (1 - t) + C.z * t
                );
                vertsAC.push(mesh.vert(p));
            }
            for (let i = 0; i < nseg; i++) {
                mesh.edge(vertsAC[i], vertsAC[i + 1]);
            }

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Ropes Parallel": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Ropes Parallel Test (||)...");

            const nseg = 10;

            const A1 = new Vec3(-5, 0, 0);
            const B1 = new Vec3(-5, 10, 0);
            const A2 = new Vec3(5, 0, 0);
            const B2 = new Vec3(5, 10, 0);

            const buildPolyline = (P0, P1) => {
                const verts = [];
                for (let i = 0; i <= nseg; i++) {
                    const t = i / nseg;
                    const p = new Vec3(
                        P0.x * (1 - t) + P1.x * t,
                        P0.y * (1 - t) + P1.y * t,
                        P0.z * (1 - t) + P1.z * t
                    );
                    verts.push(mesh.vert(p));
                }
                for (let i = 0; i < nseg; i++) {
                    mesh.edge(verts[i], verts[i + 1]);
                }
                return verts;
            };

            buildPolyline(A1, B1);
            buildPolyline(A2, B2);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Ropes V-Shape + Plate": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Ropes V-Shape + Plate Test (mismatched segments)...");

            const nsegAB = 12; // more segments on AB
            const nsegAC = 6;  // fewer segments on AC

            // Shared root
            const A = new Vec3(0, 0, 0);
            const B = new Vec3(0, 10, 0);
            const C = new Vec3(10, 10, 0);

            // First polyline A->B (fine sampling)
            const vertsAB = [];
            for (let i = 0; i <= nsegAB; i++) {
                const t = i / nsegAB;
                const p = new Vec3(
                    A.x * (1 - t) + B.x * t,
                    A.y * (1 - t) + B.y * t,
                    A.z * (1 - t) + B.z * t
                );
                vertsAB.push(mesh.vert(p));
            }
            for (let i = 0; i < nsegAB; i++) {
                mesh.edge(vertsAB[i], vertsAB[i + 1]);
            }

            // Second polyline A->C (coarser sampling, shares vertex 0 with first polyline)
            const vertsAC = [vertsAB[0]];
            for (let i = 1; i <= nsegAC; i++) {
                const t = i / nsegAC;
                const p = new Vec3(
                    A.x * (1 - t) + C.x * t,
                    A.y * (1 - t) + C.y * t,
                    A.z * (1 - t) + C.z * t
                );
                vertsAC.push(mesh.vert(p));
            }
            for (let i = 0; i < nsegAC; i++) {
                mesh.edge(vertsAC[i], vertsAC[i + 1]);
            }

            // Corners: [tipA, shared, tipB]
            const tipA = vertsAB[nsegAB];
            const shared = vertsAB[0];
            const tipB = vertsAC[nsegAC];

            // Radius and sampling parameters for strip extraction
            const r = 0.3;
            const maxPerStrip = Math.max(nsegAB, nsegAC) + 1;

            mesh.triPlateBetweenEdges([tipA, shared, tipB], r, maxPerStrip);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "Ropes Parallel + Plate": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Ropes Parallel + Plate Test (ParametricQuadPatch)...");

            const A1 = new Vec3(-5, 0, 0);   // bottom-left
            const B1 = new Vec3(-5, 10, 0);  // top-left
            const A2 = new Vec3(5, 0, 0);    // bottom-right
            const B2 = new Vec3(5, 10, 0);   // top-right

            const nTop = 10;    // along bottom edge A1-A2
            const nBottom = 6;  // along top edge B1-B2
            const nRows = 8;    // between bottom and top

            mesh.ParametricQuadPatch(nTop, nBottom, nRows, A1, B1, A2, B2);

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
        },

        "Triangulated Variations": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running Triangulated Variations Test (ParametricQuadPatch on square)...");

            const s = 5.0;
            const p00 = new Vec3(-s, 0, -s);
            const p01 = new Vec3(-s, 0, s);
            const p10 = new Vec3(s, 0, -s);
            const p11 = new Vec3(s, 0, s);

            const nTop = 8;
            const nBottom = 20;
            const nRows = 8;

            mesh.ParametricQuadPatch(nTop, nBottom, nRows, p00, p01, p10, p11);

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        }
    }
};
