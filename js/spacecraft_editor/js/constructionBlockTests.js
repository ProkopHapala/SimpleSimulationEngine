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
        },

        "ParametricParabola Dish": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running ParametricParabola Dish Test (parabolic pusher-plate)...");

            // Tuned parameters taken from MeshGenTestGUI ParametricParabola
            // params: {"ny":4,"nx":7,"nx2":25,"L":3,"R1":0.5,"R2":3,"dirMask":"0111"}
            const nRows = 4;
            const nInner = 7;   // inner ring segments
            const nOuter = 25;  // outer ring segments
            const R1 = 0.5;     // inner radius (hole)
            const R2 = 3.0;     // outer radius
            const L  = 3.0;     // axial height at R2

            // Axial shift along +Z to place the dish above a tube; tweak as needed
            const zShift = 2.0;

            const iv0 = mesh.verts.length;
            mesh.ParametricParabolaPatch(nOuter, nInner, nRows, R1, R2, L);
            const iv1 = mesh.verts.length;

            // Apply axial shift to all newly created vertices
            for (let i = iv0; i < iv1; i++) {
                mesh.verts[i].pos.z += zShift;
            }

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        },

        "TubeSheetBond Hex Ring": (engine) => {
            const mesh = ConstructionBlockTests.setup(engine);
            logger.info("Running TubeSheetBond Hex Ring Test (recoil damper outer tube)...");

            // Parameters mirrored from MeshGenTestGUI TubeSheetBond defaults
            const dirMaskStr = "1111";
            const dirMask = parseInt(dirMaskStr, 2);
            const twist = 0.0;

            const clearance = 0.1; // radial gap between outer tube inner radius and inner triangular tube

            const n1 = { x: 2, y: 6 };
            const n2 = { x: 3, y: 6 };
            const UVmin = { x: 0, y: 0 };
            const UVmax = { x: 1, y: 1 };

            const R1 = 1.732;
            const L1 = 1.0;
            const R2 = 1.0;
            const L2 = 2.0;
            const Rs1 = { x: R1, y: R1 };
            const Rs2 = { x: R2, y: R2 };

            const offX = -0.5;
            const offY = 0.5;
            const du = offX / (n2.x - 1.0);
            const dv = offY / (n2.y);
            const UVmin2 = { x: du,     y: dv };
            const UVmax2 = { x: 1 + du, y: 1 + dv };

            const Rcut = 2.0;
            const dR = 0.01;
            const iFamily = 0;
            const stickTypes = { x: 1, y: 0, z: 0, w: 0 };

            // First tube
            const iv0_1 = mesh.verts.length;
            mesh.TubeSheet(n1, UVmin, UVmax, Rs1, L1, dirMask, twist, stickTypes);
            const nVerts1 = n1.x * n1.y;
            const sel1 = new Array(nVerts1);
            for (let i = 0; i < nVerts1; i++) sel1[i] = iv0_1 + i;

            // Second tube with shifted UV window


            const iv0_2 = mesh.verts.length;
            mesh.TubeSheet(n2, UVmin2, UVmax2, Rs2, L2, dirMask, twist, stickTypes);
            const nVerts2 = n2.x * n2.y;
            const sel2 = new Array(nVerts2);
            for (let i = 0; i < nVerts2; i++) sel2[i] = iv0_2 + i;

            // Bond-length analysis and connections
            const groups = [];
            const uniqLs = mesh.mapBondLengthsFromVertexSelections(sel1, sel2, Rcut, dR, groups);
            if (uniqLs && uniqLs.length > 0 && groups && groups.length > 0) {
                const k = Math.max(0, Math.min(iFamily, uniqLs.length - 1));
                const fam = groups[k];
                if (fam && fam.length > 0) {
                    mesh.addEdgesFromPairs(fam, stickTypes.x);
                }
            }

            // Inner triangular tube (single TubeSheet with 3-sided cross-section) placed inside the outer hex ring
            const nTri = { x: 10, y: 3 };
            const Rtri = Math.max(0.0, R2 - clearance);
            const RsTri = { x: Rtri, y: Rtri };
            // Rotate inner triangle by half a hex segment (30deg) via UV shift: dv = 1/(2*3) = 1/6
            const dvTri = 1.0 / (2.0 * nTri.y)*0.5;
            const axshift = -0.9;
            const UVminTri = { x: axshift, y: dvTri };
            const UVmaxTri = { x: 1 + axshift, y: 1 + dvTri };
            const dirMaskTri = parseInt("1111", 2);
            const iv0_tri = mesh.verts.length;
            mesh.TubeSheet(nTri, UVminTri, UVmaxTri, RsTri, 20.0, dirMaskTri, 0.0, stickTypes);

            // --- Parabolic pusher-plate / plasma nozzle dish above the damper ---
            // Tuned parameters taken from MeshGenTestGUI ParametricParabola
            // params: {"ny":4,"nx":7,"nx2":25,"L":3,"R1":0.5,"R2":3,"dirMask":"0111"}
            const nRows  = 4;
            const nInner = 7;   // inner ring segments
            const nOuter = 25;  // outer ring segments
            const R1dish = 1.732; // inner radius (hole)
            const R2dish = 1.732*6; // outer radius
            const Ldish  = 3.0*6; // axial height at R2

            // Place dish just above the end of the tubes along +Z
            const zShiftDish = L2*0.5;
            const iv0_dish = mesh.verts.length;
            mesh.ParametricParabolaPatch(nOuter, nInner, nRows, R1dish, R2dish, Ldish);
            const iv1_dish = mesh.verts.length;
            for (let i = iv0_dish; i < iv1_dish; i++) {
                mesh.verts[i].pos.z += zShiftDish;
            }

            if (window.renderer) window.renderer.updateGeometry(mesh);
            logger.info(`Test Complete. Verts: ${mesh.verts.length}, Edges: ${mesh.edges.length}`);
        }
    }
};
