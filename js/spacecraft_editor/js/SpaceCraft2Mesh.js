import { Vec3 } from '../../common_js/Vec3.js';
import { Girder, Ring } from './SpaceCraft.js';

/**
 * BuildCraft_blocks_js
 * A modern block-based spacecraft generator for JS.
 * Replaces each Node with a block (cube) and each Girder with a bridged truss.
 */
export function BuildCraft_blocks_js(mesh, craft) {
    mesh.clear();

    const logV = (lvl, msg) => {
        try {
            if (typeof logger !== 'undefined' && logger.uiVerbosity >= lvl) logger.info(msg);
        } catch (e) { /* ignore logging errors */ }
    };

    const defaultNodeSize = 1.0;
    const stickTypes = { x: 1, y: 1, z: 1, w: 1 };
    const bridgeStickTypes = { x: 1, y: 1, z: 1, w: 1 };

    // 1. Generate Node Blocks
    for (const node of craft.nodes) {
        const p = new Vec3(node.pos[0], node.pos[1], node.pos[2]);
        const iv0 = mesh.verts.length / 3;
        const ie0 = mesh.edges.length / 4;
        const ic0 = (mesh.chunks) ? mesh.chunks.length : 0;

        // Using addCube as the block primitive
        const chunkRange = mesh.addCube(p, node.size || defaultNodeSize, true);
        
        node.pointRange = { x: iv0, y: mesh.verts.length / 3 };
        node.stickRange = { x: ie0, y: mesh.edges.length / 4 };
        node.chunkRange = chunkRange;
        logV(1, `[BuildCraft] Node id=${node.id} pos=${node.pos} size=${node.size || defaultNodeSize} verts=${node.pointRange.y - iv0}`);
    }

    // 2. Generate Girder Bridges
    for (const girder of craft.girders) {
        const nA = girder.nodeA;
        const nB = girder.nodeB;
        const pA = new Vec3(nA.pos[0], nA.pos[1], nA.pos[2]);
        const pB = new Vec3(nB.pos[0], nB.pos[1], nB.pos[2]);
        
        const iv0 = mesh.verts.length / 3;
        const ie0 = mesh.edges.length / 4;

        // Bridge facing faces of the two node cubes
        mesh.bridgeFacingPolygons(pA, pB, nA.chunkRange, nB.chunkRange, girder.nseg || 3, bridgeStickTypes, bridgeStickTypes);

        girder.pointRange = { x: iv0, y: mesh.verts.length / 3 };
        girder.stickRange = { x: ie0, y: mesh.edges.length / 4 };
        logV(1, `[BuildCraft] Girder id=${girder.id} nseg=${girder.nseg || 3} nodes=(${nA.id},${nB.id}) verts=${girder.pointRange.y - iv0}`);
    }

    // 3. Generate Ropes
    for (const rope of craft.ropes) {
        const pA = new Vec3(rope.nodeA.pos[0], rope.nodeA.pos[1], rope.nodeA.pos[2]);
        const pB = new Vec3(rope.nodeB.pos[0], rope.nodeB.pos[1], rope.nodeB.pos[2]);
        
        const iv0 = mesh.verts.length / 3;
        const ie0 = mesh.edges.length / 4;

        // Simple rope generator (just a line or rope primitive if available)
        mesh.rope(pA, new Vec3().setSub(pB, pA).normalize(), new Vec3().setSub(pB, pA).norm(), rope.nseg || 10);

        rope.pointRange = { x: iv0, y: mesh.verts.length / 3 };
        rope.stickRange = { x: ie0, y: mesh.edges.length / 4 };
        logV(1, `[BuildCraft] Rope id=${rope.id} nodes=(${rope.nodeA.id},${rope.nodeB.id}) nseg=${rope.nseg || 10} verts=${rope.pointRange.y - iv0}`);
    }

    // 4. Generate Plates (Radiators/Shields) spanning between two bounds (girders/ropes) with span coefs
    if (craft.plates) {
        for (const plate of craft.plates) {
            const { boundA, boundB, spanA, spanB } = plate;
            if (!boundA || !boundB) continue;

            const getEnds = (b) => {
                if (b.nodeA && b.nodeB) {
                    const pA = new Vec3(b.nodeA.pos[0], b.nodeA.pos[1], b.nodeA.pos[2]);
                    const pB = new Vec3(b.nodeB.pos[0], b.nodeB.pos[1], b.nodeB.pos[2]);
                    return { pA, pB };
                }
                return null;
            };

            const endsA = getEnds(boundA);
            const endsB = getEnds(boundB);
            if (!endsA || !endsB) continue;

            const lerp = (p0, p1, c) => new Vec3().setAddMul(p0, new Vec3().setSub(p1, p0), c);
            const a0 = lerp(endsA.pA, endsA.pB, spanA[0]);
            const a1 = lerp(endsA.pA, endsA.pB, spanA[1]);
            const b0 = lerp(endsB.pA, endsB.pB, spanB[0]);
            const b1 = lerp(endsB.pA, endsB.pB, spanB[1]);

            const iv0 = mesh.verts.length / 3;
            const ie0 = mesh.edges.length / 4;

            // quad corners: a0-a1-b1-b0 (consistent winding)
            mesh.ParametricQuadPatch(plate.nx || 2, plate.ny || 2, plate.nz || 1, a0, a1, b1, b0);

            plate.pointRange = { x: iv0, y: mesh.verts.length / 3 };
            plate.stickRange = { x: ie0, y: mesh.edges.length / 4 };
            logV(1, `[BuildCraft] ${plate.kind || 'Plate'} id=${plate.id} boundA=${plate.boundA.id} boundB=${plate.boundB.id} spanA=${spanA} spanB=${spanB} verts=${plate.pointRange.y - iv0}`);
        }
    }

    // 5. Generate Rings
    for (const ring of craft.rings) {
        const iv0 = mesh.verts.length / 3;
        const ie0 = mesh.edges.length / 4;
        
        // Using mesh.wheel from MeshGenerators
        const pos = new Vec3(ring.pos[0], ring.pos[1], ring.pos[2]);
        const p1 = new Vec3(pos.x, pos.y + ring.R, pos.z); // axis hint
        const ax = new Vec3(0, 0, 1);
        mesh.wheel(pos, p1, ax, ring.nseg, ring.wh, stickTypes);

        ring.pointRange = { x: iv0, y: mesh.verts.length / 3 };
        ring.stickRange = { x: ie0, y: mesh.edges.length / 4 };
        logV(1, `[BuildCraft] Ring id=${ring.id} pos=${ring.pos} R=${ring.R} nseg=${ring.nseg} verts=${ring.pointRange.y - iv0}`);
    }

    // 6. Generate Slider Anchors and Visualization
    const sliderPathStickType = 2; // Different material for visualization
    for (const slider of craft.sliders) {
        const boundTo = slider.boundTo;
        if (!boundTo || boundTo.pointRange.x < 0) continue;

        const iv0 = mesh.verts.length / 3;
        const ie0 = mesh.edges.length / 4;

        // Calculate position along path
        // For now, simple linear interpolation between endpoints of girder
        let p;
        if (boundTo instanceof Girder) {
            const pA = new Vec3(boundTo.nodeA.pos[0], boundTo.nodeA.pos[1], boundTo.nodeA.pos[2]);
            const pB = new Vec3(boundTo.nodeB.pos[0], boundTo.nodeB.pos[1], boundTo.nodeB.pos[2]);
            p = new Vec3().setAddMul(pA, new Vec3().setSub(pB, pA), slider.calong);
        } else if (boundTo instanceof Ring) {
            // Placeholder for ring position
            p = new Vec3(boundTo.pos[0], boundTo.pos[1], boundTo.pos[2]);
        } else {
            continue;
        }

        // Generate anchor point (different color/stick type)
        const ivert = mesh.vert(p);
        slider.ivert = ivert;
        
        // Visualization: draw path on the girder/ring in a different color
        if (boundTo instanceof Girder) {
            const pA = new Vec3(boundTo.nodeA.pos[0], boundTo.nodeA.pos[1], boundTo.nodeA.pos[2]);
            const pB = new Vec3(boundTo.nodeB.pos[0], boundTo.nodeB.pos[1], boundTo.nodeB.pos[2]);
            // Draw a line segment parallel to the girder
            const offset = 0.2;
            const pA_off = new Vec3(pA.x, pA.y + offset, pA.z);
            const pB_off = new Vec3(pB.x, pB.y + offset, pB.z);
            const ivA = mesh.vert(pA_off);
            const ivB = mesh.vert(pB_off);
            mesh.edge(ivA, ivB, sliderPathStickType);
        }

        slider.pointRange = { x: iv0, y: mesh.verts.length / 3 };
        slider.stickRange = { x: ie0, y: mesh.edges.length / 4 };
        logV(1, `[BuildCraft] Slider id=${slider.id} boundTo=${slider.boundTo.id} calong=${slider.calong} verts=${slider.pointRange.y - iv0}`);
    }
}
