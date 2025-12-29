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
        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;
        const ic0 = (mesh.chunks) ? mesh.chunks.length : 0;

        // Using addCube as the block primitive
        const chunkRange = mesh.addCube(p, node.size || defaultNodeSize, true);
        
        node.pointRange = { x: iv0, y: mesh.verts.length };
        node.stickRange = { x: ie0, y: mesh.edges.length };
        node.chunkRange = chunkRange;
        logV(1, `[BuildCraft] Node id=${node.id} pos=${node.pos} size=${node.size || defaultNodeSize} verts=${node.pointRange.y - iv0}`);
    }

    // 2. Generate Girder Bridges
    for (const girder of craft.girders) {
        const nA = girder.nodeA;
        const nB = girder.nodeB;
        const pA = new Vec3(nA.pos[0], nA.pos[1], nA.pos[2]);
        const pB = new Vec3(nB.pos[0], nB.pos[1], nB.pos[2]);
        
        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;

        // Bridge facing faces of the two node cubes
        mesh.bridgeFacingPolygons(pA, pB, nA.chunkRange, nB.chunkRange, girder.nseg || 3, bridgeStickTypes, bridgeStickTypes);

        girder.pointRange = { x: iv0, y: mesh.verts.length };
        girder.stickRange = { x: ie0, y: mesh.edges.length };
        logV(1, `[BuildCraft] Girder id=${girder.id} nseg=${girder.nseg || 3} nodes=(${nA.id},${nB.id}) verts=${girder.pointRange.y - iv0}`);
    }

    // 3. Generate Ropes
    for (const rope of craft.ropes) {
        const pA = new Vec3(rope.nodeA.pos[0], rope.nodeA.pos[1], rope.nodeA.pos[2]);
        const pB = new Vec3(rope.nodeB.pos[0], rope.nodeB.pos[1], rope.nodeB.pos[2]);
        
        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;

        // Simple rope generator (just a line or rope primitive if available)
        mesh.rope(pA, new Vec3().setSub(pB, pA).normalize(), new Vec3().setSub(pB, pA).norm(), rope.nseg || 10);

        rope.pointRange = { x: iv0, y: mesh.verts.length };
        rope.stickRange = { x: ie0, y: mesh.edges.length };
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
                    return { pA, pB, sA: (b.nodeA.size || defaultNodeSize) * 0.5, sB: (b.nodeB.size || defaultNodeSize) * 0.5 };
                }
                return null;
            };

            const endsA = getEnds(boundA);
            const endsB = getEnds(boundB);
            if (!endsA || !endsB) continue;

            const lerp = (p0, p1, c) => new Vec3().setAddMul(p0, new Vec3().setSub(p1, p0), c);
            const getPerp = (v) => v.getSomePerp();

            const dirA = new Vec3().setSub(endsA.pB, endsA.pA);
            if (dirA.norm() === 0) continue;
            dirA.normalize();
            const dirB = new Vec3().setSub(endsB.pB, endsB.pA);
            const lenB = dirB.norm();
            if (lenB > 0) dirB.normalize();
            let side = (lenB === 0) ? getPerp(dirA) : new Vec3(dirB.x, dirB.y, dirB.z);
            let up = new Vec3().setCross(dirA, side);
            if (up.norm() === 0) { up = getPerp(dirA); } else { up.normalize(); }

            // Side directions for snapping to square corners
            const sideA = (lenB === 0) ? getPerp(dirA) : new Vec3(dirB.x, dirB.y, dirB.z);
            const sideB = new Vec3(dirA.x, dirA.y, dirA.z);
            if (sideA.norm() === 0) sideA.set(getPerp(dirA));
            else sideA.normalize();
            if (sideB.norm() === 0) sideB.set(getPerp(dirA));
            else sideB.normalize();

            const offA0 = ((plate.upA !== false) ? 1 : -1) * ((1 - spanA[0]) * endsA.sA + spanA[0] * endsA.sB);
            const offA1 = ((plate.upA !== false) ? 1 : -1) * ((1 - spanA[1]) * endsA.sA + spanA[1] * endsA.sB);
            const offB0 = ((plate.upB !== false) ? 1 : -1) * ((1 - spanB[0]) * endsB.sA + spanB[0] * endsB.sB);
            const offB1 = ((plate.upB !== false) ? 1 : -1) * ((1 - spanB[1]) * endsB.sA + spanB[1] * endsB.sB);
            const sideOff = (plate.sideOffset !== undefined) ? plate.sideOffset : 0;

            const a0 = lerp(endsA.pA, endsA.pB, spanA[0]).addMul(up, offA0).addMul(sideA, offA0 + sideOff);
            const a1 = lerp(endsA.pA, endsA.pB, spanA[1]).addMul(up, offA1).addMul(sideA, offA1 + sideOff);
            const b0 = lerp(endsB.pA, endsB.pB, spanB[0]).addMul(up, offB0).addMul(sideB, offB0 + sideOff);
            const b1 = lerp(endsB.pA, endsB.pB, spanB[1]).addMul(up, offB1).addMul(sideB, offB1 + sideOff);

            const iv0 = mesh.verts.length;
            const ie0 = mesh.edges.length;

            const nTop = plate.nx || 2;
            const nBottom = plate.ny || 2;
            const nRows = plate.nz || 1;

            // quad corners: a0-a1-b0-b1 (span order preserved; avoids needing to reverse spans)
            mesh.ParametricQuadPatch(nTop, nBottom, nRows, a0, a1, b0, b1);

            // reconstruct rowCounts (mirrors MeshesUV.ParametricQuadPatch) to weld
            const rowCounts = [];
            for (let i = 0; i < nRows; i++) {
                const t = (nRows === 1) ? 0 : i / (nRows - 1);
                const cnt = Math.round(nTop + (nBottom - nTop) * t);
                rowCounts.push(Math.max(1, cnt));
            }
            const rowOffsets = [];
            let acc = 0;
            for (const c of rowCounts) { rowOffsets.push(acc); acc += c; }
            const firstRow = { x: iv0 + rowOffsets[0], y: iv0 + rowOffsets[0] + rowCounts[0] };
            const lastRow  = { x: iv0 + rowOffsets[rowCounts.length - 1], y: iv0 + rowOffsets[rowCounts.length - 1] + rowCounts[rowCounts.length - 1] };
            // ParametricQuadPatch vertex order per row: left->right. Collect column edges (u=0 : a0->b0, u=1 : a1->b1)
            const edgeA_idxs = [];
            const edgeB_idxs = [];
            for (let r = 0; r < rowCounts.length; r++) {
                const base = iv0 + rowOffsets[r];
                edgeA_idxs.push(base);                          // first in row
                edgeB_idxs.push(base + rowCounts[r] - 1);       // last in row
            }
            const edgeA = { x: edgeA_idxs[0], y: edgeA_idxs[edgeA_idxs.length - 1] + 1 };
            const edgeB = { x: edgeB_idxs[0], y: edgeB_idxs[edgeB_idxs.length - 1] + 1 };

            // optional welding to girders/ropes
            if (plate.weldDist && plate.weldDist > 0) {
                const weldType = stickTypes.w;
                const listRange = (r) => { const arr = []; for (let i = r.x; i < r.y; i++) arr.push(i); return arr; };
                const w1 = mesh.weldListToRange(edgeA_idxs, boundA.pointRange, plate.weldDist, weldType);
                const w2 = mesh.weldListToRange(edgeB_idxs, boundB.pointRange, plate.weldDist, weldType);
                logV(1, `[BuildCraft] Weld1 a0->b0->A nb=${w1} plateEdgeA=${JSON.stringify(edgeA_idxs)} girderA=${JSON.stringify(listRange(boundA.pointRange))} R=${plate.weldDist}`);
                logV(1, `[BuildCraft] Weld2 a1->b1->B nb=${w2} plateEdgeB=${JSON.stringify(edgeB_idxs)} girderB=${JSON.stringify(listRange(boundB.pointRange))} R=${plate.weldDist}`);
            }

            plate.pointRange = { x: iv0, y: mesh.verts.length };
            plate.stickRange = { x: ie0, y: mesh.edges.length };
            logV(1, `[BuildCraft] ${plate.kind || 'Plate'} id=${plate.id} boundA=${plate.boundA.id} boundB=${plate.boundB.id} spanA=${spanA} spanB=${spanB} nx=${nTop} ny=${nBottom} nz=${nRows} upA=${plate.upA !== false} upB=${plate.upB !== false} sideOff=${sideOff} weld=${plate.weldDist || 0} a0=${a0.x.toFixed(3)},${a0.y.toFixed(3)},${a0.z.toFixed(3)} a1=${a1.x.toFixed(3)},${a1.y.toFixed(3)},${a1.z.toFixed(3)} b0=${b0.x.toFixed(3)},${b0.y.toFixed(3)},${b0.z.toFixed(3)} b1=${b1.x.toFixed(3)},${b1.y.toFixed(3)},${b1.z.toFixed(3)} firstRow=[${firstRow.x}-${firstRow.y}) lastRow=[${lastRow.x}-${lastRow.y}) edgeA=[${edgeA.x}-${edgeA.y}) edgeB=[${edgeB.x}-${edgeB.y}) verts=${plate.pointRange.y - iv0}`);
        }
    }

    // 5. Generate Rings
    for (const ring of craft.rings) {
        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;
        
        // Using mesh.wheel from MeshGenerators
        const pos = new Vec3(ring.pos[0], ring.pos[1], ring.pos[2]);
        const p1 = new Vec3(pos.x, pos.y + ring.R, pos.z); // axis hint
        const ax = new Vec3(0, 0, 1);
        mesh.wheel(pos, p1, ax, ring.nseg, ring.wh, stickTypes);

        ring.pointRange = { x: iv0, y: mesh.verts.length };
        ring.stickRange = { x: ie0, y: mesh.edges.length };
        logV(1, `[BuildCraft] Ring id=${ring.id} pos=${ring.pos} R=${ring.R} nseg=${ring.nseg} verts=${ring.pointRange.y - iv0}`);
    }

    // 6. Generate Slider Anchors and Visualization
    const sliderPathStickType = 2; // Different material for visualization
    for (const slider of craft.sliders) {
        const boundTo = slider.boundTo;
        if (!boundTo || boundTo.pointRange.x < 0) continue;

        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;

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

        slider.pointRange = { x: iv0, y: mesh.verts.length };
        slider.stickRange = { x: ie0, y: mesh.edges.length };
        logV(1, `[BuildCraft] Slider id=${slider.id} boundTo=${slider.boundTo.id} calong=${slider.calong} verts=${slider.pointRange.y - iv0}`);
    }
}
