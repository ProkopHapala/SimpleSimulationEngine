import { Vec3 } from '../../common_js/Vec3.js';
import { Girder, Ring } from './SpaceCraft.js';
import { SDF_Cylinder } from '../../common_js/SDfuncs.js';
import { logger } from '../../common_js/Logger.js';

const getRailEndpoints = (rail) => {
    if (rail instanceof Girder) {
        const pA = new Vec3(rail.nodeA.pos[0], rail.nodeA.pos[1], rail.nodeA.pos[2]);
        const pB = new Vec3(rail.nodeB.pos[0], rail.nodeB.pos[1], rail.nodeB.pos[2]);
        return { pA, pB, closed: false };
    } else if (rail instanceof Ring) {
        const pos = new Vec3(rail.pos[0] ?? rail.pos.x ?? 0, rail.pos[1] ?? rail.pos.y ?? 0, rail.pos[2] ?? rail.pos.z ?? 0);
        const rDir = new Vec3(rail.dir[0] ?? rail.dir.x ?? 1, rail.dir[1] ?? rail.dir.y ?? 0, rail.dir[2] ?? rail.dir.z ?? 0);
        const p1 = new Vec3().setV(pos).addMul(rDir, rail.R);
        const dir = new Vec3().setSub(p1, pos); // radius orientation hint
        return { pA: pos, pB: pos, closed: true, ringDir: dir };
    }
    return null;
};

const sortPathByProjection = (mesh, psIdx, origin, dir) => {
    const d = new Vec3(dir.x, dir.y, dir.z).normalize();
    return psIdx.slice().sort((a, b) => {
        const pa = mesh.verts[a].pos;
        const pb = mesh.verts[b].pos;
        const ta = (pa.x - origin.x) * d.x + (pa.y - origin.y) * d.y + (pa.z - origin.z) * d.z;
        const tb = (pb.x - origin.x) * d.x + (pb.y - origin.y) * d.y + (pb.z - origin.z) * d.z;
        return ta - tb;
    });
};

const cornerDirs = [
    { a: 1, b: 1 }, { a: 1, b: -1 }, { a: -1, b: 1 }, { a: -1, b: -1 },
];

const edgeOffset = (rail, side, offMag) => {
    const ep = getRailEndpoints(rail);
    const dir = rail instanceof Ring ? ep.ringDir : new Vec3().setSub(ep.pB, ep.pA);
    const dirN = new Vec3(dir.x, dir.y, dir.z);
    const ldir = dirN.normalize();
    if (ldir === 0) { dirN.set(1, 0, 0); }
    let offset = new Vec3();
    let sideVec = new Vec3();
    if (rail instanceof Girder) {
        const railUp = rail.up || [0, 1, 0];
        const up = new Vec3(railUp[0] ?? railUp.x ?? 0, railUp[1] ?? railUp.y ?? 1, railUp[2] ?? railUp.z ?? 0);
        let lup = up.normalize();
        if (lup < 1e-6) { up.set(0, 1, 0); lup = 1.0; }
        sideVec = new Vec3().setCross(dirN, up);
        let lside = sideVec.normalize();
        if (lside < 1e-6) { sideVec.set(0, 0, 1); lside = 1.0; }
        const cd = cornerDirs[side % 4];
        offset = new Vec3().setV(up).mulScalar(cd.b * offMag).addMul(sideVec, cd.a * offMag);
    }
    return { offset, dirN, ep, ldir: ldir || dirN.norm(), lup: rail instanceof Girder ? 1.0 : 0, lside: rail instanceof Girder ? 1.0 : 0 };
};

const buildPathByStrides = (rail, side) => {
    if (!rail || rail.pointRange.x < 0) {
        logger.warn(`[SpaceCraft2Mesh] buildPathByStrides: rail is null or has invalid pointRange`);
        return { ps: [], closed: false };
    }
    const ps = rail.sideToPath(side);
    const closed = (rail instanceof Ring);
    logger.info(`[SpaceCraft2Mesh] buildPathByStrides: rail=${rail.constructor.name} side=${side} n=${ps.length} closed=${closed}`);
    return { ps, closed };
};

const buildPathOnEdge = (mesh, rail, side, radius, useStrides = false) => {
    logger.info(`[SpaceCraft2Mesh] buildPathOnEdge: rail=${rail?.constructor.name} side=${side} useStrides=${useStrides}`);
    if (useStrides) return buildPathByStrides(rail, side);
    if (!rail || rail.pointRange.x < 0) return { ps: [], closed: false };
    if (rail instanceof Ring) {
        return buildPathByStrides(rail, side); // Default to strides for Ring if SDF fails or as fallback
    }
    const eo = edgeOffset(rail, side, radius);
    const { offset, ep } = eo;
    const p0 = new Vec3(ep.pA.x, ep.pA.y, ep.pA.z).add(offset);
    const p1 = new Vec3(ep.pB.x, ep.pB.y, ep.pB.z).add(offset);
    const selRadius = 0.05; 
    const sdf = SDF_Cylinder(p0, p1, selRadius, true);
    mesh.setSelectionKind('vert'); mesh.applySelection([], 'vert', false);
    mesh.selectVertsBySDF(sdf, 0.05, true);
    let psIdx = mesh.selection.vec.filter(i => i >= rail.pointRange.x && i < rail.pointRange.y);
    if (psIdx.length === 0) {
        psIdx = []; for (let i = rail.pointRange.x; i < rail.pointRange.y; i++) psIdx.push(i);
    }
    const ps = sortPathByProjection(mesh, psIdx, p0, new Vec3().setSub(p1, p0));
    return { ps, closed: false };
};

const buildPathStride = (mesh, rail, side, radius) => {
    if (!rail || rail.pointRange.x < 0) return { ps: [], closed: false };
    const ps = rail.sideToPath(side);
    const closed = (rail instanceof Ring);
    return { ps, closed };
};

/**
 * BuildCraft_blocks_js
 * A modern block-based spacecraft generator for JS.
 * Replaces each Node with a block (cube) and each Girder with a bridged truss.
 */
export function BuildCraft_blocks_js(mesh, craft) {
    mesh.clear();

    const logV = (lvl, msg) => {
        try { if (typeof logger !== 'undefined' && logger.uiVerbosity >= lvl) logger.info(msg); } catch (e) {}
    };

    const defaultNodeSize = 1.0;
    const stickTypes = { x: 1, y: 1, z: 1, w: 1 };
    const bridgeStickTypes = { x: 1, y: 1, z: 1, w: 1 };

    // 1. Generate Node Blocks
    for (const node of craft.nodes) {
        const p = new Vec3(node.pos[0], node.pos[1], node.pos[2]);
        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;
        const chunkRange = mesh.addCube(p, node.size || defaultNodeSize, true);
        node.pointRange = { x: iv0, y: mesh.verts.length };
        node.stickRange = { x: ie0, y: mesh.edges.length };
        node.chunkRange = chunkRange;
    }

    // 2. Generate Girder Bridges
    for (const girder of craft.girders) {
        const nA = girder.nodeA; const nB = girder.nodeB;
        const pA = new Vec3(nA.pos[0], nA.pos[1], nA.pos[2]);
        const pB = new Vec3(nB.pos[0], nB.pos[1], nB.pos[2]);
        const iv0 = mesh.verts.length; const ie0 = mesh.edges.length;
        mesh.bridgeFacingPolygons(pA, pB, nA.chunkRange, nB.chunkRange, girder.nseg || 3, bridgeStickTypes, bridgeStickTypes);
        girder.pointRange = { x: iv0, y: mesh.verts.length };
        girder.stickRange = { x: ie0, y: mesh.edges.length };
        logger.info(`[BuildCraft] Girder ${girder.id}: pointRange=[${girder.pointRange.x}, ${girder.pointRange.y}] nverts=${girder.pointRange.y - girder.pointRange.x}`);
    }

    // 3. Generate Ropes
    for (const rope of craft.ropes) {
        const pA = new Vec3(rope.nodeA.pos[0], rope.nodeA.pos[1], rope.nodeA.pos[2]);
        const pB = new Vec3(rope.nodeB.pos[0], rope.nodeB.pos[1], rope.nodeB.pos[2]);
        const iv0 = mesh.verts.length; const ie0 = mesh.edges.length;
        mesh.rope(pA, new Vec3().setSub(pB, pA).normalize(), new Vec3().setSub(pB, pA).norm(), rope.nseg || 10, rope.type);
        rope.pointRange = { x: iv0, y: mesh.verts.length };
        rope.stickRange = { x: ie0, y: mesh.edges.length };
    }

    // 4. Generate Plates
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
            const endsA = getEnds(boundA); const endsB = getEnds(boundB);
            if (!endsA || !endsB) continue;
            const lerp = (p0, p1, c) => new Vec3().setAddMul(p0, new Vec3().setSub(p1, p0), c);
            const getPerp = (v) => v.getSomePerp();
            const dirA = new Vec3().setSub(endsA.pB, endsA.pA); if (dirA.norm() === 0) continue; dirA.normalize();
            const dirB = new Vec3().setSub(endsB.pB, endsB.pA); const lenB = dirB.norm(); if (lenB > 0) dirB.normalize();
            let up = new Vec3().setCross(dirA, (lenB === 0) ? getPerp(dirA) : new Vec3(dirB.x, dirB.y, dirB.z));
            if (up.norm() === 0) { up = getPerp(dirA); } else { up.normalize(); }
            const sideA = (lenB === 0) ? getPerp(dirA) : new Vec3(dirB.x, dirB.y, dirB.z);
            const sideB = new Vec3(dirA.x, dirA.y, dirA.z);
            if (sideA.norm() === 0) sideA.set(getPerp(dirA)); else sideA.normalize();
            if (sideB.norm() === 0) sideB.set(getPerp(dirA)); else sideB.normalize();
            const offA0 = ((plate.upA !== false) ? 1 : -1) * ((1 - spanA[0]) * endsA.sA + spanA[0] * endsA.sB);
            const offA1 = ((plate.upA !== false) ? 1 : -1) * ((1 - spanA[1]) * endsA.sA + spanA[1] * endsA.sB);
            const offB0 = ((plate.upB !== false) ? 1 : -1) * ((1 - spanB[0]) * endsB.sA + spanB[0] * endsB.sB);
            const offB1 = ((plate.upB !== false) ? 1 : -1) * ((1 - spanB[1]) * endsB.sA + spanB[1] * endsB.sB);
            const sideOff = (plate.sideOffset !== undefined) ? plate.sideOffset : 0;
            const a0 = lerp(endsA.pA, endsA.pB, spanA[0]).addMul(up, offA0).addMul(sideA, offA0 + sideOff);
            const a1 = lerp(endsA.pA, endsA.pB, spanA[1]).addMul(up, offA1).addMul(sideA, offA1 + sideOff);
            const b0 = lerp(endsB.pA, endsB.pB, spanB[0]).addMul(up, offB0).addMul(sideB, offB0 + sideOff);
            const b1 = lerp(endsB.pA, endsB.pB, spanB[1]).addMul(up, offB1).addMul(sideB, offB1 + sideOff);
            const iv0 = mesh.verts.length; const ie0 = mesh.edges.length;
            mesh.ParametricQuadPatch(plate.nx || 2, plate.ny || 2, plate.nz || 1, a0, a1, b0, b1);
            if (plate.weldDist && plate.weldDist > 0) {
                const rowCounts = []; for (let i = 0; i < (plate.nz || 1); i++) {
                    const t = ((plate.nz || 1) === 1) ? 0 : i / ((plate.nz || 1) - 1);
                    rowCounts.push(Math.max(1, Math.round((plate.nx || 2) + ((plate.ny || 2) - (plate.nx || 2)) * t)));
                }
                const ea_idx = [], eb_idx = []; let acc = 0;
                for (const c of rowCounts) { ea_idx.push(iv0 + acc); eb_idx.push(iv0 + acc + c - 1); acc += c; }
                mesh.weldListToRange(ea_idx, boundA.pointRange, plate.weldDist, stickTypes.w);
                mesh.weldListToRange(eb_idx, boundB.pointRange, plate.weldDist, stickTypes.w);
            }
            plate.pointRange = { x: iv0, y: mesh.verts.length }; plate.stickRange = { x: ie0, y: mesh.edges.length };
        }
    }

    // 5. Generate Rings (hull part)
    for (let i = 0; i < craft.rings.length; i++) {
        const ring = craft.rings[i];
        const iv0 = mesh.verts.length; const ie0 = mesh.edges.length;
        const center = new Vec3().setV(ring.pos);
        let ax = new Vec3().setV(ring.dir);
        const lax = ax.normalize();
        if (lax < 1e-9) { ax = new Vec3(0, 0, 1); }
        const up = ring.up ? new Vec3().setV(ring.up) : null;
        const R = Number(ring.R ?? 5.0);
        if (!(R > 1e-9)) {
            logger.error(`[BuildCraft] Ring ${ring.id} has non-positive radius R=${R}, skipping mesh.`);
            continue;
        }
        const rimDir = up ? new Vec3().setV(up) : ax.getSomeOrtho();
        rimDir.makeOrtho(ax);
        const lrd = rimDir.normalize();
        if (lrd < 1e-9) {
            logger.error(`[BuildCraft] Ring ${ring.id} failed to find rim direction, skipping mesh.`);
            continue;
        }
        const p1 = new Vec3().setV(center).addMul(rimDir, R);
        const phase = Number(ring.phase ?? 0.0);
        mesh.wheel(center, p1, ax, Number(ring.nseg ?? 16), { x: Number(ring.wh?.x ?? ring.wh?.[0] ?? 0.4), y: Number(ring.wh?.y ?? ring.wh?.[1] ?? 0.4) }, { x: 1, y: 1, z: 1, w: 1 }, phase);
        ring.pointRange = { x: iv0, y: mesh.verts.length }; ring.stickRange = { x: ie0, y: mesh.edges.length };
        logger.info(`[BuildCraft] Ring ${ring.id}: pointRange=[${ring.pointRange.x}, ${ring.pointRange.y}] nverts=${ring.pointRange.y - ring.pointRange.x} phase=${phase.toFixed(3)}`);
    }

    // 6. Generate Sliders (hull part - joining logic)
    for (const slider of craft.sliders) {
        const path = slider.path;
        const rail = path?.rail;
        if (!rail || rail.pointRange.x < 0) {
            logger.warn(`[BuildCraft] Slider ${slider.id} rail is invalid or mesh not built yet`);
            continue;
        }

        // If slider has a pAttach, we create a new vertex on the girder and weld it
        // This is used by RingAttached to ensure the slider is exactly on the path point
        if (slider.pAttach) {
            const pa = new Vec3().setV(slider.pAttach);
            if (isFinite(pa.x) && isFinite(pa.y) && isFinite(pa.z)) {
                // 1. Create a NEW vertex at the exact interpolation point
                const vNew = mesh.vert(pa);
                slider.ivert = vNew;
                
                // 2. Weld this new vertex to the sliding component (usually a girder)
                // This ensures the ring is physically connected to the girder in the simulation.
                const slidingComp = slider.sliding || rail;
                if (slidingComp.pointRange && slidingComp.pointRange.x >= 0) {
                    mesh.vert_weld(pa, slider.weldDist || 0.25, slidingComp.pointRange.x, slidingComp.pointRange.y, stickTypes.w);
                }
            }
        } else {
            // Standard sliding logic using existing vertices
            const slidingComp = slider.sliding || rail;
            const slidingRange = slidingComp.pointRange || { x: -1, y: -1 };
            const pSlide = (Number.isInteger(slider.slidingVertId) && slider.slidingVertId >= 0 && slider.slidingVertId < mesh.verts.length) 
                ? mesh.verts[slider.slidingVertId].pos 
                : (Number.isInteger(slidingRange.x) && slidingRange.x >= 0 && slidingRange.x < mesh.verts.length ? mesh.verts[slidingRange.x].pos : null);
            
            if (!pSlide) continue;
            
            const weldR = slider.weldDist || 0.25;
            const nearestIdx = mesh.vert_collapse(pSlide, weldR, rail.pointRange.x, rail.pointRange.y);
            slider.ivert = (nearestIdx >= 0) ? nearestIdx : ((slider.slidingVertId >= 0) ? slider.slidingVertId : slidingRange.x);
        }

        // Build/Update Path indices if they are empty
        // This allows multiple sliders to share the same pre-computed path indices
        if (path.ps.length === 0) {
            const pData = buildPathOnEdge(mesh, rail, path.side, path.radius || 0.35, path.methodFlag);
            path.ps = pData.ps;
            path.closed = pData.closed;
        }
        
        logger.info(`[BuildCraft] Slider ${slider.id}: rail=${rail.constructor.name} pathN=${path.ps.length} closed=${path.closed}`);
    }
}

/**
 * BuildCraft_aux_js
 * Generates auxiliary/debug visualization (slider paths, anchor lines).
 */
export function BuildCraft_aux_js(auxMesh, craft, hullMesh) {
    auxMesh.clear();
    const sliderPathStickType = 5; // Cyan
    const anchorStickType = 2;     // Green
    // Draw all unique paths using cubic B-spline for smoothness
    for (const path of craft.paths) {
        if (path.ps.length < 2) continue;
        const rail = path.rail;
        if (!rail || rail.pointRange.x < 0) continue;

        // Number of intermediate segments to draw for smoothness
        const drawSegments = path.ps.length * 4; 
        let pPrev = path.interpolate(0, hullMesh.verts);
        for (let j = 1; j <= drawSegments; j++) {
            const t = j / drawSegments;
            const pCurr = path.interpolate(t, hullMesh.verts);
            auxMesh.edge(auxMesh.vert(pPrev), auxMesh.vert(pCurr), sliderPathStickType);
            pPrev = pCurr;
        }
    }

    for (const slider of craft.sliders) {
        const path = slider.path;
        if (!path || path.ps.length < 2) continue;
        
        // Draw Anchor using smooth interpolation
        const tRaw = (slider.calong !== undefined && slider.calong !== null) ? slider.calong : 0.5;
        const pInterp = path.interpolate(tRaw, hullMesh.verts);
        
        if (Number.isInteger(slider.ivert) && slider.ivert >= 0 && slider.ivert < hullMesh.verts.length) {
            auxMesh.edge(auxMesh.vert(hullMesh.verts[slider.ivert].pos), auxMesh.vert(pInterp), anchorStickType);
        }
    }
}
