import { Vec3 } from '../../common_js/Vec3.js';
import { Girder, Ring } from './SpaceCraft.js';
import { SDF_Cylinder } from '../../common_js/SDfuncs.js';

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
    if (!rail || rail.pointRange.x < 0) return { ps: [], closed: false };
    const ps = [];
    const n = rail.pointRange.y - rail.pointRange.x;
    if (rail instanceof Girder) {
        const offset = side % 4;
        for (let i = rail.pointRange.x + offset; i < rail.pointRange.y; i += 4) ps.push(i);
        return { ps, closed: false };
    } else if (rail instanceof Ring) {
        const offset = side % 4;
        for (let i = rail.pointRange.x + offset; i < rail.pointRange.y; i += 4) ps.push(i);
        return { ps, closed: true };
    }
    return { ps: [], closed: false };
};

const buildPathOnEdge = (mesh, rail, side, radius, useStrides = false) => {
    if (useStrides) return buildPathByStrides(rail, side);
    if (!rail || rail.pointRange.x < 0) return { ps: [], closed: false };
    if (rail instanceof Ring) {
        const ep = getRailEndpoints(rail);
        const p0 = new Vec3(ep.pA.x, ep.pA.y, ep.pA.z);
        const p1 = new Vec3(p0.x, p0.y + (ep.ringDir ? ep.ringDir.norm() : rail.R), p0.z);
        const sdf = SDF_Cylinder(p0, p1, radius, true);
        mesh.setSelectionKind('vert'); mesh.applySelection([], 'vert', false);
        mesh.selectVertsBySDF(sdf, 0.0, true);
        const psIdx = mesh.selection.vec.filter(i => i >= rail.pointRange.x && i < rail.pointRange.y);
        const ps = sortPathByProjection(mesh, psIdx, p0, new Vec3().setSub(p1, p0));
        return { ps, closed: true };
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
    if (!(rail instanceof Girder) || rail.pointRange.x < 0) return buildPathOnEdge(mesh, rail, side, radius);
    const nseg = rail.nseg || Math.max(1, (rail.pointRange.y - rail.pointRange.x) - 1);
    const eo = edgeOffset(rail, side, radius);
    const { offset, ep } = eo;
    const pA = new Vec3(ep.pA.x, ep.pA.y, ep.pA.z).add(offset);
    const pB = new Vec3(ep.pB.x, ep.pB.y, ep.pB.z).add(offset);
    const dir = new Vec3().setSub(pB, pA);
    const ps = [];
    for (let j = 0; j <= nseg; j++) {
        const t = j / nseg;
        const target = new Vec3(pA.x, pA.y, pA.z).addMul(dir, t);
        let best = -1; let bestd2 = 1e30;
        for (let i = rail.pointRange.x; i < rail.pointRange.y; i++) {
            const v = mesh.verts[i].pos;
            const dx = v.x - target.x; const dy = v.y - target.y; const dz = v.z - target.z;
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < bestd2) { bestd2 = d2; best = i; }
        }
        if (best >= 0) ps.push(best);
    }
    return { ps, closed: false };
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
    }

    // 3. Generate Ropes
    for (const rope of craft.ropes) {
        const pA = new Vec3(rope.nodeA.pos[0], rope.nodeA.pos[1], rope.nodeA.pos[2]);
        const pB = new Vec3(rope.nodeB.pos[0], rope.nodeB.pos[1], rope.nodeB.pos[2]);
        const iv0 = mesh.verts.length; const ie0 = mesh.edges.length;
        mesh.rope(pA, new Vec3().setSub(pB, pA).normalize(), new Vec3().setSub(pB, pA).norm(), rope.nseg || 10);
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
        const center = new Vec3(Number(ring.pos?.[0] ?? ring.pos?.x ?? 0), Number(ring.pos?.[1] ?? ring.pos?.y ?? 0), Number(ring.pos?.[2] ?? ring.pos?.z ?? 0));
        const up = new Vec3(Number(ring.up?.[0] ?? ring.up?.x ?? 0), Number(ring.up?.[1] ?? ring.up?.y ?? 0), Number(ring.up?.[2] ?? ring.up?.z ?? 1));
        const dir = new Vec3(Number(ring.dir?.[0] ?? ring.dir?.x ?? 1), Number(ring.dir?.[1] ?? ring.dir?.y ?? 0), Number(ring.dir?.[2] ?? ring.dir?.z ?? 0));
        const p1 = new Vec3().setV(center).addMul(dir, Number(ring.R ?? 5.0));
        mesh.wheel(center, p1, up, Number(ring.nseg ?? 16), { x: Number(ring.wh?.x ?? ring.wh?.[0] ?? 0.4), y: Number(ring.wh?.y ?? ring.wh?.[1] ?? 0.4) }, { x: 5, y: 5, z: 5, w: 5 });
        ring.pointRange = { x: iv0, y: mesh.verts.length }; ring.stickRange = { x: ie0, y: mesh.edges.length };
    }

    // 6. Generate Sliders (hull part - joining logic)
    for (const slider of craft.sliders) {
        const rail = slider.rail; if (!rail || rail.pointRange.x < 0) continue;
        const slidingComp = slider.sliding || rail;
        const slidingRange = slidingComp.pointRange || { x: -1, y: -1 };
        const pSlide = (Number.isInteger(slider.slidingVertId) && slider.slidingVertId >= 0 && slider.slidingVertId < mesh.verts.length) 
            ? mesh.verts[slider.slidingVertId].pos 
            : (Number.isInteger(slidingRange.x) && slidingRange.x >= 0 && slidingRange.x < mesh.verts.length ? mesh.verts[slidingRange.x].pos : null);
        if (!pSlide) continue;
        const weldR = slider.weldDist || 0.25;
        if (slider.joinStrategy === 'weld') {
            const iv = (slider.slidingVertId >= 0) ? slider.slidingVertId : slidingRange.x;
            mesh.vert_weld(pSlide, weldR, rail.pointRange.x, rail.pointRange.y, stickTypes.w);
            slider.ivert = iv;
        } else {
            const nearestIdx = mesh.vert_collapse(pSlide, weldR, rail.pointRange.x, rail.pointRange.y);
            slider.ivert = (nearestIdx >= 0) ? nearestIdx : ((slider.slidingVertId >= 0) ? slider.slidingVertId : slidingRange.x);
        }
        // Pre-build path for Aux visualization to reuse
        slider.path = buildPathOnEdge(mesh, rail, slider.side, slider.pathRadius || 0.35, slider.methodFlag);
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
    for (const slider of craft.sliders) {
        const rail = slider.rail; if (!rail || rail.pointRange.x < 0) continue;
        const path = (slider.path && slider.path.ps) ? slider.path : { ps: [], closed: false };
        if (path.ps.length < 2) continue;
        // Draw Path
        for (let i = 0; i < path.ps.length - (path.closed ? 0 : 1); i++) {
            const iA = path.ps[i], iB = path.ps[(i + 1) % path.ps.length];
            if (Number.isInteger(iA) && Number.isInteger(iB) && iA >= 0 && iA < hullMesh.verts.length && iB >= 0 && iB < hullMesh.verts.length) {
                auxMesh.edge(auxMesh.vert(hullMesh.verts[iA].pos), auxMesh.vert(hullMesh.verts[iB].pos), sliderPathStickType);
            }
        }
        // Draw Anchor
        const tRaw = (slider.calong !== undefined && slider.calong !== null) ? slider.calong : 0.5;
        const tClamped = Math.max(0, Math.min(path.ps.length - 1 - 1e-6, tRaw * (path.ps.length - 1)));
        const i0 = Math.floor(tClamped), tFrac = tClamped - i0;
        const iA = path.ps[i0], iB = path.ps[(i0 + 1) % path.ps.length];
        
        if (Number.isInteger(iA) && Number.isInteger(iB) && iA >= 0 && iA < hullMesh.verts.length && iB >= 0 && iB < hullMesh.verts.length) {
            const pA = hullMesh.verts[iA].pos, pB = hullMesh.verts[iB].pos;
            const pInterp = new Vec3().setV(pA).addMul(new Vec3().setSub(pB, pA), tFrac);
            if (Number.isInteger(slider.ivert) && slider.ivert >= 0 && slider.ivert < hullMesh.verts.length) {
                auxMesh.edge(auxMesh.vert(hullMesh.verts[slider.ivert].pos), auxMesh.vert(pInterp), anchorStickType);
            }
        }
    }
}
