import { Vec3 } from '../../common_js/Vec3.js';
import { Girder, Ring } from './SpaceCraft.js';
import { SDF_Cylinder } from '../../common_js/SDfuncs.js';

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

    // 6. Generate Sliders (rail + sliding vertex)
    const sliderPathStickType = 2; // Different material for visualization

    const getRailEndpoints = (rail) => {
        if (rail instanceof Girder) {
            const pA = new Vec3(rail.nodeA.pos[0], rail.nodeA.pos[1], rail.nodeA.pos[2]);
            const pB = new Vec3(rail.nodeB.pos[0], rail.nodeB.pos[1], rail.nodeB.pos[2]);
            return { pA, pB, closed: false };
        } else if (rail instanceof Ring) {
            const pos = new Vec3(rail.pos[0], rail.pos[1], rail.pos[2]);
            const p1 = new Vec3(pos.x, pos.y + rail.R, pos.z);
            const dir = new Vec3().setSub(p1, pos); // axis hint
            return { pA: pos, pB: pos, closed: true, ringDir: dir };
        }
        return null;
    };

    const sortPathByProjection = (psIdx, origin, dir) => {
        const d = new Vec3(dir.x, dir.y, dir.z).normalize();
        return psIdx.slice().sort((a, b) => {
            const pa = mesh.verts[a].pos;
            const pb = mesh.verts[b].pos;
            const ta = (pa.x - origin.x) * d.x + (pa.y - origin.y) * d.y + (pa.z - origin.z) * d.z;
            const tb = (pb.x - origin.x) * d.x + (pb.y - origin.y) * d.y + (pb.z - origin.z) * d.z;
            return ta - tb;
        });
    };

    // corner selector for a girder side: side=0..3 => (+side,+up), (+side,-up), (-side,+up), (-side,-up)
    const cornerDirs = [
        { a: 1, b: 1 },
        { a: 1, b: -1 },
        { a: -1, b: 1 },
        { a: -1, b: -1 },
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

    const buildPathOnEdge = (rail, side, radius) => {
        if (!rail || rail.pointRange.x < 0) return { ps: [], closed: false };
        if (rail instanceof Ring) {
            // For ring, still use SDF around circumference (offset unused)
            const ep = getRailEndpoints(rail);
            const p0 = new Vec3(ep.pA.x, ep.pA.y, ep.pA.z);
            const p1 = new Vec3(p0.x, p0.y + (ep.ringDir ? ep.ringDir.norm() : rail.R), p0.z);
            const sdf = SDF_Cylinder(p0, p1, radius, true);
            mesh.setSelectionKind('vert'); mesh.applySelection([], 'vert', false);
            const nsel = mesh.selectVertsBySDF(sdf, 0.0, true);
            logV(1, `[SliderPath] Ring rail=${rail.id} side=${side} radius=${radius} selAll=${nsel}`);
            const psIdx = mesh.selection.vec.filter(i => i >= rail.pointRange.x && i < rail.pointRange.y);
            logV(1, `[SliderPath] Ring rail=${rail.id} filtered=${psIdx.length} range=[${rail.pointRange.x},${rail.pointRange.y}) verts=${JSON.stringify(psIdx)}`);
            const ps = sortPathByProjection(psIdx, p0, new Vec3().setSub(p1, p0));
            return { ps, closed: true };
        }
        const eo = edgeOffset(rail, side, radius);
        const { offset, ep } = eo;
        const p0 = new Vec3(ep.pA.x, ep.pA.y, ep.pA.z).add(offset);
        const p1 = new Vec3(ep.pB.x, ep.pB.y, ep.pB.z).add(offset);
        // radius is too large (0.35 selects all verts of 1.0x1.0 square girder). 
        // We need a very small radius to pick only one edge.
        const selRadius = 0.05; 
        const sdf = SDF_Cylinder(p0, p1, selRadius, true);

        mesh.setSelectionKind('vert');
        mesh.applySelection([], 'vert', false);
        const nsel = mesh.selectVertsBySDF(sdf, 0.05, true); // use small positive threshold
        const fmt = (v) => Number(v).toFixed(3);
        logV(1, `[SliderPath] Girder rail=${rail.id} side=${side} radius=${radius} selAll=${nsel} p0=(${fmt(p0.x)},${fmt(p0.y)},${fmt(p0.z)}) p1=(${fmt(p1.x)},${fmt(p1.y)},${fmt(p1.z)}) offset=(${fmt(offset.x)},${fmt(offset.y)},${fmt(offset.z)}) ldir=${eo.ldir.toFixed(3)} lup=${eo.lup} lside=${eo.lside}`);
        let psIdx = mesh.selection.vec.filter(i => i >= rail.pointRange.x && i < rail.pointRange.y);
        logV(1, `[SliderPath] Girder rail=${rail.id} filtered=${psIdx.length} range=[${rail.pointRange.x},${rail.pointRange.y}) verts=${JSON.stringify(psIdx)}`);
        
        if (psIdx.length === 0) {
            // fallback: use whole pointRange
            psIdx = [];
            for (let i = rail.pointRange.x; i < rail.pointRange.y; i++) psIdx.push(i);
            logV(1, `[SliderPath] Girder rail=${rail.id} FALLBACK full range count=${psIdx.length}`);
        }
        const ps = sortPathByProjection(psIdx, p0, new Vec3().setSub(p1, p0));
        return { ps, closed: false };
    };

    const buildPathStride = (rail, side, radius) => {
        // approximate edgePathVert(i,j): pick closest rail vertex to desired corner positions along girder segments
        if (!(rail instanceof Girder) || rail.pointRange.x < 0) return buildPathOnEdge(rail, side, radius);
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
            // find nearest rail vertex within pointRange to target
            let best = -1; let bestd2 = 1e30;
            for (let i = rail.pointRange.x; i < rail.pointRange.y; i++) {
                const v = mesh.verts[i].pos;
                const dx = v.x - target.x; const dy = v.y - target.y; const dz = v.z - target.z;
                const d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < bestd2) { bestd2 = d2; best = i; }
            }
            if (best >= 0) ps.push(best);
        }
        logV(1, `[SliderPath] Stride rail=${rail.id} side=${side} nseg=${nseg} psN=${ps.length} range=[${rail.pointRange.x},${rail.pointRange.y}) offset=(${offset.x.toFixed(3)},${offset.y.toFixed(3)},${offset.z.toFixed(3)}) verts=${JSON.stringify(ps)}`);
        if (ps.length === 0) {
            for (let i = rail.pointRange.x; i < rail.pointRange.y; i++) ps.push(i);
            logV(1, `[SliderPath] Stride FALLBACK rail=${rail.id} using full range count=${ps.length}`);
        }
        return { ps, closed: false };
    };

    const findNearestOnPath = (path, p) => {
        if (path.ps.length === 0) return { cur: 0, idx: -1 };
        let best = -1; let bestd2 = 1e30;
        for (let i = 0; i < path.ps.length; i++) {
            const vp = mesh.verts[path.ps[i]].pos;
            const dx = vp.x - p.x; const dy = vp.y - p.y; const dz = vp.z - p.z;
            const d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < bestd2) { bestd2 = d2; best = i; }
        }
        return { cur: best, idx: best };
    };

    for (const slider of craft.sliders) {
        const rail = slider.rail;
        if (!rail || rail.pointRange.x < 0) { logV(1, `[BuildCraft] Slider id=${slider.id} rail missing range`); continue; }

        const slidingComp = slider.sliding || rail; // fallback: rail itself
        const slidingRange = slidingComp.pointRange || { x: -1, y: -1 };
        const ivSlide = (slider.slidingVertId >= 0 && slider.slidingVertId < mesh.verts.length)
            ? slider.slidingVertId
            : (slidingRange.x >= 0 ? slidingRange.x : -1);
        if (ivSlide < 0) continue;

        const useStride = !!slider.methodFlag;
        const radius = slider.pathRadius || 0.35;
        let path = useStride ? buildPathStride(rail, slider.side || 0, radius) : buildPathOnEdge(rail, slider.side || 0, radius);
        if (path.ps.length < 2) { // fallback: try SDF if stride failed, or stride if SDF failed
            const alt = useStride ? buildPathOnEdge(rail, slider.side || 0, radius * 1.5) : buildPathStride(rail, slider.side || 0, radius * 1.5);
            if (alt.ps.length > path.ps.length) path = alt;
        }
        if (path.ps.length < 2) { logV(1, `[BuildCraft] Slider id=${slider.id} path too short rail=${rail.id} side=${slider.side} method=${useStride ? 'stride' : 'sdf'} selN=${path.ps.length}`); continue; }
        path.closed = !!path.closed;

        // initialize cur from calong if provided, else nearest
        const pSlide = mesh.verts[ivSlide].pos;
        const tRaw = (slider.calong !== undefined && slider.calong !== null && Number.isFinite(slider.calong))
            ? slider.calong * (path.ps.length - 1)
            : null;
        let i0, tFrac;
        if (tRaw !== null) {
            const tClamped = Math.max(0, Math.min(path.ps.length - 1 - 1e-6, tRaw));
            i0 = Math.floor(tClamped);
            tFrac = tClamped - i0;
        } else {
            const nearest = findNearestOnPath(path, pSlide);
            i0 = Math.max(0, Math.min(path.ps.length - 2, nearest.idx));
            tFrac = 0.0;
        }
        path.cur = (i0 || 0) + tFrac;
        slider.path = path;

        const iv0 = mesh.verts.length;
        const ie0 = mesh.edges.length;

        // Create explicit slider vertex if different from sliding vertex (visualization). Otherwise reuse.
        const ivert = ivSlide;
        slider.ivert = ivert;

        // weld slider vertex to the two path verts around cur
        const iA = path.ps[i0];
        const iB = path.ps[(i0 + 1) % path.ps.length];
        const weldR = slider.weldDist || 0.25;
        mesh.weldListToRange([ivert], { x: iA, y: iA + 1 }, weldR, stickTypes.w);
        mesh.weldListToRange([ivert], { x: iB, y: iB + 1 }, weldR, stickTypes.w);

        // Visual link from slider vertex to interpolated position on path (cyan line)
        const pA = mesh.verts[iA].pos;
        const pB = mesh.verts[iB].pos;
        const pInterp = new Vec3().setV(pA).addMul(new Vec3().setSub(pB, pA), tFrac || 0.0);
        const ivInterp = mesh.vert(pInterp);
        const linkStickType = 2; // cyan-ish
        mesh.edge(ivert, ivInterp, linkStickType);

        // Visualization of the whole path (green line)
        const sliderPathStickType = 5; // green
        for (let i = 0; i < path.ps.length - (path.closed ? 0 : 1); i++) {
            const ia = path.ps[i];
            const ib = path.ps[(i + 1) % path.ps.length];
            mesh.edge(ia, ib, sliderPathStickType);
        }

        slider.pointRange = { x: iv0, y: mesh.verts.length };
        slider.stickRange = { x: ie0, y: mesh.edges.length };
        logV(1, `[BuildCraft] Slider id=${slider.id} rail=${rail.id} sliding=${slidingComp.id ?? 'self'} pathN=${path.ps.length} cur=${path.cur.toFixed(3)} method=${useStride ? 'stride' : 'sdf'} weldR=${weldR} verts=${slider.pointRange.y - iv0} ivSlide=${ivSlide} iA=${iA} iB=${iB} tFrac=${(tFrac||0).toFixed(3)}`);
    }
}
