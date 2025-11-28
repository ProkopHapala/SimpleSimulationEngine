
class MeshBuilder {
    // DEPENDENCY: This class assumes a global 'logger' object exists.
    // In browser: window.logger
    // In Node.js: global.logger
    constructor() {

        this.verts = [];   // Vertices: Array of {pos: Vec3, nor: Vec3, uv: {x, y}}
        this.edges = [];   // Edges:    Array of {x: int, y: int, z: type, w: type2}
        this.chunks = [];  // Chunks:   Array of {x: stripStart, y: edgeStart, z: count, w: type}
        this.strips = [];  // Strips:   Flat array of vertex/edge indices for faces        
        this.blocks = [];  // Blocks:   [ivert_start, iedge_start, ichunk_start]

        // Edge adjacency tracking for bevel and other operations
        this.edgesOfVerts = null;  // Map from vertex index to edge indices

        // Temporary vectors (reused to avoid GC)
        this._tmp1 = new Vec3();
        this._tmp2 = new Vec3();
        this._tmp3 = new Vec3();

        // --- Selection state (JS analogue of Selection + SelectionBanks) ---
        // Use SelectionBanks from Selection.js to manage multiple banks.
        const nBanks = 8; // UI currently exposes 8 banks (0..7)
        this.selectionBanks = new SelectionBanks(nBanks);
        this.selection = this.selectionBanks.curSelection;
        this.selection.kind = 'vert';
        this.currentSelectionBank = this.selectionBanks.icurSelection;
    }

    clear() {
        this.verts = [];
        this.edges = [];
        this.chunks = [];
        this.strips = [];
        this.blocks = [];
        logger.info("MeshBuilder cleared.");
        // Clearing mesh also clears working selection
        this.clearSelection();
    }

    // --- Selection helpers ---

    clearSelection() {
        if (this.selection) {
            this.selection.clear();
        }
    }

    setSelectionKind(kind) {
        if (kind === 'vert' || kind === 'edge') {
            if (this.selection) {
                this.selection.kind = kind;
            }
        }
    }

    /**
     * Replace or augment current selection with given indices.
     * @param {Array<number>} indices
     * @param {string} kind 'vert' or 'edge'
     * @param {boolean} additive if true, keep existing selection and add; otherwise replace
     */
    applySelection(indices, kind = null, additive = false) {
        if (!this.selection) return;
        if (!additive) {
            this.selection.clear();
        }
        if (kind) {
            this.setSelectionKind(kind);
        }
        for (const i of indices) {
            if (Number.isInteger(i) && i >= 0) {
                this.selection.add(i);
            }
        }
    }

    /**
     * Subtract given indices from current selection.
     * @param {Array<number>} indices
     */
    subtractSelection(indices) {
        if (!this.selection) return;
        for (const i of indices) {
            if (Number.isInteger(i) && i >= 0) {
                this.selection.remove(i);
            }
        }
    }

    /**
     * Get current selection as a sorted array of indices.
     */
    getSelectionArray() {
        if (!this.selection) return [];
        return this.selection.toSortedArray();
    }

    /**
     * JS analogue of Mesh::Builder2::selectVertsBySDF.
     *
     * Select vertices for which the provided signed-distance-like function
     * `sdf(pos)` is below the given threshold. Stores the result into the
     * current Selection bank.
     *
     * @param {(pos: Vec3) => number} sdf       Distance / SDF function on vertex positions.
     * @param {number} threshold               Inclusive threshold (default 0.0).
     * @param {boolean} clear                  If true, clear existing selection first.
     * @returns {number}                       Number of selected vertices.
     */
    selectVertsBySDF(sdf, threshold = 0.0, clear = true) {
        if (!this.selection || typeof sdf !== 'function') return 0;
        if (clear) {
            this.selection.clear();
        }
        this.selection.kind = 'vert';

        // Mirror C++: Selection::selectByPredicate over the verts container.
        return this.selection.selectByPredicate(this.verts, (vert) => {
            return sdf(vert.pos) < threshold;
        });
    }

    /**
     * Subtract vertices that satisfy an SDF from the current selection.
     *
     * This is the complement of selectVertsBySDF: instead of replacing or
     * extending the selection, it removes all vertex indices i for which
     * sdf(verts[i].pos) < threshold.
     *
     * Internally this builds a temporary Set of indices and uses the
     * Selection.subtract method from Selection.js to do the actual removal.
     *
     * @param {(pos: Vec3) => number} sdf       Distance / SDF function on vertex positions.
     * @param {number} threshold               Inclusive threshold (default 0.0).
     * @returns {number}                       Number of vertices removed from selection.
     */
    subtractVertsBySDF(sdf, threshold = 0.0) {
        if (!this.selection || typeof sdf !== 'function') return 0;

        const toRemove = new Set();
        const n = this.verts.length;
        for (let i = 0; i < n; i++) {
            const v = this.verts[i];
            if (!v || !v.pos) continue;
            if (sdf(v.pos) < threshold) {
                toRemove.add(i);
            }
        }

        return this.selection.subtract(toRemove);
    }

    /**
     * JS analogue of Mesh::Builder2::selectEdgesBySDF.
     *
     * Select edges whose *both endpoints* satisfy sdf(pos) < threshold.
     * Writes into the current Selection bank with kind = 'edge'.
     *
     * @param {(pos: Vec3) => number} sdf       Distance / SDF function on vertex positions.
     * @param {number} threshold                Inclusive threshold (default 0.0).
     * @param {boolean} clear                   If true, clear existing selection first.
     * @returns {number}                        Number of selected edges.
     */
    selectEdgesBySDF(sdf, threshold = 0.0, clear = true) {
        if (!this.selection || typeof sdf !== 'function') return 0;
        if (clear) {
            this.selection.clear();
        }
        this.selection.kind = 'edge';

        return this.selection.selectByPredicate(this.edges, (edge) => {
            const pA = this.verts[edge.x].pos;
            const pB = this.verts[edge.y].pos;
            return (sdf(pA) < threshold) && (sdf(pB) < threshold);
        });
    }

    /**
     * Save current working selection into a bank index.
     */
    saveSelectionToBank(bankIndex) {
        if (!this.selectionBanks || !this.selection) return;
        const n = this.selectionBanks.selections.length;
        const i = Math.max(0, Math.min(n - 1, bankIndex | 0));
        const bank = this.selectionBanks.selections[i];
        bank.kind = this.selection.kind;
        bank.clear();
        // Copy current selection items into bank
        for (const id of this.selection.vec) {
            if (id >= 0) bank.add(id);
        }
        this.currentSelectionBank = i;
    }

    /**
     * Load selection from bank index into working selection.
     */
    loadSelectionFromBank(bankIndex) {
        if (!this.selectionBanks) return;
        const n = this.selectionBanks.selections.length;
        const i = Math.max(0, Math.min(n - 1, bankIndex | 0));
        const bank = this.selectionBanks.selections[i];
        this.selectionBanks.icurSelection = i;
        this.selectionBanks.curSelection = bank;
        this.selection = bank;
        this.currentSelectionBank = i;
    }

    /**
     * Clear a specific selection bank.
     */
    clearSelectionBank(bankIndex) {
        if (!this.selectionBanks) return;
        const n = this.selectionBanks.selections.length;
        const i = Math.max(0, Math.min(n - 1, bankIndex | 0));
        const bank = this.selectionBanks.selections[i];
        bank.clear();
    }

    // --- Basic Primitives ---

    vert(pos) {
        const x = pos.x !== undefined ? pos.x : pos[0];
        const y = pos.y !== undefined ? pos.y : pos[1];
        const z = pos.z !== undefined ? pos.z : pos[2];
        const v = {
            pos: new Vec3(x, y, z),
            nor: new Vec3(0, 0, 1),
            uv: { x: 0, y: 0 }
        };
        this.verts.push(v);
        const idx = this.verts.length - 1;
        if (logger.verb(4)) logger.debug(`Vert[${idx}]: ${v.pos.x}, ${v.pos.y}, ${v.pos.z}`);
        return idx;
    }

    edge(a, b, type = -1, type2 = 0) {
        this.edges.push({ x: a, y: b, z: type, w: type2 });
        const idx = this.edges.length - 1;
        if (logger.verb(4)) logger.debug(`Edge[${idx}]: ${a} -> ${b} (t=${type})`);
        return idx;
    }

    chunk(data) {
        // data = {x: stripStart, y: edgeStart, z: count, w: type}
        this.chunks.push(data);
        return this.chunks.length - 1;
    }

    block() {
        const b = {
            ivert: this.verts.length,
            iedge: this.edges.length,
            ichunk: this.chunks.length
        };
        this.blocks.push(b);
        return this.blocks.length - 1;
    }

    // --- Utility Functions ---

    getCOG(n, ivs) {
        const cog = new Vec3();
        for (let i = 0; i < n; i++) {
            cog.add(this.verts[ivs[i]].pos);
        }
        cog.mulScalar(1.0 / n);
        return cog;
    }

    getChunkCOG(ichunk) {
        const ch = this.chunks[ichunk];
        const n = ch.z;
        const i0 = ch.x;
        const ivs = [];
        for (let i = 0; i < n; i++) {
            ivs.push(this.strips[i0 + i]);
        }
        return this.getCOG(n, ivs);
    }

    polygonNormal(ichunk) {
        const ch = this.chunks[ichunk];
        const n = ch.z;
        const i0 = ch.x;
        const ivs = [];
        for (let i = 0; i < n; i++) {
            ivs.push(this.strips[i0 + i]);
        }

        const nrm = new Vec3(0, 0, 0);
        let a = this.verts[ivs[0]].pos;
        let b = this.verts[ivs[1]].pos;
        const tmp1 = new Vec3();
        const tmp2 = new Vec3();
        const tmpCross = new Vec3();

        for (let j = 1; j < n - 1; j++) {
            const c = this.verts[ivs[j + 1]].pos;
            // cross(b-a, c-a)
            tmp1.setSub(b, a);
            tmp2.setSub(c, a);
            tmpCross.setCross(tmp1, tmp2);
            nrm.add(tmpCross);
            a = b;
            b = c;
        }
        nrm.normalize();
        return nrm;
    }

    // --- Alignment Function ---

    /**
     * Align two polygons rotationally to minimize twist
     * Based on MeshBuilder2.cpp::alling_polygons
     */
    alling_polygons(n, ivs1, ivs2, ipiv = 0) {
        const cog1 = this.getCOG(n, ivs1);
        const cog2 = this.getCOG(n, ivs2);

        // Compute axis
        const ax = this._tmp1.setSub(cog2, cog1);
        ax.normalize();

        // Compute orthogonal basis
        const u = this._tmp2.setSub(this.verts[ivs1[ipiv]].pos, cog1);
        u.makeOrthoU(ax);
        u.normalize();

        const v = this._tmp3.setCross(ax, u);
        v.normalize();

        // Project points to UV plane
        const uv1 = new Array(n);
        const uv2 = new Array(n);

        for (let i = 0; i < n; i++) {
            const d1 = new Vec3().setSub(this.verts[ivs1[i]].pos, cog1);
            const d2 = new Vec3().setSub(this.verts[ivs2[i]].pos, cog2);

            uv1[i] = new Vec3(d1.dot(u), d1.dot(v), 0);
            uv1[i].normalize();

            uv2[i] = new Vec3(d2.dot(u), d2.dot(v), 0);
            uv2[i].normalize();
        }

        // Match points by nearest UV projection
        // Store original ivs2 values before remapping
        const originalIvs2 = ivs2.slice();

        for (let i = 0; i < n; i++) {
            let dmin = -Infinity; // FIX: Initialize to -Infinity, not -1.0
            let jbest = 0;
            const uvi = uv1[i];
            for (let j = 0; j < n; j++) {
                const d = uvi.x * uv2[j].x + uvi.y * uv2[j].y; // dot product
                if (d > dmin) {
                    dmin = d;
                    jbest = j;
                }
            }
            // Map ivs1[i] to the vertex that was originally at ivs2[jbest]
            ivs2[i] = originalIvs2[jbest];
            logger.debug(`alling_polygons: Map ${i} -> ${jbest} (d=${dmin})`);
        }
    }

    // --- Bridge Quads (Volumetric Truss Generation) ---

    /**
     * Bridge two quad faces with a volumetric truss
     * Based on MeshBuilder2.cpp::bridge_quads
     * 
     * @param {Object} q1 - First quad {x, y, z, w} (vertex indices)
     * @param {Object} q2 - Second quad {x, y, z, w} (vertex indices)
     * @param {number} nseg - Number of segments
     * @param {Object} stickTypes - {x: longitudinal, y: ring, z: spiral, w: internal}
     * @param {Object} mask - {x, y, z, w} controls which edges to generate
     * @param {boolean} bAlling - Whether to align quads rotationally
     */
    bridge_quads(q1, q2, nseg, stickTypes = { x: -1, y: -1, z: -1, w: -1 }, mask = { x: 1, y: 1, z: 1, w: 1 }, bAlling = true) {
        logger.debug(`bridge_quads: nseg=${nseg}, q1=[${q1.x},${q1.y},${q1.z},${q1.w}], q2=[${q2.x},${q2.y},${q2.z},${q2.w}]`);

        // Align polygons if requested
        if (bAlling) {
            const ivs1 = [q1.x, q1.y, q1.z, q1.w];
            const ivs2 = [q2.x, q2.y, q2.z, q2.w];
            this.alling_polygons(4, ivs1, ivs2, 0);
            q2.x = ivs2[0]; q2.y = ivs2[1]; q2.z = ivs2[2]; q2.w = ivs2[3];
        }

        // Get vertex positions
        const A1 = this.verts[q1.x].pos;
        const B1 = this.verts[q1.y].pos;
        const C1 = this.verts[q1.z].pos;
        const D1 = this.verts[q1.w].pos;

        const A2 = this.verts[q2.x].pos;
        const B2 = this.verts[q2.y].pos;
        const C2 = this.verts[q2.z].pos;
        const D2 = this.verts[q2.w].pos;

        const i00start = this.verts.length;
        const dc = 1.0 / nseg;

        let oA = q1.x, oB = q1.y, oC = q1.z, oD = q1.w;

        for (let i = 0; i < nseg; i++) {
            let iA, iB, iC, iD;

            if (i < (nseg - 1)) {
                // Create intermediate rings
                const c = (i + 1) * dc;
                const mc = 1 - c;

                iA = this.vert(new Vec3().setLincomb(mc, A1, c, A2));
                iB = this.vert(new Vec3().setLincomb(mc, B1, c, B2));
                iC = this.vert(new Vec3().setLincomb(mc, C1, c, C2));
                iD = this.vert(new Vec3().setLincomb(mc, D1, c, D2));

                // Ring edges
                this.edge(iA, iB, stickTypes.y);
                this.edge(iB, iC, stickTypes.y);
                this.edge(iC, iD, stickTypes.y);
                this.edge(iD, iA, stickTypes.y);
            } else {
                // Last segment connects to q2
                iA = q2.x;
                iB = q2.y;
                iC = q2.z;
                iD = q2.w;
            }

            // Longitudinal edges
            this.edge(oA, iA, stickTypes.x);
            this.edge(oB, iB, stickTypes.x);
            this.edge(oC, iC, stickTypes.x);
            this.edge(oD, iD, stickTypes.x);

            // Spiral edges (diagonal bracing)
            if (mask.x) {
                this.edge(oA, iB, stickTypes.z);
                this.edge(oB, iC, stickTypes.z);
                this.edge(oC, iD, stickTypes.z);
                this.edge(oD, iA, stickTypes.z);
            }
            if (mask.y) {
                this.edge(oB, iA, stickTypes.z);
                this.edge(oC, iB, stickTypes.z);
                this.edge(oD, iC, stickTypes.z);
                this.edge(oA, iD, stickTypes.z);
            }

            // Internal edges (cross-bracing)
            if (mask.z) {
                this.edge(oA, iC, stickTypes.z);
                this.edge(oB, iD, stickTypes.z);
            }
            if (mask.w) {
                this.edge(oC, iA, stickTypes.z);
                this.edge(oD, iB, stickTypes.z);
            }

            oA = iA; oB = iB; oC = iC; oD = iD;
        }

        return i00start;
    }

    findMostFacingNormal(hray, chrange, cosMin = 0.0, bTwoSide = false, distWeight = 0.0, ray0 = new Vec3(0, 0, 0)) {
        // chrange is {x: start, y: end} (exclusive end)
        const nch = chrange.y - chrange.x;
        // const chs = []; for(let i=0; i<nch; i++) chs.push(chrange.x + i);

        let ibest = -1;
        let cmax = -1.0;
        const bDist = Math.abs(distWeight) > 1e-100;

        logger.debug(`findMostFacingNormal: hray=${hray.toString()}, nch=${nch}, range=[${chrange.x}, ${chrange.y}]`);

        for (let i = 0; i < nch; i++) {
            const ich = chrange.x + i;
            const nr = this.polygonNormal(ich);
            let c = hray.dot(nr);

            if (bTwoSide) c = Math.abs(c);

            logger.debug(`  Chunk ${ich}: nr=${nr.toString()}, c=${c}`);

            if (c > cosMin) {
                if (bDist) {
                    const p = this.getChunkCOG(ich);
                    const r = new Vec3().setSub(p, ray0).norm();
                    c -= distWeight * r;
                }
                if (c > cmax) {
                    ibest = ich;
                    cmax = c;
                    logger.debug(`  New best: ${ibest} (c=${cmax})`);
                }
            }
        }
        return ibest;
    }

    bridgeFacingPolygons(p1, p2, chr1, chr2, nseg, stickTypes = { x: -1, y: -1, z: -1, w: -1 }, mask = { x: 1, y: 1, z: 1, w: 1 }) {
        const hray = new Vec3().setSub(p2, p1);
        hray.normalize(); // normalize() returns length, but modifies in place. We ignore return value.

        // Find best facing chunks
        // Note: In C++, findMostFacingNormal takes ray0 as the *other* point to minimize distance
        const ich1 = this.findMostFacingNormal(hray, chr1, 0.0, true, 1e-6, p2);
        const ich2 = this.findMostFacingNormal(hray, chr2, 0.0, true, 1e-6, p1);

        logger.info(`bridgeFacingPolygons: ich1=${ich1}, ich2=${ich2}`);

        if (ich1 < 0 || ich2 < 0) {
            logger.error(`bridgeFacingPolygons: Could not find facing polygons. ich1=${ich1}, ich2=${ich2}`);
            return;
        }

        // Get vertices for the chunks
        const getChunkVerts = (ich) => {
            const ch = this.chunks[ich];
            const ivs = [];
            for (let i = 0; i < ch.z; i++) {
                ivs.push(this.strips[ch.x + i]);
            }
            // Assuming quads for now as per bridge_quads requirement
            if (ivs.length !== 4) {
                logger.warn(`bridgeFacingPolygons: Chunk ${ich} has ${ivs.length} vertices, expected 4.`);
            }
            return { x: ivs[0], y: ivs[1], z: ivs[2], w: ivs[3] };
        };

        const q1 = getChunkVerts(ich1);
        const q2 = getChunkVerts(ich2);

        this.bridge_quads(q1, q2, nseg, stickTypes, mask, true);
    }

    // --- Face Generation Methods ---

    /**
     * Generate a box face (simple quad)
     * Based on MeshBuilder2.cpp::snapBoxFace
     */
    snapBoxFace(p0, rot, La, Lb) {
        // rot is a Mat3 {a, b, c} where c is normal
        const ha = La * 0.5;
        const hb = Lb * 0.5;

        const v0 = this.vert(new Vec3().setAddMul(p0, rot.a, -ha).addMul(rot.b, -hb));
        const v1 = this.vert(new Vec3().setAddMul(p0, rot.a, ha).addMul(rot.b, -hb));
        const v2 = this.vert(new Vec3().setAddMul(p0, rot.a, ha).addMul(rot.b, hb));
        const v3 = this.vert(new Vec3().setAddMul(p0, rot.a, -ha).addMul(rot.b, hb));

        // Store as chunk
        const stripStart = this.strips.length;
        this.strips.push(v0, v1, v2, v3);
        this.chunk({ x: stripStart, y: 0, z: 4, w: 1 }); // type 1 = face

        return { x: v0, y: v1, z: v2, w: v3 };
    }

    /**
     * Generate a prism face (for 2-edge fork)
     * Based on MeshBuilder2.cpp::snapPrismFace
     */
    snapPrismFace(p0, rot, La, Lb, h, Lbh) {
        const ha = La * 0.5;
        const v0 = this.vert(new Vec3().setAddMul(p0, rot.a, ha).addMul(rot.b, -Lb));
        const v1 = this.vert(new Vec3().setAddMul(p0, rot.c, h).addMul(rot.b, -(Lb - Lbh)));
        const v2 = this.vert(new Vec3().setAddMul(p0, rot.a, -ha).addMul(rot.b, -Lb));
        const stripStart = this.strips.length;
        this.strips.push(v0, v1, v2);
        this.chunk({ x: stripStart, y: 0, z: 3, w: 1 });
        return { x: v0, y: v1, z: v2, w: v2 }; // Degenerate quad
    }

    /**
     * Generate a frustum face (for 3+ edge fork)
     * Based on MeshBuilder2.cpp::snapFrustrumFace
     */
    snapFrustrumFace(p0, rot, La, Lb, h, Lbh, Lah) {
        const v0 = this.vert(new Vec3().setAddMul(p0, rot.a, La).addMul(rot.b, -Lb));
        const v1 = this.vert(new Vec3().setAddMul(p0, rot.a, (La - Lah)).addMul(rot.b, -(Lb - Lbh)).addMul(rot.c, h));
        const v2 = this.vert(new Vec3().setAddMul(p0, rot.a, -(La - Lah)).addMul(rot.b, -(Lb - Lbh)).addMul(rot.c, h));
        const v3 = this.vert(new Vec3().setAddMul(p0, rot.a, -La).addMul(rot.b, -Lb));
        const stripStart = this.strips.length;
        this.strips.push(v0, v1, v2, v3);
        this.chunk({ x: stripStart, y: 0, z: 4, w: 1 });
        return { x: v0, y: v1, z: v2, w: v3 };
    }

    // --- Simple Generators ---

    line(p0, p1, type = -1) {
        const ia = this.vert(p0);
        const ib = this.vert(p1);
        this.edge(ia, ib, type);
    }

    addCube(pos, size, bFaces = true) {
        const s = size * 0.5;
        const x = pos.x !== undefined ? pos.x : pos[0];
        const y = pos.y !== undefined ? pos.y : pos[1];
        const z = pos.z !== undefined ? pos.z : pos[2];
        // 8 vertices
        const v0 = this.vert(new Vec3(x - s, y - s, z - s));
        const v1 = this.vert(new Vec3(x + s, y - s, z - s));
        const v2 = this.vert(new Vec3(x + s, y + s, z - s));
        const v3 = this.vert(new Vec3(x - s, y + s, z - s));
        const v4 = this.vert(new Vec3(x - s, y - s, z + s));
        const v5 = this.vert(new Vec3(x + s, y - s, z + s));
        const v6 = this.vert(new Vec3(x + s, y + s, z + s));
        const v7 = this.vert(new Vec3(x - s, y + s, z + s));
        // 12 edges
        this.edge(v0, v1); this.edge(v1, v2); this.edge(v2, v3); this.edge(v3, v0);
        this.edge(v4, v5); this.edge(v5, v6); this.edge(v6, v7); this.edge(v7, v4);
        this.edge(v0, v4); this.edge(v1, v5); this.edge(v2, v6); this.edge(v3, v7);
        if (bFaces) {
            const ichStart = this.chunks.length;
            this.chunk({ x: this.strips.length, y: 0, z: 4, w: 1 }); this.strips.push(v4, v5, v6, v7);  // Front (z+)
            this.chunk({ x: this.strips.length, y: 0, z: 4, w: 1 }); this.strips.push(v1, v0, v3, v2); // Back (z-)
            this.chunk({ x: this.strips.length, y: 0, z: 4, w: 1 }); this.strips.push(v5, v1, v2, v6); // Right (x+)
            this.chunk({ x: this.strips.length, y: 0, z: 4, w: 1 }); this.strips.push(v0, v4, v7, v3); // Left (x-)
            this.chunk({ x: this.strips.length, y: 0, z: 4, w: 1 }); this.strips.push(v3, v7, v6, v2); // Top (y+)
            this.chunk({ x: this.strips.length, y: 0, z: 4, w: 1 }); this.strips.push(v0, v1, v5, v4); // Bottom (y-)
            return { x: ichStart, y: this.chunks.length };
        }
        return { x: -1, y: -1 };
    }

    /**
     * Creates a girder in the truss structure.
     * Ported from MeshBuilder2.cpp::girder1
     */
    girder1(p0, p1, up, n, width, stickTypes = { x: 1, y: 1, z: 1, w: 1 }, bCaps = false) {
        const dir = new Vec3().setSub(p1, p0);
        const length = dir.normalize(); // Returns length, normalizes in place

        const side = new Vec3();

        // up.makeOrthoU(dir); // C++: up.makeOrthoU(dir);
        // In JS Vec3.js: makeOrthoU(a) -> this -= a * (this . a)
        // But we need to make sure 'up' is not modified in place if it's passed by reference?
        // In JS objects are passed by reference. So we should clone 'up' if we don't want to modify it.
        // However, the C++ code modifies 'up' in place (passed by value? No, Vec3d is usually by value in C++ unless &).
        // In C++ signature: Vec3d up (by value). So we should clone it.

        const upVec = up.clone();
        upVec.makeOrthoU(dir);
        upVec.normalize(); // Ensure it's unit length after orthogonalization

        side.setCross(dir, upVec);
        side.normalize();

        const dl = length / (2 * n + 1);
        const dnp = 4;
        const i00start = this.verts.length;
        let i00 = i00start;

        // Pre-calculate offsets
        // side*-width
        const wSideNeg = new Vec3().setMul(side, new Vec3(-width, -width, -width));
        const wSidePos = new Vec3().setMul(side, new Vec3(width, width, width));
        const wUpNeg = new Vec3().setMul(upVec, new Vec3(-width, -width, -width));
        const wUpPos = new Vec3().setMul(upVec, new Vec3(width, width, width));

        for (let i = 0; i < n; i++) {
            const i01 = i00 + 1;
            const i10 = i00 + 2;
            const i11 = i00 + 3;

            // Vertices
            // vert( p0 + side*-width + dir*(dl*(1+2*i  )) );
            let d = dl * (1 + 2 * i);
            let pBase = new Vec3().setAddMul(p0, dir, d);
            this.vert(new Vec3().setAdd(pBase, wSideNeg));

            // vert( p0 + side*+width + dir*(dl*(1+2*i  )) );
            this.vert(new Vec3().setAdd(pBase, wSidePos));

            // vert( p0 + up  *-width + dir*(dl*(1+2*i+1)) );
            d = dl * (1 + 2 * i + 1);
            pBase.setAddMul(p0, dir, d);
            this.vert(new Vec3().setAdd(pBase, wUpNeg));

            // vert( p0 + up  *+width + dir*(dl*(1+2*i+1)) );
            this.vert(new Vec3().setAdd(pBase, wUpPos));

            // Edges
            this.edge(i00, i01, stickTypes.y);
            this.edge(i10, i11, stickTypes.y);
            this.edge(i00, i10, stickTypes.z);
            this.edge(i00, i11, stickTypes.z);
            this.edge(i01, i10, stickTypes.z);
            this.edge(i01, i11, stickTypes.z);

            if (i < (n - 1)) {
                this.edge(i10, i00 + dnp, stickTypes.w);
                this.edge(i10, i01 + dnp, stickTypes.w);
                this.edge(i11, i00 + dnp, stickTypes.w);
                this.edge(i11, i01 + dnp, stickTypes.w);

                this.edge(i00, i00 + dnp, stickTypes.x);
                this.edge(i01, i01 + dnp, stickTypes.x);
                this.edge(i10, i10 + dnp, stickTypes.x);
                this.edge(i11, i11 + dnp, stickTypes.x);
            }
            i00 += dnp;
        }

        if (bCaps) {
            const ip0 = this.vert(p0);
            const ip1 = this.vert(p1);

            this.edge(i00start + 0, ip0, stickTypes.x);
            this.edge(i00start + 1, ip0, stickTypes.x);
            this.edge(i00start + 2, ip0, stickTypes.x);
            this.edge(i00start + 3, ip0, stickTypes.x);

            const i00end = i00 - dnp;
            this.edge(i00end + 0, ip1, stickTypes.x);
            this.edge(i00end + 1, ip1, stickTypes.x);
            this.edge(i00end + 2, ip1, stickTypes.x);
            this.edge(i00end + 3, ip1, stickTypes.x);
        }

        return i00;
    }

    // --- Edge Adjacency ---

    /**
     * Build edge-to-vertex adjacency map
     * This creates a mapping from each vertex to the edges connected to it
     */
    build_edgesOfVerts() {
        const nv = this.verts.length;
        const ne = this.edges.length;

        // Initialize adjacency map
        this.edgesOfVerts = new Array(nv);
        for (let i = 0; i < nv; i++) {
            this.edgesOfVerts[i] = [];
        }

        // Build adjacency list
        for (let ie = 0; ie < ne; ie++) {
            const e = this.edges[ie];
            this.edgesOfVerts[e.x].push(ie);
            this.edgesOfVerts[e.y].push(ie);
        }

        logger.debug(`build_edgesOfVerts: ${nv} vertices, ${ne} edges`);
    }

    /**
     * Get the other vertex of an edge
     * @param {number} ie - Edge index
     * @param {number} iv - Known vertex index  
     * @returns {number} The other vertex index
     */
    getOtherEdgeVert(ie, iv) {
        const e = this.edges[ie];
        return (e.x === iv) ? e.y : e.x;
    }

    /**
     * Load neighbor vertices and edges for a vertex
     * @param {number} iv - Vertex index
     * @param {Array} ivs - Output array for neighbor vertex indices (can be null)
     * @param {Array} ies - Output array for edge indices (can be null)
     * @param {number} n - Max number to load (-1 for all)
     * @returns {number} Number of neighbors loaded
     */
    loadNeighbours(iv, ivs, ies, n = -1) {
        if (!this.edgesOfVerts) {
            logger.error("loadNeighbours: edgesOfVerts not built. Call build_edgesOfVerts() first.");
            return 0;
        }

        const edgeList = this.edgesOfVerts[iv];
        const total = edgeList.length;
        if (n === -1 || n > total) n = total;

        for (let i = 0; i < n; i++) {
            const ie = edgeList[i];
            if (ies) ies[i] = ie;
            if (ivs) ivs[i] = this.getOtherEdgeVert(ie, iv);
        }

        return n;
    }

    /**
     * Sort vertex edges by angle around a normal vector
     * @param {Vec3} p - Center point
     * @param {Vec3} nor - Normal vector
     * @param {number} n - Number of edges
     * @param {Array} ies - Edge indices (will be sorted in place)
     */
    sortVertEdgesByNormal(p, nor, n, ies) {
        if (n <= 1) return;

        const angles = new Array(n);
        let u = null, v = null;

        // Calculate angles for each edge
        for (let i = 0; i < n; i++) {
            const ie = ies[i];
            const e = this.edges[ie];
            const iv = (e.x !== e.x) ? e.y : ((e.x === e.x) ? this.getOtherEdgeVert(ie, e.x) : e.y);
            const di = new Vec3().setSub(this.verts[iv].pos, p);

            if (i === 0) {
                // First edge defines the basis
                di.makeOrtho(nor);
                u = di.clone();
                u.normalize();
                v = new Vec3().setCross(nor, u);
                v.normalize();
                angles[i] = 0;
            } else {
                // Calculate angle using atan2
                const x = di.dot(u);
                const y = di.dot(v);
                angles[i] = Math.atan2(y, x);
            }
        }

        // Sort indices by angles
        const indices = Array.from({ length: n }, (_, i) => i);
        indices.sort((a, b) => angles[a] - angles[b]);

        // Reorder ies array
        const tempIes = [...ies];
        for (let i = 0; i < n; i++) {
            ies[i] = tempIes[indices[i]];
        }
    }

    /**
     * Bevel a single vertex
     * @param {number} iv - Vertex index to bevel
     * @param {number} L - Bevel distance (polygon radius)
     * @param {number} h - Height offset along normal
     * @param {boolean} bPoly - Whether to create polygon face
     * @param {boolean} bEdgeWedge - Whether to connect original vertex to new vertices
     * @param {Vec3} nor - Normal vector (optional, uses vertex normal if not provided)
     * @returns {number} Number of edges around beveled vertex
     */
    bevel_vert(iv, L, h, bPoly = false, bEdgeWedge = false, nor = null) {
        // Bounds check
        if (iv < 0 || iv >= this.verts.length) {
            logger.error(`bevel_vert: vertex index ${iv} out of bounds [0, ${this.verts.length})`);
            return -1;
        }

        if (!this.edgesOfVerts) {
            logger.error("bevel_vert: edgesOfVerts not built. Call build_edgesOfVerts() first.");
            return -1;
        }

        // Get edges connected to this vertex
        const edgeList = this.edgesOfVerts[iv];
        const ne = edgeList.length;

        logger.debug(`bevel_vert: iv=${iv}, ne=${ne}, L=${L}, h=${h}`);

        if (ne < 1) return 0;

        const ivs = new Array(ne);
        const ies = new Array(ne);
        this.loadNeighbours(iv, ivs, ies, ne);

        const p = this.verts[iv].pos.clone(); // Copy to avoid invalidation

        // Use provided normal or vertex normal
        if (!nor || nor.norm2() < 1e-9) {
            nor = this.verts[iv].nor.clone();
        } else {
            nor = nor.clone();
        }

        // Special case: fewer than 3 edges - create edge instead of polygon
        if (ne < 3) {
            let u = new Vec3();
            if (ne === 1) {
                const d1 = new Vec3().setSub(this.verts[ivs[0]].pos, p);
                d1.normalize();
                u.setCross(d1, nor);
            } else { // ne === 2
                const d1 = new Vec3().setSub(this.verts[ivs[0]].pos, p);
                d1.normalize();
                const d2 = new Vec3().setSub(this.verts[ivs[1]].pos, p);
                d2.normalize();
                const d = new Vec3().setSub(d1, d2);
                u.setCross(d, nor);
                u.normalize();
            }

            const iv1 = this.vert(new Vec3().setAdd(p, new Vec3().setMul(u, new Vec3(L, L, L))).addMul(nor, h));
            const iv2 = this.vert(new Vec3().setAdd(p, new Vec3().setMul(u, new Vec3(-L, -L, -L))).addMul(nor, h));
            this.edge(iv1, iv2);
            if (bEdgeWedge) {
                this.edge(iv, iv1);
                this.edge(iv, iv2);
            }
            return ne;
        }

        // Sort edges by angle around normal
        this.sortVertEdgesByNormal(p, nor, ne, ies);
        // Re-sort vertex indices to match
        for (let i = 0; i < ne; i++) {
            ivs[i] = this.getOtherEdgeVert(ies[i], iv);
        }

        const centralPoint = new Vec3().setAddMul(p, nor, h);
        const newVerts = new Array(ne);

        // Create vertices of the n-gon
        for (let i = 0; i < ne; i++) {
            // Get direction from central vertex to neighbor
            const dir = new Vec3().setSub(this.verts[ivs[i]].pos, p);
            dir.normalize();
            dir.makeOrtho(nor);
            dir.normalize();

            // Previous and next directions (for averaging to get corner directions)
            const prevDir = new Vec3().setSub(this.verts[ivs[(i - 1 + ne) % ne]].pos, p);
            prevDir.makeOrtho(nor);
            prevDir.normalize();

            // Average direction for corner
            const cornerDir = new Vec3().setAdd(prevDir, dir);
            cornerDir.normalize();

            // Position of new vertex
            const newPos = new Vec3().setAddMul(centralPoint, cornerDir, L);
            newVerts[i] = this.vert(newPos);
        }

        // Connect the new vertices to form the n-gon
        const newEdges = new Array(ne);
        for (let i = 0; i < ne; i++) {
            const next = (i + 1) % ne;
            newEdges[i] = this.edge(newVerts[i], newVerts[next]);
        }

        // Create polygon face if requested
        if (bPoly) {
            this.polygon(ne, newVerts);
        }

        // Connect original vertex to new vertices if requested
        if (bEdgeWedge) {
            for (let i = 0; i < ne; i++) {
                this.edge(iv, newVerts[i]);
            }
        }

        return ne;
    }

    /**
     * Create a polygon chunk from vertex indices
     * @param {number} n - Number of vertices
     * @param {Array} ivs - Vertex indices
     */
    polygon(n, ivs) {
        const stripStart = this.strips.length;
        for (let i = 0; i < n; i++) {
            this.strips.push(ivs[i]);
        }
        this.chunk({ x: stripStart, y: 0, z: n, w: 1 }); // type 1 = polygon face
    }

    // --- Export Functions ---

    getVertsFloat32() {
        const arr = new Float32Array(this.verts.length * 3);
        for (let i = 0; i < this.verts.length; i++) {
            arr[i * 3] = this.verts[i].pos.x;
            arr[i * 3 + 1] = this.verts[i].pos.y;
            arr[i * 3 + 2] = this.verts[i].pos.z;
        }
        return arr;
    }

    getEdgesInt32() {
        const arr = new Int32Array(this.edges.length * 2);
        for (let i = 0; i < this.edges.length; i++) {
            arr[i * 2] = this.edges[i].x;
            arr[i * 2 + 1] = this.edges[i].y;
        }
        return arr;
    }

    getChunkStrip(ichunk) {
        const ch = this.chunks[ichunk];
        return this.strips.slice(ch.x, ch.x + ch.z);
    }

    // --- Logging ---

    logOperation(name, args) {
        // Simple logging for now, can be expanded to JSON or other formats
        const argsStr = JSON.stringify(args);
        logger.info(`OP: ${name} ${argsStr}`);
    }

    // --- Validation ---

    /**
     * Validate mesh for NaN, Inf, and coordinate ranges
     * Returns object with validation results and statistics
     */
    validateMesh() {
        const stats = {
            valid: true,
            vertCount: this.verts.length,
            edgeCount: this.edges.length,
            nanCount: 0,
            infCount: 0,
            min: { x: Infinity, y: Infinity, z: Infinity },
            max: { x: -Infinity, y: -Infinity, z: -Infinity },
            errors: []
        };

        if (this.verts.length === 0) {
            stats.valid = false;
            stats.errors.push("Mesh has no vertices");
            return stats;
        }

        // Check each vertex for NaN/Inf and calculate bounds
        for (let i = 0; i < this.verts.length; i++) {
            const v = this.verts[i].pos;

            // Check for NaN
            if (isNaN(v.x) || isNaN(v.y) || isNaN(v.z)) {
                stats.nanCount++;
                stats.valid = false;
                if (stats.errors.length < 10) { // Limit error reporting
                    stats.errors.push(`Vertex ${i} has NaN: (${v.x}, ${v.y}, ${v.z})`);
                }
                continue; // Skip bounds calculation for NaN vertices
            }

            // Check for Inf
            if (!isFinite(v.x) || !isFinite(v.y) || !isFinite(v.z)) {
                stats.infCount++;
                stats.valid = false;
                if (stats.errors.length < 10) {
                    stats.errors.push(`Vertex ${i} has Inf: (${v.x}, ${v.y}, ${v.z})`);
                }
                continue;
            }

            // Update bounds
            stats.min.x = Math.min(stats.min.x, v.x);
            stats.min.y = Math.min(stats.min.y, v.y);
            stats.min.z = Math.min(stats.min.z, v.z);
            stats.max.x = Math.max(stats.max.x, v.x);
            stats.max.y = Math.max(stats.max.y, v.y);
            stats.max.z = Math.max(stats.max.z, v.z);
        }

        // Calculate spans
        if (stats.valid || (stats.nanCount + stats.infCount < this.verts.length)) {
            stats.span = {
                x: stats.max.x - stats.min.x,
                y: stats.max.y - stats.min.y,
                z: stats.max.z - stats.min.z
            };

            // Check for zero span (all vertices in a plane or line)
            const epsilon = 1e-10;
            if (Math.abs(stats.span.x) < epsilon && Math.abs(stats.span.y) < epsilon && Math.abs(stats.span.z) < epsilon) {
                stats.errors.push(`All vertices are at the same point: (${stats.min.x}, ${stats.min.y}, ${stats.min.z})`);
                logger.warn(`All vertices collapsed to a point`);
            } else if (Math.abs(stats.span.x) < epsilon || Math.abs(stats.span.y) < epsilon || Math.abs(stats.span.z) < epsilon) {
                logger.warn(`Mesh is degenerate (collapsed in one dimension): span=(${stats.span.x}, ${stats.span.y}, ${stats.span.z})`);
            }
        }

        // Log summary
        if (stats.nanCount > 0) {
            logger.error(`Mesh validation failed: ${stats.nanCount} vertices with NaN`);
        }
        if (stats.infCount > 0) {
            logger.error(`Mesh validation failed: ${stats.infCount} vertices with Inf`);
        }
        if (stats.valid) {
            logger.info(`Mesh validation passed: ${stats.vertCount} verts, bounds: [${stats.min.x.toFixed(2)},${stats.max.x.toFixed(2)}] x [${stats.min.y.toFixed(2)},${stats.max.y.toFixed(2)}] x [${stats.min.z.toFixed(2)},${stats.max.z.toFixed(2)}]`);
        }

        return stats;
    }

    // --- OBJ Export ---

    toObjString() {
        let str = "# SimpleSimulationEngine MeshBuilder JS OBJ export\n";

        // Vertices
        str += `# Vertices: ${this.verts.length}\n`;
        for (const v of this.verts) {
            str += `v ${v.pos.x} ${v.pos.y} ${v.pos.z}\n`;
        }

        // Normals (optional, but good to have)
        str += `# Normals: ${this.verts.length}\n`;
        for (const v of this.verts) {
            str += `vn ${v.nor.x} ${v.nor.y} ${v.nor.z}\n`;
        }

        // UVs
        str += `# UVs: ${this.verts.length}\n`;
        for (const v of this.verts) {
            str += `vt ${v.uv.x} ${v.uv.y}\n`;
        }

        // Faces (from chunks)
        // We iterate over chunks to reconstruct faces
        // Assuming chunks are faces (type 1)
        str += `# Faces\n`;
        for (const ch of this.chunks) {
            // ch.w is type. 1 = face.
            // if (ch.w === 1) { 
            // For now, let's assume all chunks that look like faces are faces.
            // Or just dump all strips as faces if they form valid polygons.

            const ivs = this.strips.slice(ch.x, ch.x + ch.z);
            if (ivs.length >= 3) {
                str += "f";
                for (const idx of ivs) {
                    // OBJ indices are 1-based
                    const i = idx + 1;
                    // f v/vt/vn
                    str += ` ${i}/${i}/${i}`;
                }
                str += "\n";
            }
            // }
        }

        // Edges (as lines)
        str += `# Edges: ${this.edges.length}\n`;
        for (const e of this.edges) {
            // l v1 v2
            str += `l ${e.x + 1} ${e.y + 1}\n`;
        }

        return str;
    }
}

// Export for both module and non-module usage
if (typeof module !== 'undefined' && module.exports) {
    // Node.js
    module.exports = { MeshBuilder };
} else if (typeof window !== 'undefined') {
    window.MeshBuilder = MeshBuilder;
}
