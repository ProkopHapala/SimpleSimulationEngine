import { Vec3 } from './Vec3.js';
import { SelectionBanks } from './Selection.js';
import { logger } from './Logger.js';

export class MeshBuilder {
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

        // Topology / de-duplication controls (see SpaceCraftConstructionProblems.md §1.7).
        // By default these checks are disabled to preserve legacy behavior and performance.
        // They can be enabled by generators that care about topological uniqueness.
        this.bCheckVertExist = false;      // If true, attempt to detect duplicate vertices within RvertCollapse.
        this.bVertExistError = true;       // If true and duplicate is found, throw; otherwise reuse/skip.
        this.bVertExistSkip  = true;       // In non-error mode, return existing vertex index instead of adding.
        this.RvertCollapse   = 1e-3;       // Distance tolerance for detecting duplicate vertices.

        this.bCheckEdgeExist = false;      // If true, attempt to detect duplicate edges (undirected).
        this.bEdgeExistError = true;       // If true and duplicate is found, throw; otherwise reuse/skip.
        this.bEdgeExistSkip  = true;       // In non-error mode, return existing edge index instead of adding.

        // Internal map for edge de-duplication. Key is a packed integer based on (i0,i1).
        // NOTE: For now we assume meshes are small enough that i0 * EDGE_KEY_STRIDE + i1
        //       fits safely into JS Number without precision issues.
        this._edgeMap = new Map();   // key:number -> edge index
        this.EDGE_KEY_STRIDE = 65536; // 2^16; max supported vertex index ~2^16

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
        // Also clear any cached de-duplication maps; topology will be rebuilt from scratch.
        if (this._edgeMap) this._edgeMap.clear();
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

    _edgeKey(a, b) {
        const i0 = a < b ? a : b;
        const i1 = a < b ? b : a;
        // Pack into a single integer-like Number. This is safe as long as
        // i0,i1 << EDGE_KEY_STRIDE and overall key < 2^53.
        return i0 * this.EDGE_KEY_STRIDE + i1;
    }

    /**
     * Find index of an existing vertex within distance Rmax of point p0.
     * Returns -1 if none is found.
     * This is a simple brute-force search over all verts, analogous in spirit to Mesh::Builder2::findVert, and isacceptable for current small meshes. 
     * Later we can replace the internals with a spatial hash / BVH / sweep structure without changing callers.
     */
    findVert(p0, Rmax = this.RvertCollapse, iStart = 0) {
        if (!this.verts || this.verts.length === 0) return -1;
        const R2 = Rmax * Rmax;
        let r2min = R2;
        let imin = -1;
        const n = this.verts.length;
        for (let i = iStart; i < n; i++) {
            const p = this.verts[i].pos;
            if (!p) continue;
            const dx = p.x - p0.x;
            const dy = p.y - p0.y;
            const dz = p.z - p0.z;
            const r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < r2min) {
                r2min = r2;
                imin = i;
            }
        }
        return imin;
    }

    scanDuplicateVerts(Rmax = this.RvertCollapse) {
        const n = this.verts.length;
        if (n === 0) { logger.info("scanDuplicateVerts: mesh has no vertices."); return 0; }
        let dupCount = 0;
        for (let i = 0; i < n; i++) {
            const p0 = this.verts[i].pos;
            if (!p0) continue;
            let start = i + 1;
            while (true) {
                const j = this.findVert(p0, Rmax, start);
                if (j < 0) break;
                const p1 = this.verts[j].pos;
                logger.info( `  dupVert i=${i} and j=${j} | v[${i}]{${p0.x},${p0.y},${p0.z}} v[${j}]{${p1.x},${p1.y},${p1.z}}` );
                dupCount++;
                start = j + 1;
            }                     // just skip this i
        }
        logger.info(`scanDuplicateVerts: n=${n}, Rmax=${Rmax}, duplicates=${dupCount}`);
        return dupCount;
    }

    scanDuplicateEdges() {
        const n = this.edges.length;
        if (n === 0) { logger.info("scanDuplicateEdges: mesh has no edges."); return 0; }
        // Reuse the builder's edge map as a workspace to detect duplicates.
        if (this._edgeMap) this._edgeMap.clear();
        let dupCount = 0;
        for (let i = 0; i < n; i++) {
            const e   = this.edges[i];
            const key = this._edgeKey(e.x,e.y);
            if (this._edgeMap.has(key)) {
                const first = this._edgeMap.get(key);
                logger.info(`  dupEdge ie=${first} and ie=${i} between verts (${e.x}, ${e.y})`);
                dupCount++;
            } else {
                this._edgeMap.set(key, i);
            }
        }
        logger.info(`scanDuplicateEdges: n=${n}, duplicates=${dupCount}`);
        return dupCount;
    }

    getEdgeLength(e) {
        const va = this.verts[e.x];
        const vb = this.verts[e.y];
        if (!va || !vb || !va.pos || !vb.pos) return -1;
        const pa = va.pos;
        const pb = vb.pos;
        const dx = pa.x - pb.x;
        const dy = pa.y - pb.y;
        const dz = pa.z - pb.z;
        return Math.hypot(dx, dy, dz);
    }

    /**
     * Insert a bond length into a sorted unique-length list with tolerance.
     *
     * @param {number} L         Measured bond length.
     * @param {*}      key       Identifier for the bond (edge index, pair, etc.).
     * @param {number[]} uniqLs  Sorted array of representative lengths (will be mutated).
     * @param {Array[]} groups   Optional parallel array of arrays; groups[i] collects keys
     *                           for uniqLs[i]. If null/undefined, keys are ignored.
     * @param {number} dR        Tolerance for merging lengths (default 1e-6).
     * @returns {number}         Index in uniqLs where L was merged/inserted.
     */
    registerBondLength( L, key, uniqLs, groups = null, dR = 1e-6) {   
        let lo = 0;
        let hi = uniqLs.length;
        while (lo < hi) {
            const mid = (lo + hi) >> 1;
            const v = uniqLs[mid];
            const diff = L - v;
            if (Math.abs(diff) <= dR) {
                if (groups) {
                    if (!groups[mid]) groups[mid] = [];
                    groups[mid].push(key);
                }
                return mid;
            }
            if (diff < 0) { hi = mid; } else { lo = mid + 1; }
        }
        const idx = lo;
        uniqLs.splice(idx, 0, L);
        if (groups) { groups.splice(idx, 0, [key]);}
        return idx;
    }

    /**
     * Build a list of unique bond lengths from explicit edges.
     * @param {Array|undefined|null} edgeList  Optional list of edges/bonds.
     * @param {number}               dR        Length merging tolerance.
     * @param {Array[]|null}         groups    Optional output groups parallel to uniqLs.
     * @returns {number[]}                      Sorted unique bond lengths.
     */
    mapBondLengthsFromEdges(edgeList, dR = 1e-6, groups = null) {
        const uniqLs = [];
        for (const item of edgeList) {
            const e = this.edges[item];
            if (!e) continue;
            const L = this.getEdgeLength(e);
            if (L < 0) continue;
            this.registerBondLength(L, item, uniqLs, groups, dR);
        }
        return uniqLs;
    }

    /**
     * Build a list of unique bond lengths from two vertex selections.
     *
     * This does **not** create edges; it just treats all candidate pairs
     * (ia in selA, ib in selB) whose distance is below Rmax as potential
     * bonds and feeds them through insertBondLength.
     *
     * Keys stored into groups (if provided) are simple 2-element arrays
     * [ia, ib].
     *
     * @param {number[]} selA   List of vertex indices (first set).
     * @param {number[]} selB   List of vertex indices (second set).
     * @param {number}   Rmax   Cutoff distance for candidates (Infinity = no cutoff).
     * @param {number}   dR     Length merging tolerance.
     * @param {Array[]|null} groups Optional output groups parallel to uniqLs.
     * @returns {number[]}        Sorted unique bond lengths.
     */
    mapBondLengthsFromVertexSelections(selA, selB, Rmax = Infinity, dR = 1e-6, groups = null) {
        const uniqLs = [];
        if (!Array.isArray(selA) || !Array.isArray(selB) || selA.length === 0 || selB.length === 0) { return uniqLs;}
        const useCutoff = Number.isFinite(Rmax) && Rmax > 0;
        const R2 = useCutoff ? Rmax * Rmax : Infinity;
        for (const ia of selA) {
            const va = this.verts[ia];
            if (!va || !va.pos) continue;
            const pa = va.pos;
            for (const ib of selB) {
                const vb = this.verts[ib];
                if (!vb || !vb.pos) continue;
                const l2 = pa.dist2(vb.pos);
                if (useCutoff && l2 > R2) continue;
                const key = [ia, ib];
                this.registerBondLength(Math.sqrt(l2), key, uniqLs, groups, dR);
            }
        }

        return uniqLs;
    }

    vert(pos) {
        const x = pos.x !== undefined ? pos.x : pos[0];
        const y = pos.y !== undefined ? pos.y : pos[1];
        const z = pos.z !== undefined ? pos.z : pos[2];

        if (this.bCheckVertExist) {
            const iv = this.findVert(new Vec3(x, y, z), this.RvertCollapse);
            if (iv >= 0) {
                if (this.bVertExistError) { throw new Error(`MeshBuilder.vert: vertex already exists within RvertCollapse at iv=${iv}`); }
                if (this.bVertExistSkip ) { return iv; }
                // Otherwise fall through and create a near-duplicate explicitly.
            }
        }
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
        const idx = this.edges.length;
        if (this.bCheckEdgeExist) {
            const key = this._edgeKey(a, b);
            const existing = this._edgeMap.get(key);
            if (existing !== undefined) {
                if (this.bEdgeExistError) { throw new Error(`MeshBuilder.edge: edge already exists between ${a} and ${b} -> ie=${existing}`);  }
                if (this.bEdgeExistSkip ) { return existing; }
                // If neither error nor skip is desired, fall through and create a duplicate explicitly.
            }
            this._edgeMap.set(key, idx); // Cache new  edge for de-duplication
        }
        this.edges.push({ x: a, y: b, z: type, w: type2 });
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

    /**
     * Add a sequence of edges from a list of index pairs.
     * Each element of pairs is expected to be [i0, i1].
     * This is a generic helper so GUI/tests do not need to inline edge loops.
     */
    addEdgesFromPairs(pairs, type = -1, type2 = 0) {
        if (!pairs) return;
        for (const pair of pairs) {
            if (!pair || pair.length < 2) continue;
            const a = pair[0];
            const b = pair[1];
            this.edge(a, b, type, type2);
        }
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

    /**
     * Generate a polyline ("rope") with nSeg points starting from pStart
     * and going along the given direction for the given length.
     * Returns an array of newly created vertex indices in order.
     */
    rope(pStart, dir, length, nSeg, type = -1) {
        nSeg = nSeg | 0;
        if (nSeg <= 0 || length === 0) return [];

        const d = new Vec3(dir.x, dir.y, dir.z);
        const L = d.norm();
        if (L <= 0) return [];
        d.mulScalar(length / L);

        const p0 = new Vec3(pStart.x, pStart.y, pStart.z);
        const verts = [];
        if (nSeg === 1) {
            const iv = this.vert(p0);
            verts.push(iv);
            return verts;
        }

        const step = 1.0 / (nSeg - 1);
        let prev = -1;
        for (let i = 0; i < nSeg; i++) {
            const t = i * step;
            const p = new Vec3().setAddMul(p0, d, t);
            const iv = this.vert(p);
            verts.push(iv);
            if (prev >= 0) {
                this.edge(prev, iv, type);
            }
            prev = iv;
        }
        return verts;
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

    // --- Plates Between Ropes / Girders (JS P2 helpers) ---

    /**
     * Extract an ordered strip of vertex indices along a geometric line segment.
     * Vertices within distance r from the infinite line p0-p1 are collected and
     * sorted by projection onto the line direction.
     */
    extractStripAlongLine(p0, p1, r, maxCount = -1) {
        const dir = new Vec3().setSub(p1, p0);
        const len = dir.norm();
        if (len <= 0) return [];
        dir.mulScalar(1.0 / len);

        const r2 = r * r;
        const candidates = [];
        for (let i = 0; i < this.verts.length; i++) {
            const pos = this.verts[i].pos;
            const d = new Vec3().setSub(pos, p0);
            const t = d.dot(dir);
            const proj = new Vec3().setAddMul(p0, dir, t);
            const off = new Vec3().setSub(pos, proj);
            const dist2 = off.dot(off);
            if (dist2 <= r2) {
                candidates.push({ i, t });
            }
        }

        candidates.sort((a, b) => a.t - b.t);
        if (maxCount > 0 && candidates.length > maxCount) {
            // Keep roughly evenly spaced subset
            const out = [];
            const step = (candidates.length - 1) / (maxCount - 1);
            for (let k = 0; k < maxCount; k++) {
                const idx = Math.round(k * step);
                out.push(candidates[idx].i);
            }
            return out;
        }

        return candidates.map(c => c.i);
    }

    /**
     * Connect two ordered vertex strips with cross edges.
     * For now this only builds a ladder of edges between corresponding
     * vertices; it can be extended later to generate full plates.
     */
    plateBetweenVertStrips(strip1, strip2) {
        const n = Math.min(strip1.length, strip2.length);
        if (n < 2) return 0;
        let cnt = 0;
        for (let i = 0; i < n; i++) {
            this.edge(strip1[i], strip2[i]);
            cnt++;
        }
        return cnt;
    }

    /**
     * Triangulate the strip between two ordered vertex strips using the
     * parametric algorithm from doc/python/interpolated_vertex_count.py
     * (mode='parametric' only).
     *
     * We assume that:
     * - strip1 and strip2 are ordered along their respective polylines.
     * - Axial rope edges (along each strip) are already present.
     * - Optional cross edges (from plateBetweenVertStrips) may also exist.
     *
     * This method only adds the **diagonal edges** that choose, for each
     * quad cell, which way to split it, using index-based parametric sync
     * rather than geometric distances.
     */
    triangulateBetweenVertStrips(strip1, strip2) {
        const na = strip1.length;
        const nb = strip2.length;
        if (na < 2 || nb < 2) return 0;

        let ia = 0;
        let ib = 0;
        let added = 0;

        while (ia < na - 1 || ib < nb - 1) {
            let choice = null;

            if (ia === na - 1) {
                // A finished, must advance B
                choice = 'B';
            } else if (ib === nb - 1) {
                // B finished, must advance A
                choice = 'A';
            } else {
                // Both available: advance the side that is "behind" in param space
                const uNextA = (ia + 1) / (na - 1);
                const uNextB = (ib + 1) / (nb - 1);
                choice = (uNextA <= uNextB) ? 'A' : 'B';
            }

            if (choice === 'A') {
                // Triangle (A_curr, B_curr, A_next) -> add diagonal (A_next, B_curr)
                const iaNext = strip1[ia + 1];
                const ibCurr = strip2[ib];
                this.edge(iaNext, ibCurr);
                added++;
                ia++;
            } else {
                // Triangle (A_curr, B_curr, B_next) -> add diagonal (A_curr, B_next)
                const iaCurr = strip1[ia];
                const ibNext = strip2[ib + 1];
                this.edge(iaCurr, ibNext);
                added++;
                ib++;
            }
        }

        return added;
    }

    /**
     * High-level helper for P2: build cross connections between two
     * polylines defined by 3 or 4 corner vertices.
     *
     * - 3 vertices: treat as V/L case (shared corner + two tips).
     * - 4 vertices: treat as quad case (two opposite edges).
     */
    plateBetweenEdges(corners, r, maxPerStrip = -1) {
        if (!Array.isArray(corners)) return 0;
        const n = corners.length;
        if (n !== 3 && n !== 4) {
            logger.error(`plateBetweenEdges: expected 3 or 4 corners, got ${n}`);
            return 0;
        }

        const getPos = (iv) => this.verts[iv]?.pos;

        let strip1 = [];
        let strip2 = [];

        if (n === 3) {
            const iCenter = corners[1];
            const iA = corners[0];
            const iB = corners[2];
            const pC = getPos(iCenter);
            const pA = getPos(iA);
            const pB = getPos(iB);
            if (!pC || !pA || !pB) return 0;

            strip1 = this.extractStripAlongLine(pC, pA, r, maxPerStrip);
            strip2 = this.extractStripAlongLine(pC, pB, r, maxPerStrip);

            // In the triangle case we skip the shared corner so we get
            // two strips that can be treated like a quad strip.
            if (strip1.length > 0 && strip2.length > 0) {
                strip1 = strip1.slice(1);
                strip2 = strip2.slice(1);
            }
        } else if (n === 4) {
            const i0 = corners[0];
            const i1 = corners[1];
            const i2 = corners[2];
            const i3 = corners[3];
            const p0 = getPos(i0);
            const p1 = getPos(i1);
            const p2 = getPos(i2);
            const p3 = getPos(i3);
            if (!p0 || !p1 || !p2 || !p3) return 0;

            strip1 = this.extractStripAlongLine(p0, p1, r, maxPerStrip);
            strip2 = this.extractStripAlongLine(p3, p2, r, maxPerStrip);
        }

        if (!strip1.length || !strip2.length) {
            logger.warn('plateBetweenEdges: one or both strips are empty');
            return 0;
        }

        const used = this.plateBetweenVertStrips(strip1, strip2);
        logger.info(`plateBetweenEdges: connected ${used} pairs between strips (len1=${strip1.length}, len2=${strip2.length})`);
        return used;
    }

    /**
     * Triangulated variant of plateBetweenEdges: builds strips from corners
     * as in plateBetweenEdges, then adds diagonal edges following
     * triangulateBetweenVertStrips.
     */
    triPlateBetweenEdges(corners, r, maxPerStrip = -1) {
        if (!Array.isArray(corners)) return 0;
        const n = corners.length;
        if (n !== 3 && n !== 4) {
            logger.error(`triPlateBetweenEdges: expected 3 or 4 corners, got ${n}`);
            return 0;
        }

        const getPos = (iv) => this.verts[iv]?.pos;

        let strip1 = [];
        let strip2 = [];

        if (n === 3) {
            const iCenter = corners[1];
            const iA = corners[0];
            const iB = corners[2];
            const pC = getPos(iCenter);
            const pA = getPos(iA);
            const pB = getPos(iB);
            if (!pC || !pA || !pB) return 0;

            strip1 = this.extractStripAlongLine(pC, pA, r, maxPerStrip);
            strip2 = this.extractStripAlongLine(pC, pB, r, maxPerStrip);

            if (strip1.length > 0 && strip2.length > 0) {
                strip1 = strip1.slice(1);
                strip2 = strip2.slice(1);
            }
        } else if (n === 4) {
            const i0 = corners[0];
            const i1 = corners[1];
            const i2 = corners[2];
            const i3 = corners[3];
            const p0 = getPos(i0);
            const p1 = getPos(i1);
            const p2 = getPos(i2);
            const p3 = getPos(i3);
            if (!p0 || !p1 || !p2 || !p3) return 0;

            strip1 = this.extractStripAlongLine(p0, p1, r, maxPerStrip);
            strip2 = this.extractStripAlongLine(p3, p2, r, maxPerStrip);
        }

        if (!strip1.length || !strip2.length) {
            logger.warn('triPlateBetweenEdges: one or both strips are empty');
            return 0;
        }

        // First ensure we have a basic ladder of cross edges
        this.plateBetweenVertStrips(strip1, strip2);

        const added = this.triangulateBetweenVertStrips(strip1, strip2);
        logger.info(`triPlateBetweenEdges: added ${added} diagonals between strips (len1=${strip1.length}, len2=${strip2.length})`);
        return added;
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

// Also expose on window for legacy global-script users
if (typeof window !== 'undefined') {
    window.MeshBuilder = MeshBuilder;
}
