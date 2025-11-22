
class MeshBuilder {
    constructor() {
        // Vertices: Array of {pos: Vec3, nor: Vec3, uv: {x, y}}
        this.verts = [];
        // Edges: Array of {x: int, y: int, z: type, w: type2}
        this.edges = [];
        // Chunks: Array of {x: stripStart, y: edgeStart, z: count, w: type}
        this.chunks = [];
        // Strips: Flat array of vertex/edge indices for faces
        this.strips = [];
        // Blocks: [ivert_start, iedge_start, ichunk_start]
        this.blocks = [];

        // Temporary vectors (reused to avoid GC)
        this._tmp1 = new Vec3();
        this._tmp2 = new Vec3();
        this._tmp3 = new Vec3();
    }

    clear() {
        this.verts = [];
        this.edges = [];
        this.chunks = [];
        this.strips = [];
        this.blocks = [];
        if (window.VERBOSITY_LEVEL > 0) window.logger.info("MeshBuilder cleared.");
    }

    // --- Basic Primitives ---

    vert(pos) {
        const v = {
            pos: new Vec3(pos.x || pos[0], pos.y || pos[1], pos.z || pos[2]),
            nor: new Vec3(0, 0, 1),
            uv: { x: 0, y: 0 }
        };
        this.verts.push(v);
        const idx = this.verts.length - 1;
        if (window.VERBOSITY_LEVEL >= 3) window.logger.debug(`Vert[${idx}]: ${v.pos.x}, ${v.pos.y}, ${v.pos.z}`);
        return idx;
    }

    edge(a, b, type = -1, type2 = 0) {
        this.edges.push({ x: a, y: b, z: type, w: type2 });
        const idx = this.edges.length - 1;
        if (window.VERBOSITY_LEVEL >= 3) window.logger.debug(`Edge[${idx}]: ${a} -> ${b} (t=${type})`);
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

    // --- Alignment Function ---

    /**
     * Align two polygons rotationally to minimize twist
     * Based on MeshBuilder2.cpp::alling_polygons
     */
    alling_polygons(n, ivs1, ivs2, ipiv = 0) {
        const cog1 = this.getCOG(n, ivs1);
        const cog2 = this.getCOG(n, ivs2);

        // Compute axis
        const ax = this._tmp1.setSub(cog2, cog1).normalize();

        // Compute orthogonal basis
        const u = this._tmp2.setSub(this.verts[ivs1[ipiv]].pos, cog1);
        u.makeOrthoU(ax);
        u.normalize();

        const v = this._tmp3.setCross(ax, u).normalize();

        // Project points to UV plane
        const uv1 = new Array(n);
        const uv2 = new Array(n);

        for (let i = 0; i < n; i++) {
            const d1 = new Vec3().setSub(this.verts[ivs1[i]].pos, cog1);
            const d2 = new Vec3().setSub(this.verts[ivs2[i]].pos, cog2);
            uv1[i] = new Vec3(d1.dot(u), d1.dot(v), 0).normalize();
            uv2[i] = new Vec3(d2.dot(u), d2.dot(v), 0).normalize();
        }

        // Match points by nearest UV projection
        // Store original ivs2 values before remapping
        const originalIvs2 = ivs2.slice();

        for (let i = 0; i < n; i++) {
            let dmin = -1.0;
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
        if (window.VERBOSITY_LEVEL >= 2) {
            window.logger.info(`bridge_quads: nseg=${nseg}, q1=[${q1.x},${q1.y},${q1.z},${q1.w}], q2=[${q2.x},${q2.y},${q2.z},${q2.w}]`);
        }

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

    addCube(pos, size) {
        const s = size * 0.5;
        const x = pos.x || pos[0];
        const y = pos.y || pos[1];
        const z = pos.z || pos[2];

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
    }

    girder1(p0, p1, up, n, width, type = -1) {
        const dir = new Vec3().setSub(p1, p0);
        let prevI = this.vert(p0);

        for (let i = 1; i <= n; i++) {
            const t = i / n;
            const pos = new Vec3().setAddMul(p0, dir, t);
            const currI = this.vert(pos);
            this.edge(prevI, currI, type);
            prevI = currI;
        }
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
}
