
class MeshBuilder {
    constructor() {
        // Vertices: [x, y, z] flat array
        this.verts = [];
        // Edges: [a, b, type, type2] flat array
        this.edges = [];

        // Blocks: [ivert_start, iedge_start, itri_start, ichunk_start]
        this.blocks = [];

        // Use global window.VERBOSITY_LEVEL
    }

    clear() {
        this.verts = [];
        this.edges = [];
        this.blocks = [];
        if (window.VERBOSITY_LEVEL > 0) window.logger.info("MeshBuilder cleared.");
    }

    // --- Basic Primitives ---

    vert(pos) {
        this.verts.push(pos[0], pos[1], pos[2]);
        const idx = (this.verts.length / 3) - 1;
        if (window.VERBOSITY_LEVEL >= 3) window.logger.debug(`Vert[${idx}]: ${pos[0]}, ${pos[1]}, ${pos[2]}`);
        return idx;
    }

    edge(a, b, type = -1) {
        this.edges.push(a, b, type, 0);
        const idx = (this.edges.length / 4) - 1;
        if (window.VERBOSITY_LEVEL >= 3) window.logger.debug(`Edge[${idx}]: ${a} -> ${b} (t=${type})`);
        return idx;
    }

    block() {
        const b = {
            ivert: this.verts.length / 3,
            iedge: this.edges.length / 4
        };
        this.blocks.push(b);
        return b;
    }

    // --- Generators (Simplified from C++) ---

    // Generate a single line segment (girder segment)
    line(p0, p1, type = -1) {
        const ia = this.vert(p0);
        const ib = this.vert(p1);
        this.edge(ia, ib, type);
    }

    // --- Advanced Generators ---

    // Add a cube centered at pos with given size
    addCube(pos, size) {
        const s = size * 0.5;
        const x = pos[0], y = pos[1], z = pos[2];

        // 8 vertices
        const v0 = this.vert([x - s, y - s, z - s]);
        const v1 = this.vert([x + s, y - s, z - s]);
        const v2 = this.vert([x + s, y + s, z - s]);
        const v3 = this.vert([x - s, y + s, z - s]);
        const v4 = this.vert([x - s, y - s, z + s]);
        const v5 = this.vert([x + s, y - s, z + s]);
        const v6 = this.vert([x + s, y + s, z + s]);
        const v7 = this.vert([x - s, y + s, z + s]);

        // 12 edges
        this.edge(v0, v1); this.edge(v1, v2); this.edge(v2, v3); this.edge(v3, v0); // Back face
        this.edge(v4, v5); this.edge(v5, v6); this.edge(v6, v7); this.edge(v7, v4); // Front face
        this.edge(v0, v4); this.edge(v1, v5); this.edge(v2, v6); this.edge(v3, v7); // Connecting edges
    }

    // Generate a multi-segment girder
    // p0, p1: vec3 (start, end)
    // up: vec3 (orientation)
    // n: int (segments)
    // width: float
    girder1(p0, p1, up, n, width, type = -1) {
        // For Phase 1/2, we will just generate a single line of edges for the girder spine.
        // Later we can generate the full truss structure (longerons + cross-bracing).

        const dir = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];

        let prevI = this.vert(p0);

        for (let i = 1; i <= n; i++) {
            const t = i / n;
            const pos = [
                p0[0] + dir[0] * t,
                p0[1] + dir[1] * t,
                p0[2] + dir[2] * t
            ];
            const currI = this.vert(pos);
            this.edge(prevI, currI, type);
            prevI = currI;
        }
    }

    // Helper to get typed arrays for upload
    getVertsFloat32() {
        return new Float32Array(this.verts);
    }

    getEdgesInt32() {
        return new Int32Array(this.edges);
    }
}
