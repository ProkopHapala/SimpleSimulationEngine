// Abstract Component Classes

export class Node {
    constructor(pos, size = 1, type = 0) {
        this.pos = pos; // [x, y, z]
        this.size = size;
        this.type = type;
        this.id = -1; // Assigned by Engine
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
        this.chunkRange = { x: -1, y: -1 };
    }
}

export class Girder {
    constructor(nodeA, nodeB, nseg = 1, type = 0) {
        this.nodeA = nodeA; // Node object reference
        this.nodeB = nodeB; // Node object reference
        this.type = type;
        this.nseg = nseg; // segments along girder
        this.up = [0, 1, 0];
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }

    getVertByStride(i, j) {
        if (this.pointRange.x < 0) return -1;
        const nVertPerSeg = 4;
        return this.pointRange.x + i * nVertPerSeg + j;
    }

    sideToPath(side) {
        if (this.pointRange.x < 0) return [];
        const nVertPerSeg = 4;
        const nverts = this.pointRange.y - this.pointRange.x;
        const nrings = Math.floor(nverts / nVertPerSeg);
        const ps = [];
        const offset = side % nVertPerSeg;
        // For girders (bridged quads), we have 'nrings' cross-sections.
        // Each cross-section has 4 vertices. 
        for (let i = 0; i < nrings; i++) {
            ps.push(this.pointRange.x + i * nVertPerSeg + offset);
        }
        console.log(`[SpaceCraft] Girder ${this.id} sideToPath(${side}): nrings=${nrings} ps=[${ps[0]}...${ps[ps.length-1]}]`);
        return ps;
    }
}

export class Rope {
    constructor(nodeA, nodeB, thick, type = 0) {
        this.nodeA = nodeA;
        this.nodeB = nodeB;
        this.thick = thick;
        this.type = type;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }
}

export class Plate {
    constructor(boundA, boundB, spanA = [0, 1], spanB = [0, 1], type = 0, kind = 'Radiator', nx = 2, ny = 2, nz = 1, upA = true, upB = true, sideOffset = 0, weldDist = 0) {
        this.boundA = boundA; // Girder or Rope
        this.boundB = boundB; // Girder or Rope
        this.spanA = spanA;   // [cmin,cmax] along boundA
        this.spanB = spanB;   // [cmin,cmax] along boundB
        this.kind = kind;     // 'Radiator' or 'Shield'
        this.nx = nx; this.ny = ny; this.nz = nz; // segment counts for ParametricQuadPatch
        this.upA = upA; this.upB = upB; // which edge (upper/lower) to attach on each girder
        this.sideOffset = sideOffset; // shift along side axis to snap to girder corner
        this.weldDist = weldDist;     // weld radius to connect to nearest girder verts
        this.type = type;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }
}

export class Slider {
    constructor(rail, sliding = null, calong = 0.0, matName = null, slidingVertId = -1, side = 0, methodFlag = true) {
        // rail: structural component providing the path
        // sliding: structural component that owns the sliding vertex
        this.rail = rail;
        this.sliding = sliding;
        this.slidingVertId = slidingVertId; // vertex index on sliding component (within mesh)
        this.calong = calong;   // param along path (0..1), wraps if path is closed
        this.side = side;       // desired rail side (stride mode)
        this.methodFlag = methodFlag; // true=stride, false=SDF
        this.type = matName;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
        this.path = { ps: [], cur: 0, closed: rail instanceof Ring };
        this.ivert = -1;        // actual slider vertex id in mesh (created during build)
        this.weldDist = 0.25;   // default weld radius to connect slider vertex to rail segment
        this.pathRadius = 0.35; // SDF cylinder radius for path picking
    }
}

export class Ring {
    constructor(pos, dir, up, R, nseg, wh, matName, st) {
        this.pos = pos; // Vec3
        this.dir = dir; // Vec3 (axis)
        this.up = up;   // Vec3 (stable up)
        this.R = R;
        this.nseg = nseg;
        this.wh = wh;
        this.matName = matName;
        this.st = st;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }

    static from3Points(p1, p2, p3, nseg, wh, matName, st) {
        const { center, radius, x, y } = Vec3.circle3Point(p1, p2, p3);
        const axis = new Vec3().setCross(x, y);
        axis.normalize();
        return new Ring(center, axis, x, radius, nseg, wh, matName, st);
    }

    getVertByStride(i, j) {
        if (this.pointRange.x < 0) return -1;
        const nVertPerSeg = 4;
        return this.pointRange.x + i * nVertPerSeg + j;
    }

    sideToPath(side) {
        if (this.pointRange.x < 0) return [];
        const nVertPerSeg = 4;
        const nseg = this.nseg;
        const ps = [];
        const offset = side % nVertPerSeg;
        // The wheel generator creates 'nseg' segments, each with 4 vertices.
        // Vertex offsets within a segment:
        // 0, 1: at angle 'i * angStep'
        // 2, 3: at angle '(i+1) * angStep'
        // To form a continuous path along the rim, we can pick offset 0 or 1.
        for (let i = 0; i < nseg; i++) {
            ps.push(this.pointRange.x + i * nVertPerSeg + offset);
        }
        console.log(`[SpaceCraft] Ring ${this.id} sideToPath(${side}): nseg=${nseg} ps=[${ps[0]}...${ps[ps.length-1]}]`);
        return ps;
    }
}

export class SpaceCraft {
    constructor() {
        this.nodes = [];
        this.girders = [];
        this.ropes = [];
        this.plates = [];
        this.rings = [];
        this.sliders = [];
    }

    clear() {
        this.nodes = [];
        this.girders = [];
        this.ropes = [];
        this.plates = [];
        this.sliders = [];
        this.rings = [];
    }

    addNode(pos, size) {
        const n = new Node(pos, size);
        n.id = this.nodes.length;
        this.nodes.push(n);
        return n;
    }

    addGirder(n1, n2, nseg, type) {
        const g = new Girder(n1, n2, nseg, type);
        g.id = this.girders.length;
        this.girders.push(g);
        return g;
    }

    addRope(n1, n2, thick, type) {
        const r = new Rope(n1, n2, thick, type);
        r.id = this.ropes.length;
        this.ropes.push(r);
        return r;
    }

    addPlate(boundA, boundB, spanA, spanB, type, kind = 'Radiator', nx = 2, ny = 2, nz = 1, upA = true, upB = true, sideOffset = 0, weldDist = 0) {
        if (!this.plates) this.plates = [];
        const p = new Plate(boundA, boundB, spanA, spanB, type, kind, nx, ny, nz, upA, upB, sideOffset, weldDist);
        p.id = this.plates.length;
        this.plates.push(p);
        return p;
    }

    addSlider(rail, sliding = null, calong = 0.0, matName = null, slidingVertId = -1, side = 0, methodFlag = true) {
        const s = new Slider(rail, sliding, calong, matName, slidingVertId, side, methodFlag);
        s.id = this.sliders.length;
        this.sliders.push(s);
        return s;
    }

    addRing(pos, dir, up, R, nseg, wh, matName, st) {
        const r = new Ring(pos, dir, up, R, nseg, wh, matName, st);
        r.id = this.rings.length;
        this.rings.push(r);
        return r;
    }

    addRing3P(p1, p2, p3, nseg, wh, matName, st) {
        const r = Ring.from3Points(p1, p2, p3, nseg, wh, matName, st);
        r.id = this.rings.length;
        this.rings.push(r);
        return r;
    }

    getStructuralComponent(id, kind) {
        if (kind === 1) return this.girders[id]; // ComponetKind::Girder
        if (kind === 3) return this.rings[id];   // ComponetKind::Ring
        return null;
    }
}
