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
    constructor(rail, calong = 0.0, type = 0, sliding = null, slidingVertId = -1, side = 0, methodFlag = true) {
        // rail: structural component providing the path
        // sliding: structural component that owns the sliding vertex
        this.rail = rail;
        this.sliding = sliding;
        this.slidingVertId = slidingVertId; // vertex index on sliding component (within mesh)
        this.calong = calong;   // initial param guess along sliding component (e.g., along its girder)
        this.side = side;       // desired rail side (stride mode)
        this.methodFlag = methodFlag; // true=stride, false=SDF
        this.type = type;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
        this.path = { ps: [], cur: 0, closed: false };
        this.ivert = -1;        // actual slider vertex id in mesh (created during build)
        this.weldDist = 0.25;   // default weld radius to connect slider vertex to rail segment
        this.pathRadius = 0.35; // SDF cylinder radius for path picking
    }
}

export class Ring {
    constructor(pos, rot, R, type = 0) {
        this.pos = pos;
        this.rot = rot; // Mat3
        this.R = R;
        this.type = type;
        this.nseg = 16;
        this.wh = { x: 0.5, y: 0.5 };
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }

    sideToPath(side) {
        // Placeholder for path extraction logic
        return [];
    }
}

export class SpaceCraft {
    constructor() {
        this.nodes = [];
        this.girders = [];
        this.ropes = [];
        this.plates = [];
        this.sliders = [];
        this.rings = [];
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

    addSlider(rail, calong, type, sliding = null, slidingVertId = -1, side = 0, methodFlag = true) {
        const s = new Slider(rail, calong, type, sliding, slidingVertId, side, methodFlag);
        s.id = this.sliders.length;
        this.sliders.push(s);
        return s;
    }

    addRing(pos, rot, R, type) {
        const r = new Ring(pos, rot, R, type);
        r.id = this.rings.length;
        this.rings.push(r);
        return r;
    }
}
