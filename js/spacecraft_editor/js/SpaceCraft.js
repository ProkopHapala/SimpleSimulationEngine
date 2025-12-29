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
    constructor(boundA, boundB, spanA = [0, 1], spanB = [0, 1], type = 0, kind = 'Radiator') {
        this.boundA = boundA; // Girder or Rope
        this.boundB = boundB; // Girder or Rope
        this.spanA = spanA;   // [cmin,cmax] along boundA
        this.spanB = spanB;   // [cmin,cmax] along boundB
        this.kind = kind;     // 'Radiator' or 'Shield'
        this.type = type;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }
}

export class Slider {
    constructor(boundTo, calong, type = 0) {
        this.boundTo = boundTo; // StructuralComponent reference
        this.calong = calong;   // position along 0..1
        this.type = type;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
        this.path = { ps: [], cur: 0, closed: false };
        this.ivert = -1;
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

    addPlate(boundA, boundB, spanA, spanB, type, kind = 'Radiator') {
        if (!this.plates) this.plates = [];
        const p = new Plate(boundA, boundB, spanA, spanB, type, kind);
        p.id = this.plates.length;
        this.plates.push(p);
        return p;
    }

    addSlider(boundTo, calong, type) {
        const s = new Slider(boundTo, calong, type);
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
