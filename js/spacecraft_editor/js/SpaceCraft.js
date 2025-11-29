// Abstract Component Classes

export class Node {
    constructor(pos, type = 0) {
        this.pos = pos; // [x, y, z]
        this.type = type;
        this.id = -1; // Assigned by Engine
    }
}

export class Girder {
    constructor(nodeA, nodeB, type = 0) {
        this.nodeA = nodeA; // Node object reference
        this.nodeB = nodeB; // Node object reference
        this.type = type;
        this.nseg = 1; // Default segments
        this.id = -1;
    }
}

export class Rope {
    constructor(nodeA, nodeB, thick, type = 0) {
        this.nodeA = nodeA;
        this.nodeB = nodeB;
        this.thick = thick;
        this.type = type;
        this.id = -1;
    }
}

export class SpaceCraft {
    constructor() {
        this.nodes = [];
        this.girders = [];
        this.ropes = [];
    }

    clear() {
        this.nodes = [];
        this.girders = [];
        this.ropes = [];
    }

    addNode(pos) {
        const n = new Node(pos);
        n.id = this.nodes.length;
        this.nodes.push(n);
        return n;
    }

    addGirder(n1, n2, type) {
        const g = new Girder(n1, n2, type);
        g.id = this.girders.length;
        this.girders.push(g);
        return g;
    }
}
