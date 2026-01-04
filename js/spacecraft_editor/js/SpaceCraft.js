import { SpaceCraftWorkshop } from './SpaceCraftWorkshop.js';
import { Vec3 } from '../../common_js/Vec3.js';
import { bsplineInterpolate } from '../../common_js/SplineCubic.js';

// Abstract Component Classes

export class Node {
    constructor(pos, size = 1, type = 0) {
        this.pos = pos; // [x, y, z]
        this.size = size;
        this.type = type; // StickMaterial index
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
        this.type = type; // StickMaterial index
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

/**
 * Path represents an ordered sequence of vertices in the concrete mesh.
 * Multiple Sliders can share the same Path object to slide along the same edge.
 */
export class Path {
    constructor(rail, side = 0, methodFlag = true, radius = 0.35) {
        this.rail = rail;       // The structural component (Girder or Ring) providing the geometry
        this.side = side;       // The specific side/corner of the rail to follow
        this.methodFlag = methodFlag; // true = stride (exact vertices), false = SDF (nearest vertices)
        this.radius = radius;   // Picking radius for SDF method
        this.id = -1;           // Unique ID assigned by SpaceCraft
        this.ps = [];           // Concrete vertex indices in the mesh (populated during build)
        this.closed = (rail instanceof Ring); // Paths on rings are naturally circular
    }

    interpolate(t, mesh_verts) {
        return bsplineInterpolate(t, this.ps, this.closed, mesh_verts);
    }

    // /**
    //  * Interpolates a point along the path using cubic B-spline or linear interpolation.
    //  * @param {number} t - Interpolation parameter (0.0 to 1.0) along the whole path.
    //  * @param {Vec3[]} mesh_verts - The array of vertices from the concrete mesh.
    //  * @returns {Vec3} The interpolated position.
    //  */
    // interpolate(t, mesh_verts) {
    //     if (!this.ps || this.ps.length < 2) return new Vec3();
    //     const n = this.ps.length;
    //     const totalSegments = this.closed ? n : n - 1;
    //     const scaledT = t * totalSegments;
    //     let i = Math.floor(scaledT);
    //     let f = scaledT - i;

    //     if (i >= totalSegments) {
    //         if (this.closed) {
    //             i = i % n;
    //         } else {
    //             i = totalSegments - 1;
    //             f = 1.0;
    //         }
    //     }

    //     const getPt = (idx) => {
    //         if (this.closed) {
    //             idx = (idx % n + n) % n;
    //         } else {
    //             idx = Math.max(0, Math.min(n - 1, idx));
    //         }
    //         const vIdx = this.ps[idx];
    //         const v = mesh_verts ? mesh_verts[vIdx] : null;
    //         // MeshBuilder stores verts as objects with a `pos` Vec3; fall back to zero vec if missing
    //         const p = (v && v.pos) ? v.pos : v;
    //         if (!p || !isFinite(p.x) || !isFinite(p.y) || !isFinite(p.z)) {
    //             return new Vec3(); // safe fallback to avoid NaNs
    //         }
    //         return p;
    //     };

    //     // Cubic B-spline interpolation for smooth movement
    //     const p0 = getPt(i - 1);
    //     const p1 = getPt(i);
    //     const p2 = getPt(i + 1);
    //     const p3 = getPt(i + 2);

    //     const f2 = f * f;
    //     const f3 = f2 * f;

    //     // B-spline basis functions
    //     const b0 = (1 - 3 * f + 3 * f2 - f3) / 6.0;
    //     const b1 = (4 - 6 * f2 + 3 * f3) / 6.0;
    //     const b2 = (1 + 3 * f + 3 * f2 - 3 * f3) / 6.0;
    //     const b3 = f3 / 6.0;

    //     const res = new Vec3();
    //     res.x = p0.x * b0 + p1.x * b1 + p2.x * b2 + p3.x * b3;
    //     res.y = p0.y * b0 + p1.y * b1 + p2.y * b2 + p3.y * b3;
    //     res.z = p0.z * b0 + p1.z * b1 + p2.z * b2 + p3.z * b3;
    //     return res;
    // }
}

/**
 * Slider attaches a vertex to a Path.
 * It can either use an existing vertex from a 'sliding' component or generate a new one.
 */
export class Slider {
    constructor(path, sliding = null, calong = 0.0, matName = null, slidingVertId = -1) {
        this.path = path;       // The shared Path object this slider moves on
        this.sliding = sliding; // Optional: The component whose vertex is being attached (usually a Girder)
        this.slidingVertId = slidingVertId; // Optional: Specific vertex index if provided
        this.calong = calong;   // Interpolation parameter (0.0 to 1.0) along the path
        this.type = matName;    // Material name
        this.id = -1;           // Unique ID
        this.ivert = -1;        // The concrete vertex index in the mesh for this slider
        this.pAttach = null;    // IMPORTANT: If set, a NEW vertex is created at this [x,y,z] and welded to 'sliding'
        this.weldDist = 0.25;   // Tolerance for welding to the sliding component
    }
}

export class Ring {
    constructor(pos, dir, up, R, nseg, wh, matName, st, phase = 0.0) {
        this.pos = pos;
        this.dir = dir;
        this.up = up;
        this.R = R;
        this.nseg = nseg;
        this.wh = wh;
        this.matName = matName;
        this.st = st;
        this.phase = phase;
        this.id = -1;
        this.pointRange = { x: -1, y: -1 };
        this.stickRange = { x: -1, y: -1 };
    }

    static from3Points(p1, p2, p3, nseg, wh, matName, st, phase = 0.0) {
        const { center, radius, x, y } = Vec3.circle3Point(p1, p2, p3);
        const axis = new Vec3().setCross(x, y);
        axis.normalize();
        return new Ring(center, axis, x, radius, nseg, wh, matName, st, phase);
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
        this.paths = [];
        this.sliders = [];
        this.workshop = new SpaceCraftWorkshop();
    }

    clear() {
        this.nodes = [];
        this.girders = [];
        this.ropes = [];
        this.plates = [];
        this.paths = [];
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

    addPath(rail, side, methodFlag, radius) {
        const p = new Path(rail, side, methodFlag, radius);
        p.id = this.paths.length;
        this.paths.push(p);
        return p;
    }

    addSlider(path, sliding = null, calong = 0.0, matName = null, slidingVertId = -1) {
        const s = new Slider(path, sliding, calong, matName, slidingVertId);
        s.id = this.sliders.length;
        this.sliders.push(s);
        return s;
    }

    addRing(pos, dir, up, R, nseg, wh, matName, st, phase = 0.0) {
        const r = new Ring(pos, dir, up, R, nseg, wh, matName, st, phase);
        r.id = this.rings.length;
        this.rings.push(r);
        return r;
    }

    addRing3P(p1, p2, p3, nseg, wh, matName, st, phase = 0.0) {
        const r = Ring.from3Points(p1, p2, p3, nseg, wh, matName, st, phase);
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
