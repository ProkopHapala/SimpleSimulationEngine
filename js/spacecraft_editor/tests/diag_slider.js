
const { Vec3 } = require('../../common_js/Vec3.js');
const { MeshBuilder } = require('../../common_js/MeshBuilder.js');
const { SDF_Cylinder } = require('../../common_js/SDfuncs.js');

// Mock for local testing
const Girder = class {
    constructor(id, nA, nB, nseg, up) {
        this.id = id;
        this.nodeA = nA;
        this.nodeB = nB;
        this.nseg = nseg;
        this.up = up;
        this.pointRange = { x: -1, y: -1 };
    }
};

const Node = class {
    constructor(id, pos, size) {
        this.id = id;
        this.pos = pos;
        this.size = size;
        this.chunkRange = { x: -1, y: -1 };
    }
};

const cornerDirs = [
    { a: 1, b: 1 },
    { a: 1, b: -1 },
    { a: -1, b: 1 },
    { a: -1, b: -1 },
];

function testSliderSelection() {
    const mesh = new MeshBuilder();
    mesh.clear();

    // 1. Setup mock girder
    const n1 = new Node(0, [0, 0, 0], 1.0);
    const n2 = new Node(1, [10, 0, 0], 1.0);
    
    // Add node blocks
    n1.chunkRange = mesh.addCube(new Vec3(0, 0, 0), 1.0, true);
    n2.chunkRange = mesh.addCube(new Vec3(10, 0, 0), 1.0, true);
    
    const iv0 = mesh.verts.length;
    // Bridge them
    mesh.bridgeFacingPolygons(new Vec3(0,0,0), new Vec3(10,0,0), n1.chunkRange, n2.chunkRange, 5, {x:1,y:1,z:1,w:1}, {x:1,y:1,z:1,w:1});
    const iv1 = mesh.verts.length;
    
    const girder = new Girder(0, n1, n2, 5, [0, 1, 0]);
    girder.pointRange = { x: iv0, y: iv1 };

    console.log(`Girder pointRange: [${iv0}, ${iv1}) total: ${iv1 - iv0}`);

    const side = 0; // (+side, +up) corner
    const radius = 0.35;
    
    // 2. Compute Offset
    const pA = new Vec3(0, 0, 0);
    const pB = new Vec3(10, 0, 0);
    const dirN = new Vec3().setSub(pB, pA);
    dirN.normalize();
    const up = new Vec3(0, 1, 0);
    const sideVec = new Vec3().setCross(dirN, up);
    sideVec.normalize();
    const cd = cornerDirs[side % 4];
    
    // Corner positions for side 0 should be at +/- 0.5 (half-size)
    const offMag = 0.5; 
    const offset = new Vec3().setV(up).mulScalar(cd.b * offMag).addMul(sideVec, cd.a * offMag);
    
    const p0 = new Vec3().setV(pA).add(offset);
    const p1 = new Vec3().setV(pB).add(offset);
    
    console.log(`Offset: ${offset.x}, ${offset.y}, ${offset.z}`);
    console.log(`p0: ${p0.x}, ${p0.y}, ${p0.z}`);
    console.log(`p1: ${p1.x}, ${p1.y}, ${p1.z}`);

    // 3. SDF Selection
    const selRadius = 0.05; 
    const sdf = SDF_Cylinder(p0, p1, selRadius, true);
    mesh.setSelectionKind('vert');
    mesh.applySelection([], 'vert', false);
    const nsel = mesh.selectVertsBySDF(sdf, 0.05, true);
    
    const selectionArray = mesh.selection.vec.filter(i => i >= 0);
    const psIdx = selectionArray.filter(i => i >= girder.pointRange.x && i < girder.pointRange.y);
    
    console.log(`SDF Selection (R=0.05, T=0.05) - total in mesh: ${nsel}, in girder range: ${psIdx.length}`);
    console.log(`Global selection indices: ${JSON.stringify(selectionArray)}`);
    
    console.log("Girder vertices details:");
    for (let i = girder.pointRange.x; i < girder.pointRange.y; i++) {
        const v = mesh.verts[i].pos;
        const d = sdf(v);
        console.log(`  Vert ${i}: pos=(${v.x.toFixed(2)}, ${v.y.toFixed(2)}, ${v.z.toFixed(2)}) sdfDist=${d.toFixed(4)}`);
    }

    if (psIdx.length === 0) {
        console.log("FAILED: No vertices selected from girder range");
    } else if (psIdx.length > (girder.nseg + 1)) {
        console.log("FAILED: Selected too many vertices (multiple edges?) count=" + psIdx.length);
    } else {
        console.log("SUCCESS: Selected single edge vertices count=" + psIdx.length);
    }
}

testSliderSelection();
