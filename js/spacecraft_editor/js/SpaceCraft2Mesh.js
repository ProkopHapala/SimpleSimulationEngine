import { Vec3 } from '../../common_js/Vec3.js';

// SpaceCraft2Mesh.js

export function BuildCraft_truss(mesh, craft) {
    mesh.clear();
    // 1. Nodes
    // In the C++ version, nodes might generate mesh geometry (like a connector hub).
    // For now, we just ensure they exist as vertices if they aren't already.
    // But actually, Girders need to connect to these nodes.

    // Strategy:
    // We map each Abstract Node to a Vertex Index in the MeshBuilder.
    const nodeToVertMap = new Map();

    for (const node of craft.nodes) {
        const iv = mesh.vert(node.pos);
        nodeToVertMap.set(node, iv);
    }

    // 2. Girders
    for (const girder of craft.girders) {
        const ivA = nodeToVertMap.get(girder.nodeA);
        const ivB = nodeToVertMap.get(girder.nodeB);

        // For simple visualization, just connect them with an edge
        // In the future, use mesh.girder1() to generate sub-segments
        // mesh.edge(ivA, ivB, girder.type);

        // Use girder1 for visible truss structure
        const pA = new Vec3(girder.nodeA.pos[0], girder.nodeA.pos[1], girder.nodeA.pos[2]);
        const pB = new Vec3(girder.nodeB.pos[0], girder.nodeB.pos[1], girder.nodeB.pos[2]);
        const up = new Vec3(0, 1, 0);
        const stickTypes = { x: 1, y: 1, z: 1, w: 1 }; // Default stick types

        mesh.girder1(pA, pB, up, girder.nseg, 0.1, stickTypes, true);
    }

    // 3. Ropes
    // ... similar to girders
}
