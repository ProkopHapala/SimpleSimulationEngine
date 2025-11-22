
// SpaceCraft2Mesh.js

function BuildCraft_truss(mesh, craft) {
    mesh.clear();
    mesh.block(); // Start initial block

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
        mesh.edge(ivA, ivB, girder.type);

        // Example of using sub-segments (commented out for now until we want high-res)
        // const pA = girder.nodeA.pos;
        // const pB = girder.nodeB.pos;
        // mesh.girder1(pA, pB, [0,1,0], girder.nseg, 0.1, girder.type);
    }

    // 3. Ropes
    // ... similar to girders
}
