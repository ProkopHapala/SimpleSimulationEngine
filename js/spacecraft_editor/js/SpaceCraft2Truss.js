
import { Vec3 } from '../../common_js/Vec3.js';

export function exportToTruss(sim, mesh, workshop) {
    const np = mesh.verts.length;
    const nb = mesh.edges.length;
    console.log(`exportToTruss() START nvert=${np} nedge=${nb}`);

    // sim.recalloc(np, nneighmax, nb) - assuming sim has similar API to TrussDynamics
    // For JS implementation, we might need to adjust based on how sim is implemented.
    // If sim is not yet fully implemented, we'll assume a basic structure.

    if (sim.recalloc) {
        // Find max neighbors
        const nneighs = new Int32Array(np).fill(0);
        for (const e of mesh.edges) {
            nneighs[e.x]++;
            nneighs[e.y]++;
        }
        let nneighmax = 0;
        for (let i = 0; i < np; i++) if (nneighs[i] > nneighmax) nneighmax = nneighs[i];
        sim.recalloc(np, nneighmax, nb);
    }

    // Set points and clear masses
    for (let i = 0; i < np; i++) {
        const v = mesh.verts[i];
        if (sim.setPoint) {
            sim.setPoint(i, v.pos.x, v.pos.y, v.pos.z, 0);
        } else if (sim.points) {
            // If points is a Float64Array or similar
            const offset = i * 4; // assuming [x, y, z, mass]
            sim.points[offset] = v.pos.x;
            sim.points[offset+1] = v.pos.y;
            sim.points[offset+2] = v.pos.z;
            sim.points[offset+3] = 0;
        }
    }

    // Process edges
    for (let i = 0; i < nb; i++) {
        const e = mesh.edges[i];
        const type = e.z; // edge.z is used for type in MeshBuilder.js
        if (type < 0 || type >= workshop.stickMaterials.length) {
            console.error(`ERROR: mesh.edges[${i}].type=${type} out of bounds`);
            continue;
        }

        const mat = workshop.stickMaterials[type];
        const p1 = mesh.verts[e.x].pos;
        const p2 = mesh.verts[e.y].pos;
        
        const dx = p2.x - p1.x;
        const dy = p2.y - p1.y;
        const dz = p2.z - p1.z;
        const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
        const l0 = dist * (1.0 - mat.preStrain);
        const mass = l0 * mat.linearDensity;

        // Distribute mass to endpoints
        if (sim.addMass) {
            sim.addMass(e.x, mass * 0.5);
            sim.addMass(e.y, mass * 0.5);
        } else if (sim.points) {
            sim.points[e.x * 4 + 3] += mass * 0.5;
            sim.points[e.y * 4 + 3] += mass * 0.5;
        }

        const kPush = mat.Kpush / l0;
        const kPull = mat.Kpull / l0;
        const damping = mat.damping;

        if (sim.setBond) {
            sim.setBond(i, e.x, e.y, l0, kPush, kPull, damping);
        }
        
        // Additional parameters if needed (strength/strain)
        if (sim.setBondStrain) {
            sim.setBondStrain(i, mat.Spull / mat.Kpull, mat.Spush / mat.Kpush);
        }
    }

    if (sim.checkMasses) sim.checkMasses();
    if (sim.cleanForce) sim.cleanForce();
    if (sim.cleanVel) sim.cleanVel();

    console.log("exportToTruss() DONE!");
}
