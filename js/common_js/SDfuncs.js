// JS counterparts of C++ SDfuncs.h
// Provides simple SDF helpers that return callables taking a Vec3 and
// returning a scalar distance / value, suitable for use with
// MeshBuilder.selectVertsBySDF / selectEdgesBySDF.

import { Vec3 } from './Vec3.js';

// Helper to ensure we always work with Vec3 instances
function toVec3(p) {
    if (p instanceof Vec3) return p;
    if (Array.isArray(p)) return new Vec3(p[0], p[1], p[2]);
    if (p && typeof p.x === 'number') return new Vec3(p.x, p.y, p.z);
    return new Vec3(0, 0, 0);
}

// --- SDF_point2 : squared distance to a point ---
// C++: returns (p - center).norm2()
export function SDF_point2(center) {
    const c = toVec3(center);
    return function sdf_point2(p) {
        const v = toVec3(p);
        const dx = v.x - c.x;
        const dy = v.y - c.y;
        const dz = v.z - c.z;
        return dx * dx + dy * dy + dz * dz;
    };
}

// --- SDF_Sphere ---
// C++: (p - center).norm() - radius
export function SDF_Sphere(center, radius) {
    const c = toVec3(center);
    const r = radius;
    return function sdf_sphere(p) {
        const v = toVec3(p);
        const dx = v.x - c.x;
        const dy = v.y - c.y;
        const dz = v.z - c.z;
        const dist = Math.hypot(dx, dy, dz);
        return dist - r;
    };
}

// --- SDF_AABB ---
// C++:
//   Vec3d  q = (p-center).abs() - halfSpan;
//   double l = q.maxComponent(); if(l< 0.0) return l;
//   return q.max(Vec3dZero).norm() + l;
export function SDF_AABB(center, halfSpan) {
    const c = toVec3(center);
    const h = toVec3(halfSpan);
    return function sdf_aabb(p) {
        const v = toVec3(p);
        // q = |p - c| - h
        let qx = Math.abs(v.x - c.x) - h.x;
        let qy = Math.abs(v.y - c.y) - h.y;
        let qz = Math.abs(v.z - c.z) - h.z;

        const l = Math.max(qx, qy, qz);
        if (l < 0.0) return l;

        // max(q, 0)
        qx = Math.max(qx, 0.0);
        qy = Math.max(qy, 0.0);
        qz = Math.max(qz, 0.0);
        const qlen = Math.hypot(qx, qy, qz);
        return qlen + l;
    };
}

// --- SDF_Cylinder ---
// C++ stores p0, hdir (unit), l (length), r, capped and in operator():
//   q = p - p0; h = q.dot(hdir);
//   if (h<0) or (h>l) handle caps or infinite-cylinder distance;
//   else distance to radial direction.
export function SDF_Cylinder(p0, p1, radius, capped = true) {
    const base = toVec3(p0);
    const top = toVec3(p1);
    const axis = new Vec3().setSub(top, base);
    const length = axis.normalize(); // axis becomes unit, length stored separately
    const r = radius;
    const bCapped = !!capped;

    return function sdf_cylinder(p) {
        const v = toVec3(p);
        const qx = v.x - base.x;
        const qy = v.y - base.y;
        const qz = v.z - base.z;
        const h = qx * axis.x + qy * axis.y + qz * axis.z;
        let dist;

        if (h < 0.0) {
            if (bCapped) {
                // distance to sphere around p0 minus radius
                const dx = qx;
                const dy = qy;
                const dz = qz;
                dist = Math.hypot(dx, dy, dz) - r;
            } else {
                dist = -h; // distance along axis outside infinite cylinder
            }
        } else if (h > length) {
            if (bCapped) {
                // point projected beyond top cap at p0 + axis*length
                const px = base.x + axis.x * length;
                const py = base.y + axis.y * length;
                const pz = base.z + axis.z * length;
                const dx = v.x - px;
                const dy = v.y - py;
                const dz = v.z - pz;
                dist = Math.hypot(dx, dy, dz) - r;
            } else {
                dist = h - length;
            }
        } else {
            // between caps: distance to axis minus radius
            const px = base.x + axis.x * h;
            const py = base.y + axis.y * h;
            const pz = base.z + axis.z * h;
            const dx = v.x - px;
            const dy = v.y - py;
            const dz = v.z - pz;
            dist = Math.hypot(dx, dy, dz) - r;
        }
        return dist;
    };
}

// --- SDF_ScreenRectWorld ---
// Screen-space rectangle SDF operating on world-space positions.
//
// This mirrors the inline implementation that was previously in GUI.js:
//  - project world pos with camera into NDC
//  - treat [minXN,maxXN]x[minYN,maxYN] as the selection box in NDC
//  - return negative inside, positive outside, ~0 on the boundary.
//
// NOTE (current implementation vs desired architecture):
//   - Right now this uses THREE.Vector3().project(camera) internally, which
//     makes this helper depend on THREE being present in the environment.
//   - A more engine-agnostic version should take a precomputed view-projection
//     matrix (and, if needed, camera origin / ray parameters) and perform the
//     projection with plain matrix math on Vec3Ref, removing the hard THREE
//     dependency. For orthographic projection (which SpaceCraft editor uses),
//     this reduces to a simple linear transform to clip-space/NDC.
//   - This function is already structured so that the expensive rectangle and
//     NDC bounds setup is done once when SDF_ScreenRectWorld(...) is called;
//     the returned closure only does per-vertex projection + a few branches.
//
// Parameters:
//   camera : THREE.Camera (must support .project on a THREE.Vector3)
//   rect   : DOMRect of the canvas getBoundingClientRect()
//   minX,maxX,minY,maxY : screen-space box in client coordinates
export function SDF_ScreenRectWorld(camera, rect, minX, maxX, minY, maxY) {
    if (!camera || !rect) {
        // Degenerate SDF: everything is "far away" so nothing is selected
        return function sdf_disable() { return 1e9; };
    }

    const width = rect.width;
    const height = rect.height;
    const left = rect.left;
    const top = rect.top;

    const toNDC = (sx, sy) => {
        const xndc = ((sx - left) / width) * 2.0 - 1.0;
        const yndc = 1.0 - 2.0 * ((sy - top) / height);
        return { x: xndc, y: yndc };
    };

    const p00 = toNDC(minX, minY);
    const p11 = toNDC(maxX, maxY);
    const minXN = Math.min(p00.x, p11.x);
    const maxXN = Math.max(p00.x, p11.x);
    const minYN = Math.min(p00.y, p11.y);
    const maxYN = Math.max(p00.y, p11.y);

    // Use THREE.Vector3 if available (expected in browser builds).
    const tmp = (typeof THREE !== 'undefined' && THREE.Vector3)
        ? new THREE.Vector3()
        : null;

    return function sdf_screenRect(pos) {
        if (!tmp) {
            // No projection available (e.g. Node without THREE): disable.
            return 1e9;
        }

        tmp.set(pos.x, pos.y, pos.z).project(camera);
        const x = tmp.x;
        const y = tmp.y;

        const dx = (x < minXN) ? (minXN - x) : (x > maxXN ? x - maxXN : 0.0);
        const dy = (y < minYN) ? (minYN - y) : (y > maxYN ? y - maxYN : 0.0);
        if (dx === 0.0 && dy === 0.0) return -1.0;
        return Math.hypot(dx, dy);
    };
}

// Aggregate export for convenience and legacy global wiring
export const SDfuncs = {
    SDF_point2,
    SDF_Sphere,
    SDF_AABB,
    SDF_Cylinder,
    SDF_ScreenRectWorld,
};

// Also expose on window for legacy global-script users
if (typeof window !== 'undefined') {
    window.SDfuncs = SDfuncs;
}
