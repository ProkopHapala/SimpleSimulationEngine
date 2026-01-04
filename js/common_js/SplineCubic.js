import { Vec3 } from './Vec3.js';

/**
 * Interpolates a point along a path using cubic B-spline.
 * Pure function: no reliance on `this`.
 * @param {number} t          - Interpolation parameter (0..1) along the whole path.
 * @param {number[]} ps       - Array of vertex indices forming the path.
 * @param {boolean} closed    - Whether the path is closed (wrap indices).
 * @param {Object[]} verts    - MeshBuilder vertex array (entries may be {pos:Vec3} or Vec3).
 * @returns {Vec3}            - Interpolated position (Vec3).
 */
export function bsplineInterpolate(t, ps, closed, verts) {
    if (!ps || ps.length < 2 || !verts) return new Vec3();
    const n = ps.length;
    const totalSegments = closed ? n : n - 1;
    if (totalSegments <= 0) return new Vec3();

    const scaledT = t * totalSegments;
    let i = Math.floor(scaledT);
    let f = scaledT - i;

    if (i >= totalSegments) {
        if (closed) {
            i = i % n;
        } else {
            i = totalSegments - 1;
            f = 1.0;
        }
    }

    const getPt = (idx) => {
        if (closed) {
            idx = (idx % n + n) % n;
        } else {
            idx = Math.max(0, Math.min(n - 1, idx));
        }
        const vIdx = ps[idx];
        const v = verts[vIdx];
        const p = (v && v.pos) ? v.pos : v;
        if (!p || !isFinite(p.x) || !isFinite(p.y) || !isFinite(p.z)) {
            return new Vec3(); // safe fallback
        }
        return p;
    };

    // Cubic B-spline control points
    const p0 = getPt(i - 1);
    const p1 = getPt(i);
    const p2 = getPt(i + 1);
    const p3 = getPt(i + 2);

    const f2 = f * f;
    const f3 = f2 * f;

    // B-spline basis functions
    const b0 = (1 - 3*f  + 3 * f2 - f3     ) / 6.0;
    const b1 = (4 - 6*f2 + 3 * f3          ) / 6.0;
    const b2 = (1 + 3*f  + 3 * f2 - 3 * f3 ) / 6.0;
    const b3 = f3                            / 6.0;

    const res = new Vec3();
    res.x = p0.x * b0 + p1.x * b1 + p2.x * b2 + p3.x * b3;
    res.y = p0.y * b0 + p1.y * b1 + p2.y * b2 + p3.y * b3;
    res.z = p0.z * b0 + p1.z * b1 + p2.z * b2 + p3.z * b3;
    return res;
}