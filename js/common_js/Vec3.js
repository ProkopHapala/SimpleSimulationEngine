export class Vec3 {
    constructor(x = 0.0, y = 0.0, z = 0.0) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    // ================= Basic Setters =================

    set(x, y, z) {
        this.x = x; this.y = y; this.z = z;
        return this;
    }

    setV(v) {
        if (!v) return this;
        this.x = v.x ?? v[0] ?? 0;
        this.y = v.y ?? v[1] ?? 0;
        this.z = v.z ?? v[2] ?? 0;
        return this;
    }

    clone() {
        return new Vec3(this.x, this.y, this.z);
    }

    // ================= Arithmetic (In-Place) =================

    add(v) { this.x += v.x; this.y += v.y; this.z += v.z; return this; }
    sub(v) { this.x -= v.x; this.y -= v.y; this.z -= v.z; return this; }
    mul(v) { this.x *= v.x; this.y *= v.y; this.z *= v.z; return this; }
    div(v) { this.x /= v.x; this.y /= v.y; this.z /= v.z; return this; }

    addScalar(f) { this.x += f; this.y += f; this.z += f; return this; }
    subScalar(f) { this.x -= f; this.y -= f; this.z -= f; return this; }
    mulScalar(f) { this.x *= f; this.y *= f; this.z *= f; return this; }
    divScalar(f) { this.x /= f; this.y /= f; this.z /= f; return this; }

    // ================= Fused Operations (Optimization) =================

    addMul(v, f) {
        this.x += v.x * f;
        this.y += v.y * f;
        this.z += v.z * f;
        return this;
    }

    subMul(v, f) {
        this.x -= v.x * f;
        this.y -= v.y * f;
        this.z -= v.z * f;
        return this;
    }

    setAdd(a, b) {
        this.x = a.x + b.x;
        this.y = a.y + b.y;
        this.z = a.z + b.z;
        return this;
    }

    setSub(a, b) {
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        this.z = a.z - b.z;
        return this;
    }

    setAdd(a, b) {
        this.x = a.x + b.x;
        this.y = a.y + b.y;
        this.z = a.z + b.z;
        return this;
    }

    setMul(a, b) {
        this.x = a.x * b.x;
        this.y = a.y * b.y;
        this.z = a.z * b.z;
        return this;
    }

    setAddMul(a, b, f) {
        this.x = a.x + b.x * f;
        this.y = a.y + b.y * f;
        this.z = a.z + b.z * f;
        return this;
    }

    setLincomb(fa, a, fb, b) {
        this.x = fa * a.x + fb * b.x;
        this.y = fa * a.y + fb * b.y;
        this.z = fa * a.z + fb * b.z;
        return this;
    }

    addLincomb(fa, a, fb, b) {
        this.x += fa * a.x + fb * b.x;
        this.y += fa * a.y + fb * b.y;
        this.z += fa * a.z + fb * b.z;
        return this;
    }

    setLincomb3(fa, a, fb, b, fc, c) {
        this.x = fa * a.x + fb * b.x + fc * c.x;
        this.y = fa * a.y + fb * b.y + fc * c.y;
        this.z = fa * a.z + fb * b.z + fc * c.z;
        return this;
    }

    // ================= Geometric =================

    dot(v) {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }

    norm2() {
        return this.x * this.x + this.y * this.y + this.z * this.z;
    }

    norm() {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }

    normalize() {
        const l2 = this.x * this.x + this.y * this.y + this.z * this.z;
        if (l2 < 1e-32) return 0.0;
        const l = Math.sqrt(l2);
        const invL = 1.0 / l;
        this.x *= invL;
        this.y *= invL;
        this.z *= invL;
        return l;
    }

    setCross(a, b) {
        const ax = a.x, ay = a.y, az = a.z;
        const bx = b.x, by = b.y, bz = b.z;
        this.x = ay * bz - az * by;
        this.y = az * bx - ax * bz;
        this.z = ax * by - ay * bx;
        return this;
    }

    addCross(a, b) {
        this.x += a.y * b.z - a.z * b.y;
        this.y += a.z * b.x - a.x * b.z;
        this.z += a.x * b.y - a.y * b.x;
        return this;
    }

    // ================= Advanced Geometric =================

    makeOrtho(a) {
        const dot = this.x * a.x + this.y * a.y + this.z * a.z;
        const norm2 = a.x * a.x + a.y * a.y + a.z * a.z;
        const c = dot / norm2;
        this.x -= a.x * c;
        this.y -= a.y * c;
        this.z -= a.z * c;
        return c;
    }

    makeOrthoU(a) {
        const c = this.x * a.x + this.y * a.y + this.z * a.z;
        this.x -= a.x * c;
        this.y -= a.y * c;
        this.z -= a.z * c;
        return c;
    }

    // ================= Rotation =================

    rotateCSA(ca, sa, uaxis) {
        const dot = this.x * uaxis.x + this.y * uaxis.y + this.z * uaxis.z;
        const cu = (1 - ca) * dot;

        const utx = uaxis.y * this.z - uaxis.z * this.y;
        const uty = uaxis.z * this.x - uaxis.x * this.z;
        const utz = uaxis.x * this.y - uaxis.y * this.x;

        const nx = ca * this.x + sa * utx + cu * uaxis.x;
        const ny = ca * this.y + sa * uty + cu * uaxis.y;
        const nz = ca * this.z + sa * utz + cu * uaxis.z;

        this.x = nx; this.y = ny; this.z = nz;
        return this;
    }

    rotate(angle, axis) {
        const ca = Math.cos(angle);
        const sa = Math.sin(angle);
        return this.rotateCSA(ca, sa, axis);
    }

    getSomePerp() {
        // returns a normalized vector perpendicular to this
        const ax = Math.abs(this.x), ay = Math.abs(this.y);
        const ref = (ax < 0.9 && ay < 0.9) ? new Vec3(1, 0, 0) : new Vec3(0, 0, 1);
        const perp = new Vec3().setCross(this, ref);
        const n = perp.norm();
        if (n < 1e-12) return new Vec3(0, 1, 0);
        return perp.mulScalar(1 / n);
    }
    getSomeOrtho() {
        // return a normalized vector orthogonal to this, choosing the most stable axis
        const ax = Math.abs(this.x), ay = Math.abs(this.y), az = Math.abs(this.z);
        const ref = (ax < ay)
            ? ((ax < az) ? new Vec3(1, 0, 0) : new Vec3(0, 0, 1))
            : ((ay < az) ? new Vec3(0, 1, 0) : new Vec3(0, 0, 1));
        const perp = new Vec3().setCross(this, ref);
        const n = perp.norm();
        if (n < 1e-12) return new Vec3(1, 0, 0);
        return perp.mulScalar(1 / n);
    }
    /**
     * Compute circle from 3 points
     * @param {Vec3} A - point 1
     * @param {Vec3} B - point 2
     * @param {Vec3} C - point 3
     * @returns {Object} { center: Vec3, radius: number, x: Vec3, y: Vec3 }
     */
    static circle3Point(A, B, C) {
        const vA = new Vec3().setV(A);
        const vB = new Vec3().setV(B);
        const vC = new Vec3().setV(C);
        // center = (A+B)*0.5;
        let center = new Vec3().setAdd(vA, vB).mulScalar(0.5);
        // a = B-A; double xa = a.normalize()/2;
        let x = new Vec3().setSub(vB, vA);
        let xa = x.normalize() / 2;
        // b = C-center; xc = b.dot(a); b.add_mul(a, -xc); double yc = b.normalize();
        let y = new Vec3().setSub(vC, center);
        let xc = y.dot(x);
        y.addMul(x, -xc);
        let yc = y.normalize();
        if (yc < 1e-12) {
            // nearly collinear, fallback to zero-radius at midpoint
            return { center, radius: 0, x: new Vec3(1, 0, 0), y: new Vec3(0, 1, 0) };
        }
        // y_coord = (xa^2 - yc^2 - xc^2)/(-2*yc)
        let y_coord = (xa * xa - yc * yc - xc * xc) / (-2 * yc);
        center.addMul(y, y_coord);
        let radius = Math.sqrt(xa * xa + y_coord * y_coord);
        return { center, radius, x, y };
    }

    // ================= Utility =================

    dist2(v) {
        const dx = this.x - v.x;
        const dy = this.y - v.y;
        const dz = this.z - v.z;
        return dx * dx + dy * dy + dz * dz;
    }

    toString() {
        return `(${this.x.toFixed(6)}, ${this.y.toFixed(6)}, ${this.z.toFixed(6)})`;
    }
}

// Also expose on window for legacy global-script users
if (typeof window !== 'undefined') {
    window.Vec3 = Vec3;
}
