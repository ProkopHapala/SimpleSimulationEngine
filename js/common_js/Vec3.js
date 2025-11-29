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
        this.x = v.x; this.y = v.y; this.z = v.z;
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
