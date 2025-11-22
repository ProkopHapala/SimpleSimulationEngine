
export class Vec3 extends Float64Array {
    constructor(x = 0.0, y = 0.0, z = 0.0) {
        super(3);
        this[0] = x;
        this[1] = y;
        this[2] = z;
    }

    // ================= Accessors =================
    get x() { return this[0]; } set x(v) { this[0] = v; }
    get y() { return this[1]; } set y(v) { this[1] = v; }
    get z() { return this[2]; } set z(v) { this[2] = v; }

    // ================= Basic Setters =================

    set(x, y, z) {
        this[0] = x; this[1] = y; this[2] = z;
        return this;
    }

    setV(v) {
        this[0] = v[0]; this[1] = v[1]; this[2] = v[2];
        return this;
    }

    clone() {
        return new Vec3(this[0], this[1], this[2]);
    }

    // ================= Arithmetic (In-Place) =================

    add(v) { this[0] += v[0]; this[1] += v[1]; this[2] += v[2]; return this; }
    sub(v) { this[0] -= v[0]; this[1] -= v[1]; this[2] -= v[2]; return this; }
    mul(v) { this[0] *= v[0]; this[1] *= v[1]; this[2] *= v[2]; return this; }
    div(v) { this[0] /= v[0]; this[1] /= v[1]; this[2] /= v[2]; return this; }

    addScalar(f) { this[0] += f; this[1] += f; this[2] += f; return this; }
    subScalar(f) { this[0] -= f; this[1] -= f; this[2] -= f; return this; }
    mulScalar(f) { this[0] *= f; this[1] *= f; this[2] *= f; return this; }
    divScalar(f) { this[0] /= f; this[1] /= f; this[2] /= f; return this; }

    // ================= Fused Operations (Optimization) =================

    addMul(v, f) {
        this[0] += v[0] * f;
        this[1] += v[1] * f;
        this[2] += v[2] * f;
        return this;
    }

    subMul(v, f) {
        this[0] -= v[0] * f;
        this[1] -= v[1] * f;
        this[2] -= v[2] * f;
        return this;
    }

    setAdd(a, b) {
        this[0] = a[0] + b[0];
        this[1] = a[1] + b[1];
        this[2] = a[2] + b[2];
        return this;
    }

    setSub(a, b) {
        this[0] = a[0] - b[0];
        this[1] = a[1] - b[1];
        this[2] = a[2] - b[2];
        return this;
    }

    setMul(a, b) {
        this[0] = a[0] * b[0];
        this[1] = a[1] * b[1];
        this[2] = a[2] * b[2];
        return this;
    }

    setAddMul(a, b, f) {
        this[0] = a[0] + b[0] * f;
        this[1] = a[1] + b[1] * f;
        this[2] = a[2] + b[2] * f;
        return this;
    }

    setLincomb(fa, a, fb, b) {
        this[0] = fa * a[0] + fb * b[0];
        this[1] = fa * a[1] + fb * b[1];
        this[2] = fa * a[2] + fb * b[2];
        return this;
    }

    addLincomb(fa, a, fb, b) {
        this[0] += fa * a[0] + fb * b[0];
        this[1] += fa * a[1] + fb * b[1];
        this[2] += fa * a[2] + fb * b[2];
        return this;
    }

    setLincomb3(fa, a, fb, b, fc, c) {
        this[0] = fa * a[0] + fb * b[0] + fc * c[0];
        this[1] = fa * a[1] + fb * b[1] + fc * c[1];
        this[2] = fa * a[2] + fb * b[2] + fc * c[2];
        return this;
    }

    // ================= Geometric =================

    dot(v) {
        return this[0] * v[0] + this[1] * v[1] + this[2] * v[2];
    }

    norm2() {
        return this[0] * this[0] + this[1] * this[1] + this[2] * this[2];
    }

    norm() {
        return Math.sqrt(this[0] * this[0] + this[1] * this[1] + this[2] * this[2]);
    }

    normalize() {
        const l2 = this[0] * this[0] + this[1] * this[1] + this[2] * this[2];
        if (l2 < 1e-32) return 0.0;
        const l = Math.sqrt(l2);
        const invL = 1.0 / l;
        this[0] *= invL;
        this[1] *= invL;
        this[2] *= invL;
        return l;
    }

    setCross(a, b) {
        const ax = a[0], ay = a[1], az = a[2];
        const bx = b[0], by = b[1], bz = b[2];
        this[0] = ay * bz - az * by;
        this[1] = az * bx - ax * bz;
        this[2] = ax * by - ay * bx;
        return this;
    }

    addCross(a, b) {
        this[0] += a[1] * b[2] - a[2] * b[1];
        this[1] += a[2] * b[0] - a[0] * b[2];
        this[2] += a[0] * b[1] - a[1] * b[0];
        return this;
    }

    // ================= Advanced Geometric =================

    makeOrtho(a) {
        const dot = this[0] * a[0] + this[1] * a[1] + this[2] * a[2];
        const norm2 = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
        const c = dot / norm2;
        this[0] -= a[0] * c;
        this[1] -= a[1] * c;
        this[2] -= a[2] * c;
        return c;
    }

    makeOrthoU(a) {
        const c = this[0] * a[0] + this[1] * a[1] + this[2] * a[2];
        this[0] -= a[0] * c;
        this[1] -= a[1] * c;
        this[2] -= a[2] * c;
        return c;
    }

    // ================= Rotation =================

    rotateCSA(ca, sa, uaxis) {
        const dot = this[0] * uaxis[0] + this[1] * uaxis[1] + this[2] * uaxis[2];
        const cu = (1 - ca) * dot;

        const utx = uaxis[1] * this[2] - uaxis[2] * this[1];
        const uty = uaxis[2] * this[0] - uaxis[0] * this[2];
        const utz = uaxis[0] * this[1] - uaxis[1] * this[0];

        const nx = ca * this[0] + sa * utx + cu * uaxis[0];
        const ny = ca * this[1] + sa * uty + cu * uaxis[1];
        const nz = ca * this[2] + sa * utz + cu * uaxis[2];

        this[0] = nx; this[1] = ny; this[2] = nz;
        return this;
    }

    rotate(angle, axis) {
        const ca = Math.cos(angle);
        const sa = Math.sin(angle);
        return this.rotateCSA(ca, sa, axis);
    }

    // ================= Utility =================

    dist2(v) {
        const dx = this[0] - v[0];
        const dy = this[1] - v[1];
        const dz = this[2] - v[2];
        return dx * dx + dy * dy + dz * dz;
    }

    toString() {
        return `(${this[0].toFixed(6)}, ${this[1].toFixed(6)}, ${this[2].toFixed(6)})`;
    }
}
