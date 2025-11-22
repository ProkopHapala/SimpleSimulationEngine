
import * as THREE from 'three';
import { vec3 as glVec3 } from 'gl-matrix';
import { Vec3 as Vec3Obj } from '../common_js/Vec3.js';
import { Vec3 as Vec3Arr } from '../common_js/Vec3_arr.js';
import { Vec3 as Vec3Lst } from '../common_js/Vec3_lst.js';



export class SimBase {
    constructor(n) {
        this.n = n;
        this.positions = new Array(n);
        this.velocities = new Array(n);
        this.forces = new Array(n);
    }
    init() { /* Override */ }
    step() { this.run(1, 0.01); }
    run(steps, dt) { /* Override */ }
    sync(attr) { /* Override */ }
}

// --- Vec3Obj ---
export class SimVec3Obj extends SimBase {
    init() {
        for (let i = 0; i < this.n; i++) {
            this.positions[i] = new Vec3Obj((Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100);
            this.velocities[i] = new Vec3Obj((Math.random() - 0.5), (Math.random() - 0.5), (Math.random() - 0.5));
            this.forces[i] = new Vec3Obj();
        }
    }
    sync(attr) {
        for (let i = 0; i < this.n; i++) attr.setXYZ(i, this.positions[i].x, this.positions[i].y, this.positions[i].z);
        attr.needsUpdate = true;
    }
}

export class SimVec3Obj_Manual extends SimVec3Obj {
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) { force[i].x = 0; force[i].y = 0; force[i].z = 0; }
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    const dx = pj.x - pi.x; const dy = pj.y - pi.y; const dz = pj.z - pi.z;
                    const r2 = dx * dx + dy * dy + dz * dz + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    const invR3 = invR * invR * invR;
                    const fx = dx * invR3; const fy = dy * invR3; const fz = dz * invR3;
                    fi.x += fx; fi.y += fy; fi.z += fz;
                    fj.x -= fx; fj.y -= fy; fj.z -= fz;
                }
            }
            for (let i = 0; i < n; i++) {
                vel[i].x += force[i].x * dt; vel[i].y += force[i].y * dt; vel[i].z += force[i].z * dt;
                pos[i].x += vel[i].x * dt; pos[i].y += vel[i].y * dt; pos[i].z += vel[i].z * dt;
            }
        }
    }
}

export class SimVec3Obj_Methods extends SimVec3Obj {
    constructor(n) { super(n); this.diff = new Vec3Obj(); }
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n; const diff = this.diff;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) force[i].set(0, 0, 0);
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    diff.setSub(pj, pi);
                    const r2 = diff.norm2() + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    diff.mulScalar(invR * invR * invR);
                    fi.add(diff); fj.sub(diff);
                }
            }
            for (let i = 0; i < n; i++) { vel[i].addMul(force[i], dt); pos[i].addMul(vel[i], dt); }
        }
    }
}

// --- Vec3Arr ---
export class SimVec3Arr extends SimBase {
    init() {
        for (let i = 0; i < this.n; i++) {
            this.positions[i] = new Vec3Arr((Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100);
            this.velocities[i] = new Vec3Arr((Math.random() - 0.5), (Math.random() - 0.5), (Math.random() - 0.5));
            this.forces[i] = new Vec3Arr();
        }
    }
    sync(attr) {
        for (let i = 0; i < this.n; i++) attr.setXYZ(i, this.positions[i][0], this.positions[i][1], this.positions[i][2]);
        attr.needsUpdate = true;
    }
}

export class SimVec3Arr_Manual extends SimVec3Arr {
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) { force[i][0] = 0; force[i][1] = 0; force[i][2] = 0; }
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    const dx = pj[0] - pi[0]; const dy = pj[1] - pi[1]; const dz = pj[2] - pi[2];
                    const r2 = dx * dx + dy * dy + dz * dz + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    const invR3 = invR * invR * invR;
                    const fx = dx * invR3; const fy = dy * invR3; const fz = dz * invR3;
                    fi[0] += fx; fi[1] += fy; fi[2] += fz;
                    fj[0] -= fx; fj[1] -= fy; fj[2] -= fz;
                }
            }
            for (let i = 0; i < n; i++) {
                vel[i][0] += force[i][0] * dt; vel[i][1] += force[i][1] * dt; vel[i][2] += force[i][2] * dt;
                pos[i][0] += vel[i][0] * dt; pos[i][1] += vel[i][1] * dt; pos[i][2] += vel[i][2] * dt;
            }
        }
    }
}

export class SimVec3Arr_Methods extends SimVec3Arr {
    constructor(n) { super(n); this.diff = new Vec3Arr(); }
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n; const diff = this.diff;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) force[i].set(0, 0, 0);
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    diff.setSub(pj, pi);
                    const r2 = diff.norm2() + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    diff.mulScalar(invR * invR * invR);
                    fi.add(diff); fj.sub(diff);
                }
            }
            for (let i = 0; i < n; i++) { vel[i].addMul(force[i], dt); pos[i].addMul(vel[i], dt); }
        }
    }
}

// --- Vec3Lst ---
export class SimVec3Lst extends SimBase {
    init() {
        for (let i = 0; i < this.n; i++) {
            this.positions[i] = new Vec3Lst((Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100);
            this.velocities[i] = new Vec3Lst((Math.random() - 0.5), (Math.random() - 0.5), (Math.random() - 0.5));
            this.forces[i] = new Vec3Lst();
        }
    }
    sync(attr) {
        for (let i = 0; i < this.n; i++) attr.setXYZ(i, this.positions[i][0], this.positions[i][1], this.positions[i][2]);
        attr.needsUpdate = true;
    }
}
export class SimVec3Lst_Manual extends SimVec3Lst {
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) { force[i][0] = 0; force[i][1] = 0; force[i][2] = 0; }
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    const dx = pj[0] - pi[0]; const dy = pj[1] - pi[1]; const dz = pj[2] - pi[2];
                    const r2 = dx * dx + dy * dy + dz * dz + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    const invR3 = invR * invR * invR;
                    const fx = dx * invR3; const fy = dy * invR3; const fz = dz * invR3;
                    fi[0] += fx; fi[1] += fy; fi[2] += fz;
                    fj[0] -= fx; fj[1] -= fy; fj[2] -= fz;
                }
            }
            for (let i = 0; i < n; i++) {
                vel[i][0] += force[i][0] * dt; vel[i][1] += force[i][1] * dt; vel[i][2] += force[i][2] * dt;
                pos[i][0] += vel[i][0] * dt; pos[i][1] += vel[i][1] * dt; pos[i][2] += vel[i][2] * dt;
            }
        }
    }
}
export class SimVec3Lst_Methods extends SimVec3Lst {
    constructor(n) { super(n); this.diff = new Vec3Lst(); }
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n; const diff = this.diff;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) force[i].set(0, 0, 0);
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    diff.setSub(pj, pi);
                    const r2 = diff.norm2() + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    diff.mulScalar(invR * invR * invR);
                    fi.add(diff); fj.sub(diff);
                }
            }
            for (let i = 0; i < n; i++) { vel[i].addMul(force[i], dt); pos[i].addMul(vel[i], dt); }
        }
    }
}

// --- THREE ---
export class SimTHREE extends SimBase {
    init() {
        for (let i = 0; i < this.n; i++) {
            this.positions[i] = new THREE.Vector3((Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100);
            this.velocities[i] = new THREE.Vector3((Math.random() - 0.5), (Math.random() - 0.5), (Math.random() - 0.5));
            this.forces[i] = new THREE.Vector3();
        }
    }
    sync(attr) {
        for (let i = 0; i < this.n; i++) attr.setXYZ(i, this.positions[i].x, this.positions[i].y, this.positions[i].z);
        attr.needsUpdate = true;
    }
}
export class SimTHREE_Manual extends SimTHREE {
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) { force[i].x = 0; force[i].y = 0; force[i].z = 0; }
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    const dx = pj.x - pi.x; const dy = pj.y - pi.y; const dz = pj.z - pi.z;
                    const r2 = dx * dx + dy * dy + dz * dz + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    const invR3 = invR * invR * invR;
                    const fx = dx * invR3; const fy = dy * invR3; const fz = dz * invR3;
                    fi.x += fx; fi.y += fy; fi.z += fz;
                    fj.x -= fx; fj.y -= fy; fj.z -= fz;
                }
            }
            for (let i = 0; i < n; i++) {
                vel[i].x += force[i].x * dt; vel[i].y += force[i].y * dt; vel[i].z += force[i].z * dt;
                pos[i].x += vel[i].x * dt; pos[i].y += vel[i].y * dt; pos[i].z += vel[i].z * dt;
            }
        }
    }
}
export class SimTHREE_Methods extends SimTHREE {
    constructor(n) { super(n); this.diff = new THREE.Vector3(); }
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n; const diff = this.diff;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) force[i].set(0, 0, 0);
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    diff.subVectors(pj, pi);
                    const r2 = diff.lengthSq() + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    diff.multiplyScalar(invR * invR * invR);
                    fi.add(diff); fj.sub(diff);
                }
            }
            for (let i = 0; i < n; i++) { vel[i].addScaledVector(force[i], dt); pos[i].addScaledVector(vel[i], dt); }
        }
    }
}

// --- GLM ---
export class SimGLM extends SimBase {
    init() {
        for (let i = 0; i < this.n; i++) {
            this.positions[i] = glVec3.fromValues((Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100, (Math.random() - 0.5) * 100);
            this.velocities[i] = glVec3.fromValues((Math.random() - 0.5), (Math.random() - 0.5), (Math.random() - 0.5));
            this.forces[i] = glVec3.create();
        }
    }
    sync(attr) {
        for (let i = 0; i < this.n; i++) attr.setXYZ(i, this.positions[i][0], this.positions[i][1], this.positions[i][2]);
        attr.needsUpdate = true;
    }
}
export class SimGLM_Manual extends SimGLM {
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) { force[i][0] = 0; force[i][1] = 0; force[i][2] = 0; }
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    const dx = pj[0] - pi[0]; const dy = pj[1] - pi[1]; const dz = pj[2] - pi[2];
                    const r2 = dx * dx + dy * dy + dz * dz + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    const invR3 = invR * invR * invR;
                    const fx = dx * invR3; const fy = dy * invR3; const fz = dz * invR3;
                    fi[0] += fx; fi[1] += fy; fi[2] += fz;
                    fj[0] -= fx; fj[1] -= fy; fj[2] -= fz;
                }
            }
            for (let i = 0; i < n; i++) {
                vel[i][0] += force[i][0] * dt; vel[i][1] += force[i][1] * dt; vel[i][2] += force[i][2] * dt;
                pos[i][0] += vel[i][0] * dt; pos[i][1] += vel[i][1] * dt; pos[i][2] += vel[i][2] * dt;
            }
        }
    }
}
export class SimGLM_Methods extends SimGLM {
    constructor(n) { super(n); this.diff = glVec3.create(); }
    run(steps, dt) {
        const pos = this.positions; const vel = this.velocities; const force = this.forces; const n = this.n; const diff = this.diff;
        for (let s = 0; s < steps; s++) {
            for (let i = 0; i < n; i++) glVec3.set(force[i], 0, 0, 0);
            for (let i = 0; i < n; i++) {
                const pi = pos[i]; const fi = force[i];
                for (let j = i + 1; j < n; j++) {
                    const pj = pos[j]; const fj = force[j];
                    glVec3.subtract(diff, pj, pi);
                    const r2 = glVec3.sqrLen(diff) + 1e-4;
                    const invR = 1.0 / Math.sqrt(r2);
                    glVec3.scale(diff, diff, invR * invR * invR);
                    glVec3.add(fi, fi, diff); glVec3.subtract(fj, fj, diff);
                }
            }
            for (let i = 0; i < n; i++) { glVec3.scaleAndAdd(vel[i], vel[i], force[i], dt); glVec3.scaleAndAdd(pos[i], pos[i], vel[i], dt); }
        }
    }
}
