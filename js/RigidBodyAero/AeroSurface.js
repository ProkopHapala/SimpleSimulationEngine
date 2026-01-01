import * as THREE from 'three';

const defaultParams = {
    area: 1.0,
    CD0: 0.02,
    dCD: 0.9,
    dCDS: 0.9,
    dCL: 6.28,
    dCLS: 2.82,
    sStall: 0.16,
    wStall: 0.08
};

export function geomFromString(text, z_shift = 0) {
    const lines = text.split('\n').filter(l => l.trim() && !l.trim().startsWith('#'));
    const geom = [];
    for (const line of lines) {
        const g = line.trim().split(/\s+/).map(Number);
        if (g.length >= 9 && g.every(v => Number.isFinite(v))) {
            geom.push({
                pos: new THREE.Vector3(g[0], g[1], g[2] + z_shift),
                t: new THREE.Vector3(g[3], g[4], g[5]).normalize(),
                n: new THREE.Vector3(g[6], g[7], g[8]).normalize()
            });
        }
    }
    return geom;
}

export function paramsFromString(text) {
    const lines = text.split('\n').filter(l => l.trim() && !l.trim().startsWith('#'));
    const params = [];
    for (const line of lines) {
        const a = line.trim().split(/\s+/).map(Number);
        if (a.length >= 8 && a.every(v => Number.isFinite(v))) {
            params.push({
                area: a[0], CD0: a[1], dCD: a[2], dCDS: a[3],
                dCL: a[4], dCLS: a[5], sStall: a[6], wStall: a[7]
            });
        }
    }
    return params;
}

export function buildPoints(geom, params) {
    const p = [];
    const fallback = params[0] || { ...defaultParams };
    for (let i = 0; i < geom.length; i++) {
        const g = geom[i];
        const par = params[i] || fallback;
        if (g) {
            p.push({
                pos: g.pos.clone(),
                t: g.t.clone(),
                n: g.n.clone(),
                params: { ...par }
            });
        }
    }
    return p;
}

export function computeCoefficients(params, aoaRad) {
    const ca = Math.cos(aoaRad);
    const sa = Math.sin(aoaRad);
    const abs_sa = Math.abs(sa);
    let wS = 0;
    if (abs_sa > params.sStall) {
        const x = (abs_sa - params.sStall) / params.wStall;
        wS = (x > 1) ? 1 : (3 * x * x - 2 * x * x * x);
    }
    const mS = 1 - wS;
    let CL_ca = ca;
    let CL_sa = sa;
    if (CL_ca < 0) { CL_ca = -CL_ca; CL_sa = -CL_sa; }
    const CL = (mS * params.dCL + wS * params.dCLS * CL_ca) * CL_sa;
    const CD = params.CD0 + (mS * params.dCD * abs_sa + wS * params.dCDS) * abs_sa;
    return { CL, CD };
}

export function computePointForce(point, body, windVec = new THREE.Vector3(), R = null) {
    if (!point || !body) return { force: new THREE.Vector3(), position: new THREE.Vector3() };
    const rot = R || new THREE.Matrix3().setFromMatrix4(new THREE.Matrix4().makeRotationFromQuaternion(body.quat));
    const r_world = point.pos.clone().applyMatrix3(rot);
    const p_world = body.pos.clone().add(r_world);

    const v_point = body.vel.clone().add(new THREE.Vector3().crossVectors(body.angVel, r_world));
    const v_air = windVec.clone().sub(v_point);
    const v2 = v_air.lengthSq();
    if (v2 < 1e-6) {
        return { force: new THREE.Vector3(), position: p_world };
    }

    const u_air = v_air.clone().normalize();
    const t_world = point.t.clone().applyMatrix3(rot);
    const n_world = point.n.clone().applyMatrix3(rot);
    const ca = u_air.dot(n_world);
    const sa = u_air.dot(t_world);
    const aoa = Math.atan2(sa, ca);
    const { CL, CD } = computeCoefficients(point.params, aoa);
    const pref = v2 * point.params.area;
    const liftDir = t_world.clone().addScaledVector(u_air, -sa).normalize();
    const f_total = liftDir.multiplyScalar(CL * pref).add(u_air.multiplyScalar(CD * pref));
    return { force: f_total, position: p_world };
}
