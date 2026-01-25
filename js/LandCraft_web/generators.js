// generators.js - Terrain generation algorithms ported from Shadertoy
import { NOISE_LIB, NOISE_FLAVORS } from './noise-lib.js';

const toF32 = (value) => {
    const num = Number(value);
    if (!Number.isFinite(num)) return '0.0';
    const str = num.toString();
    return str.includes('.') ? str : `${str}.0`;
};

const buildAnalyticSinNoiseFn = (scale, period, fnName = 'analyticSinNoise') => {
    const freqScale = Math.max(0.001, Number(scale ?? 1));
    const periodVal = Math.max(Number(period ?? 256), 1e-6);
    const base = freqScale * 6.28318530718 / periodVal;
    const baseLiteral = toF32(base);
    return {
        body: `fn ${fnName}(p: vec2f) -> f32 {
    let base = ${baseLiteral};
    var s = 0.0;
    s += sin(dot(p, vec2f(1.0, 0.0)) * base + 1.234);
    s += sin(dot(p, vec2f(0.0, 1.0)) * base * 1.61803399 + 2.345);
    s += sin(dot(p, normalize(vec2f(1.0, 1.0))) * base * 1.41421356 + 3.456);
    s += sin(dot(p, normalize(vec2f(1.0, -2.0))) * base * 1.90211303 + 4.321);
    s += sin(dot(p, normalize(vec2f(-2.0, 1.0))) * base * 2.23606798 + 5.678);
    return 0.5 + 0.1 * s;
}`,
        fn: fnName,
        octaveStep: 'x = mat2x2f(1.7, -0.9, 0.9, 1.7) * x;'
    };
};

const buildNoiseFlavorSnippet = (flavorKey, params) => {
    const key = flavorKey ?? 'sin';
    if (key === 'sin') {
        return buildAnalyticSinNoiseFn(params?.scale, params?.period, 'sinFlavorNoise');
    }
    const fallback = NOISE_FLAVORS[key] ?? NOISE_FLAVORS.sin;
    return { body: fallback.body, fn: fallback.fn, octaveStep: 'x *= 2.0;' };
};

const buildOctaveMatrixSnippet = (p, varName, preMatrix = null) => {
    const angle = toF32(p.rotation ?? 0);
    const grow = toF32(p.grow ?? 2);
    const anis = toF32(p.anisotropy ?? 1);
    const skew = toF32(p.skew ?? 0);
    const base = preMatrix ? `${preMatrix} * userMat` : 'userMat';
    return `let angle = ${angle} * 0.017453292519943295;
    let c = cos(angle);
    let s = sin(angle);
    let grow = ${grow};
    let anis = ${anis};
    let skew = ${skew};
    let stretchX = grow;
    let stretchY = grow * anis;
    let rot = mat2x2f(c, -s, s, c);
    let shear = mat2x2f(1.0, skew, 0.0, 1.0);
    let stretch = mat2x2f(stretchX, 0.0, 0.0, stretchY);
    let userMat = rot * shear * stretch;
    let ${varName} = ${base};`;
};

export const GENERATORS = {
    // Raw noise generators for debugging (single octave, no FBM)
    sinNoise: { id: 'sinNoise', label: 'Sin Noise (Debug)', params: [{name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4}, {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256}], wgsl: (p) => {
        const analytic = buildAnalyticSinNoiseFn(p.scale, p.period, 'sinDebugNoise');
        return `${analytic.body}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let pSample = uv * ${p.scale};
    let h = ${analytic.fn}(pSample);
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
    } },
    valueNoise: { id: 'valueNoise', label: 'Value Noise (Debug)', params: [{name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4}, {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256}], wgsl: (p) => `${NOISE_LIB.valueNoise}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let p = uv * ${p.scale};
    let h = valueNoise(p) * 0.5 + 0.5;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}` },
    simplexNoise: { id: 'simplexNoise', label: 'Simplex Noise (Debug)', params: [{name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4}, {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256}], wgsl: (p) => `${NOISE_LIB.simplexNoise}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let p = uv * ${p.scale};
    let h = simplexNoise(p) * 0.5 + 0.5;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}` },
    gradNoise: { id: 'gradNoise', label: 'Gradient Noise (Debug)', params: [{name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4}, {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256}], wgsl: (p) => `${NOISE_LIB.gradNoise}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let p = uv * ${p.scale};
    let h = gradNoise(p) * 0.5 + 0.5;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}` },
    orbitNoise: { id: 'orbitNoise', label: 'Orbit Noise (Debug)', params: [{name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4}, {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256}], wgsl: (p) => `${NOISE_LIB.orbitNoise}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let p = uv * ${p.scale};
    let h = orbitNoise(p);
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}` },
    // FBM generator with injectable periodic noise
    fbm: {
        id: 'fbm',
        label: 'FBM',
        noiseFlavor: 'sin',
        params: [
            { name: 'octaves', label: 'Octaves', type: 'int', min: 1, max: 8, default: 6 },
            { name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4 },
            { name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256 },
            { name: 'grow', label: 'Grow', type: 'float', min: 0.5, max: 3.0, default: 1.8 },
            { name: 'anisotropy', label: 'Anisotropy', type: 'float', min: 0.2, max: 2.0, default: 1.0 },
            { name: 'rotation', label: 'Rotation', type: 'float', min: 0, max: 360, default: 35 },
            { name: 'skew', label: 'Skew', type: 'float', min: -1.0, max: 1.0, default: 0.0 }
        ],
        wgsl: (p, flavorKey) => {
            const flavor = buildNoiseFlavorSnippet(flavorKey, p);
            const coordExpr = `(uv * ${p.scale})`;
            const matrixSnippet = buildOctaveMatrixSnippet(p, 'octMat');
            return `${flavor.body}
fn noise(p: vec2f) -> f32 { return ${flavor.fn}(p); }
fn fbm(p: vec2f, octaves: i32) -> f32 {
    var v = 0.0; var a = 0.5; var x = p;
${matrixSnippet}
    for(var i = 0; i < octaves; i++){
        v += noise(x) * a;
        x = octMat * x;
        a *= 0.5;
    }
    return v;
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let p = ${coordExpr};
    var h = fbm(p, ${p.octaves});
    h = h * 1.5 - 0.2;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
        }
    },
    elevated: { id: 'elevated', label: 'Elevated Terrain', noiseFlavor: 'value', params: [
        {name: 'octaves', label: 'Octaves', type: 'int', min: 1, max: 16, default: 16},
        {name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4},
        {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256},
        { name: 'grow', label: 'Grow', type: 'float', min: 0.5, max: 3.0, default: 1.7 },
        { name: 'anisotropy', label: 'Anisotropy', type: 'float', min: 0.2, max: 2.0, default: 0.9 },
        { name: 'rotation', label: 'Rotation', type: 'float', min: 0, max: 360, default: 25 },
        { name: 'skew', label: 'Skew', type: 'float', min: -1.0, max: 1.0, default: 0.2 }
    ], wgsl: (p, flavorKey) => {
        const flavor = buildNoiseFlavorSnippet(flavorKey ?? 'value', p);
        const matrixSnippet = `let baseMat = mat2x2f(0.8, -0.6, 0.6, 0.8);
${buildOctaveMatrixSnippet(p, 'octMat', 'baseMat')}`;
        return `${flavor.body}
fn noise(p: vec2f) -> f32 { return ${flavor.fn}(p); }
fn noised(p: vec2f) -> vec3f {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    let du = 6.0*f*(1.0-f);
    let a = noise(i); let b = noise(i+vec2f(1,0)); let c = noise(i+vec2f(0,1)); let d = noise(i+1.0);
    return vec3f(a+(b-a)*u.x+(c-a)*u.y+(a-b-c+d)*u.x*u.y, du*(vec2f(b-a,c-a)+(a-b-c+d)*u.yx));
}
fn terrainH(x: vec2f, octaves: i32) -> f32 {
    var p = x * 0.003 / 250.0;
    var a = 0.0; var b = 1.0; var d = vec2f(0.0);
${matrixSnippet}
    for(var i=0; i<octaves; i++) {
        let n = noised(p);
        d += n.yz;
        a += b * n.x / (1.0 + dot(d, d));
        b *= 0.5;
        p = octMat * p;
    }
    return 250.0 * 120.0 * a;
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let centered = (vec2f(id.xy) - vec2f(globals.resX, globals.resY) * 0.5) * ${p.scale};
    let raw = terrainH(centered, ${p.octaves});
    let h = raw * 0.00002 + 0.35;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
    } },
    sirenian: {
        id: 'sirenian',
        label: 'Sirenian Dawn',
        noiseFlavor: 'value',
        params: [
            { name: 'octaves', label: 'Octaves', type: 'int', min: 1, max: 5, default: 5 },
            { name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4 },
            { name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256 },
            { name: 'grow', label: 'Grow', type: 'float', min: 0.5, max: 3.0, default: 1.6 },
            { name: 'anisotropy', label: 'Anisotropy', type: 'float', min: 0.2, max: 2.0, default: 1.1 },
            { name: 'rotation', label: 'Rotation', type: 'float', min: 0, max: 360, default: 15 },
            { name: 'skew', label: 'Skew', type: 'float', min: -1.0, max: 1.0, default: 0.0 }
        ],
        wgsl: (p, flavorKey) => {
            const flavor = buildNoiseFlavorSnippet(flavorKey ?? 'value', p);
            const matrixSnippet = `let baseMat = mat2x2f(0.8, 0.6, -0.6, 0.8);
${buildOctaveMatrixSnippet(p, 'octMat', 'baseMat')}`;
            return `${flavor.body}
fn noise(p: vec2f) -> f32 { return ${flavor.fn}(p); }
fn noised(p: vec2f) -> vec3f {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    let du = 6.0*f*(1.0-f);
    let a = noise(i); let b = noise(i+vec2f(1,0)); let c = noise(i+vec2f(0,1)); let d = noise(i+1.0);
    return vec3f(a+(b-a)*u.x+(c-a)*u.y+(a-b-c+d)*u.x*u.y, du*(vec2f(b-a,c-a)+(a-b-c+d)*u.yx));
}
fn terrain(pIn: vec2f, octaves: i32) -> f32 {
    var p = pIn;
    var rz = 0.0; var z = 1.0; var d = vec2f(0.0); var zscl = -0.4; var zz = 5.0;
${matrixSnippet}
    for(var i=0; i<octaves; i++) {
        let n = noised(p);
        d += vec2f(abs(n.y), abs(n.z));
        d -= vec2f(smoothstep(-0.5, 1.5, n.y), smoothstep(-0.5, 1.5, n.z));
        zz -= 1.0;
        rz += z * n.x / (dot(d, d) + 0.85);
        z *= zscl;
        zscl *= 0.8;
        p = octMat * p;
    }
    return rz / (smoothstep(1.5, -0.5, rz) + 0.75);
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let pos = uv * ${p.scale} * 6.0;
    var h = terrain(pos, ${p.octaves});
    let dist = distance(uv, vec2f(0.5));
    h -= smoothstep(0.3, 0.5, dist);
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
        }
    }
    ,
    elevated_ref: {
        id: 'elevated_ref',
        label: 'Elevated (Ref)',
        params: [
            {name: 'octaves', label: 'Octaves', type: 'int', min: 1, max: 16, default: 16},
            {name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4},
            {name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256},
            { name: 'grow', label: 'Grow', type: 'float', min: 0.5, max: 3.0, default: 2.0 },
            { name: 'anisotropy', label: 'Anisotropy', type: 'float', min: 0.2, max: 2.0, default: 1.0 },
            { name: 'rotation', label: 'Rotation', type: 'float', min: 0, max: 360, default: 0 },
            { name: 'skew', label: 'Skew', type: 'float', min: -1.0, max: 1.0, default: 0.0 },
            { name: 'smoothNoise', label: 'Smooth Noise', type: 'int', min: 0, max: 1, default: 0 }
        ],
        extraInputs: ({ noise256 }) => [noise256],
        wgsl: (p) => {
            const matrixSnippet = `let baseMat = mat2x2f(0.8, -0.6, 0.6, 0.8);
${buildOctaveMatrixSnippet({ ...p, rotation: p.rotation ?? 0, grow: p.grow ?? 2.0, anisotropy: p.anisotropy ?? 1.0, skew: p.skew ?? 0.0 }, 'octMat', 'baseMat')}`;
            const useSmooth = Number(p.smoothNoise ?? 0) ? 1 : 0;
            return `fn noised(x: vec2f) -> vec3f {
    let f = fract(x);
    var u: vec2f;
    var du: vec2f;
    if (${useSmooth} == 0) {
        u = f*f*(3.0-2.0*f);
        du = 6.0*f*(1.0-f);
    } else {
        u = f*f*f*(f*(f*6.0-15.0)+10.0);
        du = 30.0*f*f*(f*(f-2.0)+1.0);
    }
    let pi = vec2<i32>(floor(x));
    let mask256 = vec2<i32>(255, 255);
    let a = textureLoad(in_0, (pi + vec2<i32>(0,0)) & mask256, 0).x;
    let b = textureLoad(in_0, (pi + vec2<i32>(1,0)) & mask256, 0).x;
    let c = textureLoad(in_0, (pi + vec2<i32>(0,1)) & mask256, 0).x;
    let d = textureLoad(in_0, (pi + vec2<i32>(1,1)) & mask256, 0).x;
    return vec3f(a+(b-a)*u.x+(c-a)*u.y+(a-b-c+d)*u.x*u.y, du*(vec2f(b-a,c-a)+(a-b-c+d)*u.yx));
}
fn terrainH(x: vec2f, octaves: i32) -> f32 {
    var p = x*0.003/250.0;
    var a = 0.0; var b = 1.0; var d = vec2f(0.0);
${matrixSnippet}
    for (var i=0; i<octaves; i++) {
        let n = noised(p);
        d += n.yz;
        a += b*n.x/(1.0 + dot(d,d));
        b *= 0.5;
        p = octMat * p;
    }
    if (${useSmooth} == 1) { a *= 0.9; }
    return 250.0*120.0*a;
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let centered = (vec2f(id.xy) - vec2f(globals.resX, globals.resY)*0.5) * ${p.scale};
    let raw = terrainH(centered, ${p.octaves});
    let h = raw * 0.00002 + 0.35;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
        }
    },
    sirenian_ref: {
        id: 'sirenian_ref',
        label: 'Sirenian Dawn (Ref)',
        params: [
            { name: 'octaves', label: 'Octaves', type: 'int', min: 1, max: 8, default: 5 },
            { name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4 },
            { name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256 },
            { name: 'grow', label: 'Grow', type: 'float', min: 0.5, max: 4.0, default: 2.95 },
            { name: 'anisotropy', label: 'Anisotropy', type: 'float', min: 0.2, max: 2.0, default: 1.0 },
            { name: 'rotation', label: 'Rotation', type: 'float', min: 0, max: 360, default: 0 },
            { name: 'skew', label: 'Skew', type: 'float', min: -1.0, max: 1.0, default: 0.0 }
        ],
        extraInputs: ({ noise256 }) => [noise256],
        wgsl: (p) => {
            const matrixSnippet = `let baseMat = mat2x2f(0.80, 0.60, -0.60, 0.80);
${buildOctaveMatrixSnippet({ ...p, rotation: p.rotation ?? 0, grow: p.grow ?? 2.95, anisotropy: p.anisotropy ?? 1.0, skew: p.skew ?? 0.0 }, 'octMat', 'baseMat')}`;
            return `fn noised(x: vec2f) -> vec3f {
    let p = floor(x);
    let f = fract(x);
    let u = f*f*(3.0-2.0*f);
    let ip = vec2<i32>(p);
    let mask256 = vec2<i32>(255, 255);
    let a = textureLoad(in_0, (ip + vec2<i32>(0,0)) & mask256, 0).x;
    let b = textureLoad(in_0, (ip + vec2<i32>(1,0)) & mask256, 0).x;
    let c = textureLoad(in_0, (ip + vec2<i32>(0,1)) & mask256, 0).x;
    let d = textureLoad(in_0, (ip + vec2<i32>(1,1)) & mask256, 0).x;
    let du = 6.0*f*(1.0-f);
    return vec3f(a+(b-a)*u.x+(c-a)*u.y+(a-b-c+d)*u.x*u.y, du*(vec2f(b-a,c-a)+(a-b-c+d)*u.yx));
}
fn terrain(pIn: vec2f, octaves: i32) -> f32 {
    var p = pIn;
    var rz = 0.0; var z = 1.0; var d = vec2f(0.0); var zscl = -0.4; var zz = 5.0;
${matrixSnippet}
    for (var i=0; i<octaves; i++) {
        let n = noised(p);
        d += pow(abs(n.yz), vec2f(zz));
        d -= vec2f(smoothstep(-0.5, 1.5, n.y), smoothstep(-0.5, 1.5, n.z));
        zz -= 1.0;
        rz += z*n.x/(dot(d,d)+0.85);
        z *= zscl;
        zscl *= 0.8;
        p = octMat * p;
    }
    return rz/(smoothstep(1.5, -0.5, rz)+0.75);
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    let pos = uv * ${p.scale} * 6.0;
    var h = terrain(pos, ${p.octaves});
    let dist = distance(uv, vec2f(0.5));
    h -= smoothstep(0.3, 0.5, dist);
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
        }
    },
    canyo_ref: {
        id: 'canyo_ref',
        label: 'Canyo (Ref)',
        params: [
            { name: 'scale', label: 'Scale', type: 'float', min: 1, max: 100, default: 4 },
            { name: 'period', label: 'Period', type: 'float', min: 64, max: 1024, default: 256 }
        ],
        extraInputs: ({ noise1024 }) => [noise1024],
        wgsl: (p) => {
            return `fn textureGood(uvIn: vec2f) -> vec4f {
    let uv = uvIn*1024.0 - 0.5;
    let iuv = floor(uv);
    let f = fract(uv);
    let p1 = vec2<i32>(i32(iuv.x+0.5), i32(iuv.y+0.5));
    let p2 = vec2<i32>(i32(iuv.x+1.5), i32(iuv.y+0.5));
    let p3 = vec2<i32>(i32(iuv.x+0.5), i32(iuv.y+1.5));
    let p4 = vec2<i32>(i32(iuv.x+1.5), i32(iuv.y+1.5));
    let mask1024 = vec2<i32>(1023, 1023);
    let rg1 = textureLoad(in_0, p1 & mask1024, 0);
    let rg2 = textureLoad(in_0, p2 & mask1024, 0);
    let rg3 = textureLoad(in_0, p3 & mask1024, 0);
    let rg4 = textureLoad(in_0, p4 & mask1024, 0);
    return mix(mix(rg1, rg2, f.x), mix(rg3, rg4, f.x), f.y);
}
fn terrain(q: vec2f) -> f32 {
    let th = smoothstep(0.0, 0.7, textureGood(0.001*q).x);
    let rr = smoothstep(0.1, 0.5, textureGood(2.0*0.03*q).y);
    var h = 1.9;
    h += th*7.0;
    h += 0.3*rr;
    return -h;
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let centered = (vec2f(id.xy) - vec2f(globals.resX, globals.resY)*0.5) * ${p.scale};
    let h = terrain(centered);
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
        }
    },

    procErosionHeightmap_ref: {
        id: 'procErosionHeightmap_ref',
        label: 'Procedural Erosion Heightmap (Ref)',
        params: [
            { name: 'samples', label: 'Total Samples', type: 'int', min: 1, max: 64, default: 16 },
            { name: 'range', label: 'Range', type: 'int', min: 4, max: 32, default: 16 },
            { name: 'thresholdMin', label: 'Thr Min', type: 'float', min: 0.0, max: 1.0, default: 0.0 },
            { name: 'thresholdMax', label: 'Thr Max', type: 'float', min: 0.0, max: 1.0, default: 1.0 },
            { name: 'domainScale', label: 'Domain Scale', type: 'float', min: 0.5, max: 20.0, default: 4.0 },
            { name: 'samplesPerFrame', label: 'Samples/Frame', type: 'int', min: 1, max: 16, default: 4 },
            { name: 'progressive', label: 'Progressive', type: 'bool', default: false }
        ],
        extraInputs: ({ noise1024 }) => [noise1024],
        progressive: true,
        wgsl: (p) => {
            const samples = Math.max(1, Math.min(64, Number(p.samples ?? 16) | 0));
            const samplesPerFrame = Math.max(1, Math.min(16, Number(p.samplesPerFrame ?? 4) | 0));
            const range = Math.max(1, Math.min(64, Number(p.range ?? 16) | 0));
            const thrMin = toF32(p.thresholdMin ?? 0.0);
            const thrMax = toF32(p.thresholdMax ?? 1.0);
            const domainScale = toF32(p.domainScale ?? 4.0);
            const progressive = !!p.progressive;
            const common = `fn hash11(x: f32) -> f32 { return fract(sin(x) * 43758.5453123); }
fn hash31(x: f32) -> vec3f { return vec3f(hash11(x + 0.1), hash11(x + 1.7), hash11(x + 2.9)); }

fn perlin(uv: vec2f) -> f32 {
    let dim = vec2f(textureDimensions(in_0));
    var base = (uv * ${domainScale} + vec2f(8.2813, 1.42114)) / 4.0;
    var occ = vec2f(0.0);
    var a = 1.0;
    for (var i = 0; i < 7; i++) {
        let uvTex = base * dim - vec2f(0.5);
        let iuv = floor(uvTex);
        let f = fract(uvTex);
        let p1 = vec2<i32>(i32(iuv.x + 0.5), i32(iuv.y + 0.5));
        let p2 = vec2<i32>(i32(iuv.x + 1.5), i32(iuv.y + 0.5));
        let p3 = vec2<i32>(i32(iuv.x + 0.5), i32(iuv.y + 1.5));
        let p4 = vec2<i32>(i32(iuv.x + 1.5), i32(iuv.y + 1.5));
        let mask = vec2<i32>(i32(dim.x) - 1, i32(dim.y) - 1);
        let v1 = textureLoad(in_0, p1 & mask, 0).x;
        let v2 = textureLoad(in_0, p2 & mask, 0).x;
        let v3 = textureLoad(in_0, p3 & mask, 0).x;
        let v4 = textureLoad(in_0, p4 & mask, 0).x;
        let vx0 = mix(v1, v2, f.x);
        let vx1 = mix(v3, v4, f.x);
        let v = mix(vx0, vx1, f.y);
        occ += vec2f(v, 1.0) * a;
        base *= 0.5;
        a *= 2.0;
    }
    var v = occ.x / occ.y;
    v = v * 2.0 - 1.0;
    v = tanh(v * 2.0);
    v = v * 0.5 + 0.5;
    v *= v;
    return v;
}

fn occAt(coord: vec2<i32>) -> f32 {
    let uv = vec2f(coord) / vec2f(f32(textureDimensions(in_0).x), f32(textureDimensions(in_0).y));
    return perlin(uv);
}

fn isIn(coord: vec2<i32>, threshold: f32) -> bool { return occAt(coord) > threshold; }

fn sampleDistanceToEdge(center: vec2f, threshold: f32, iRange: i32) -> f32 {
    let start = vec2<i32>(floor(center));
    let fragIsIn = isIn(start, threshold);
    let rangeF = f32(iRange);
    let maxSqr = rangeF * rangeF;
    var best = maxSqr;

    for (var dx = -iRange; dx <= iRange; dx++) {
        for (var dy = -iRange; dy <= iRange; dy++) {
            let delta = vec2f(f32(dx), f32(dy));
            let d2 = dot(delta, delta);
            if (d2 >= maxSqr) { continue; }
            if (d2 >= best) { continue; }
            let scan = start + vec2<i32>(dx, dy);
            let scanIsIn = isIn(scan, threshold);
            if (scanIsIn != fragIsIn) { best = d2; }
        }
    }

    var dist = sqrt(best);
    dist -= 0.5;
    if (fragIsIn) { dist = -dist; }
    dist /= rangeF * 2.0;
    dist = 0.5 - dist;
    dist = smoothstep(0.0, 1.0, dist);
    return dist;
}
`;            if (progressive) {
                return common + `

@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let frag = vec2f(id.xy);
    var acc = 0.0;
    let baseIdx: u32 = globals.frame * u32(${samplesPerFrame});
    for (var i: u32 = 0u; i < u32(${samplesPerFrame}); i++) {
        let sampleIdx: u32 = baseIdx + i;
        if (sampleIdx >= u32(${samples})) { break; }
        let seed = f32(sampleIdx) * 1.618033 + f32(i) * 0.7548776662;
        var n = hash31(seed);
        let jitter = n.xy - vec2f(0.5);
        let thr = mix(${thrMin}, ${thrMax}, n.z);
        let center = frag + jitter;
        acc += sampleDistanceToEdge(center, thr, ${range});
    }
    let h = acc / f32(${samplesPerFrame});
    let prevH = textureLoad(in_1, vec2<i32>(id.xy), 0).x;
    let totalFrames = (${samples} + ${samplesPerFrame} - 1) / ${samplesPerFrame};
    let currentFrame = globals.frame + 1u;
    let hAvg = (prevH * f32(currentFrame - 1u) + h) / f32(currentFrame);
    textureStore(out_1, vec2<i32>(id.xy), vec4f(hAvg, 0.0, 0.0, 1.0));
    textureStore(out_0, vec2<i32>(id.xy), vec4f(hAvg, 0.0, 0.0, 1.0));
}`;
            }
            // Non-progressive single-pass
            return common + `

@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let frag = vec2f(id.xy);
    var acc = 0.0;
    for (var i = 0; i < ${samples}; i++) {
        let seed = f32(i) * 1.618033 + f32(i) * 0.7548776662;
        var n = hash31(seed);
        let jitter = n.xy - vec2f(0.5);
        let thr = mix(${thrMin}, ${thrMax}, n.z);
        let center = frag + jitter;
        acc += sampleDistanceToEdge(center, thr, ${range});
    }
    let h = acc / f32(${samples});
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, 0.0, 0.0, 1.0));
}`;
        }
    }
};
