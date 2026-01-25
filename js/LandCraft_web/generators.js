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
};
