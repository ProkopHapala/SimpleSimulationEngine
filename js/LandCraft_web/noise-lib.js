// noise-lib.js - Periodic noise functions for terrain generation
export const NOISE_LIB = {
// Periodic coordinate wrapping helper
wrapCoord: `fn wrapCoord(p: vec2f, period: f32) -> vec2f {
    // WGSL lacks a vec2 mod; do scalar wrap
    let px = p.x - floor(p.x / period) * period;
    let py = p.y - floor(p.y / period) * period;
    return vec2f(px, py);
}`,

// 1. Sin-based periodic noise (simple, always periodic)
sinNoise: `
fn sinHash(p: vec2f) -> f32 { return fract(sin(dot(p, vec2f(12.9898, 78.233))) * 43758.5453); }
fn sinNoise(p: vec2f) -> f32 {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    return mix(mix(sinHash(i), sinHash(i+vec2f(1,0)), u.x), mix(sinHash(i+vec2f(0,1)), sinHash(i+1.0), u.x), u.y);
}`,

// 2. Value noise (periodic via integer hash)
valueNoise: `
fn valueHash(p: vec2f) -> f32 {
    let i = vec2<i32>(floor(p));
    var n: u32 = bitcast<u32>(i.x * 374761393 + i.y * 668265263);
    n = (n << 13u) ^ n;
    n = n * (n * n * 15731u + 789221u) + 1376312589u;
    return -1.0 + 2.0 * f32(n & 0x0fffffffu) / f32(0x0fffffffu);
}
fn valueNoise(p: vec2f) -> f32 {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    return mix(mix(valueHash(i), valueHash(i+vec2f(1,0)), u.x), mix(valueHash(i+vec2f(0,1)), valueHash(i+1.0), u.x), u.y);
}`,

// 3. Simplex noise (periodic via hash)
simplexNoise: `
fn simplexHash(p: vec2f) -> vec2f {
    let h = vec2f(dot(p, vec2f(127.1, 311.7)), dot(p, vec2f(269.5, 183.3)));
    return -1.0 + 2.0 * fract(sin(h) * 43758.5453);
}
fn simplexNoise(p: vec2f) -> f32 {
    let K1 = 0.366025404; let K2 = 0.211324865;
    let i = floor(p + (p.x + p.y) * K1); let a = p - i + (i.x + i.y) * K2;
    let m = step(a.y, a.x); let o = vec2f(m, 1.0 - m);
    let b = a - o + K2; let c = a - 1.0 + 2.0 * K2;
    let h = max(0.5 - vec3f(dot(a, a), dot(b, b), dot(c, c)), vec3f(0.0));
    let n = h * h * h * h * vec3f(dot(a, simplexHash(i + vec2f(0.0))), dot(b, simplexHash(i + o)), dot(c, simplexHash(i + vec2f(1.0))));
    return dot(n, vec3f(70.0));
}`,

// 4. Gradient noise (periodic via integer hash)
gradNoise: `
fn gradHash(p: vec2f) -> vec2f {
    let i = vec2<i32>(floor(p));
    var n: u32 = bitcast<u32>(i.x + i.y * 11111);
    n = (n << 13u) ^ n;
    n = (n * (n * n * 15731u + 789221u) + 1376312589u) >> 16u;
    n = n & 7u;
    let gr = vec2f(f32(n & 1u), f32(n >> 1u)) * 2.0 - 1.0;
    if(n >= 6u) { return vec2f(0.0, gr.x); }
    if(n >= 4u) { return vec2f(gr.x, 0.0); }
    return gr;
}
fn gradNoise(p: vec2f) -> f32 {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    return mix(mix(dot(gradHash(i), f), dot(gradHash(i+vec2f(1,0)), f-vec2f(1,0)), u.x), mix(dot(gradHash(i+vec2f(0,1)), f-vec2f(0,1)), dot(gradHash(i+1.0), f-vec2f(1,1)), u.x), u.y);
}`,

// 5. Orbit noise (periodic via hash)
orbitNoise: `
fn orbitHash(p: vec2f) -> vec4f {
    let px: u32 = bitcast<u32>(floor(p.x));
    let py: u32 = bitcast<u32>(floor(p.y));
    var n: u32 = py + 374761393u + px * 3266489917u;
    n = 2246822519u * (n ^ (n >> 15u));
    n = 3266489917u * (n ^ (n >> 13u));
    let rz = vec4u(n, n * 16807u, n * n * 48271u, n * n * 69621u);
    return vec4f(rz & vec4u(0x7fffffffu)) / f32(0x7fffffff);
}
fn nuttall(x: f32, w: f32) -> f32 {
    let pi = 3.14159265358979; if (abs(x) > w) { return 0.0; }
    return 0.365 - 0.5 * cos(pi * x / w + pi) + 0.135 * cos(2.0 * pi * x / w);
}
fn orbitNoise(p: vec2f) -> f32 {
    let ip = floor(p); let fp = fract(p); var rz = 0.0;
    for (var j = -1; j <= 2; j++) { for (var i = -1; i <= 2; i++) {
        let dp = vec2f(f32(i), f32(j));
        let rn = orbitHash(dp + ip) - 0.5;
        let op = fp - dp + rn.zw * 0.75;
        rz += nuttall(length(op), 1.85) * dot(rn.xy * 1.7, op);
    }}
    return rz * 0.5 + 0.5;
}`
};


// Noise flavor presets for mix-in
export const NOISE_FLAVORS = {
    sin: { label: 'Sin', body: NOISE_LIB.sinNoise, fn: 'sinNoise' },
    value: { label: 'Value', body: NOISE_LIB.valueNoise, fn: 'valueNoise' },
    simplex: { label: 'Simplex', body: NOISE_LIB.simplexNoise, fn: 'simplexNoise' },
    grad: { label: 'Gradient', body: NOISE_LIB.gradNoise, fn: 'gradNoise' },
    orbit: { label: 'Orbit', body: NOISE_LIB.orbitNoise, fn: 'orbitNoise' }
};