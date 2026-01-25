// noise-texture.js - Noise texture generator with mirror-repeat addressing
export function createNoiseTextureGenerator(gpu, size = 256) {
    const wgsl = `
fn hash(p: vec2f) -> f32 { return fract(sin(dot(p, vec2f(12.9898, 78.233))) * 43758.5453); }
fn noise(p: vec2f) -> f32 {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    return mix(mix(hash(i), hash(i+vec2f(1,0)), u.x), mix(hash(i+vec2f(0,1)), hash(i+1.0), u.x), u.y);
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let p = vec2f(id.xy) / ${size}.0;
    let h = noise(p * 4.0) * 0.5 + 0.5;
    textureStore(out_0, vec2<i32>(id.xy), vec4f(h, h, h, 1.0));
}`;
    const algo = new (await import('./compute.js')).MapAlgorithm(gpu, wgsl, [], [new (await import('./gpu-core.js')).VisualField(gpu, 'NoiseTexture', 'rgba8unorm')]);
    return algo;
}
