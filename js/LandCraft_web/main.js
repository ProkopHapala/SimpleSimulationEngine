import { GPUContext, Field, VisualField } from './gpu-core.js';
import { MapAlgorithm } from './compute.js';
import { ViewportRenderer, LineRenderer } from './renderers.js';
import { TextRenderer } from './text-renderer.js';

// --- Shaders ---

const WGSL_GEN = `
fn hash(p: vec2f) -> f32 { return fract(sin(dot(p, vec2f(12.9898, 78.233))) * 43758.5453); }
fn noise(p: vec2f) -> f32 {
    let i = floor(p); let f = fract(p); let u = f*f*(3.0-2.0*f);
    return mix(mix(hash(i), hash(i+vec2f(1,0)), u.x), mix(hash(i+vec2f(0,1)), hash(i+1.0), u.x), u.y);
}
fn fbm(p: vec2f) -> f32 {
    var v=0.0; var a=0.5; var x=p;
    for(var i=0;i<6;i++){ v+=noise(x)*a; x*=2.0; a*=0.5; }
    return v;
}

@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let uv = vec2f(id.xy) / vec2f(globals.resX, globals.resY);
    
    var h = fbm(uv * 4.0 + vec2f(globals.time * 0.05, 0.0)); // Slowly moving noise
    h = h * 1.5 - 0.2;
    
    textureStore(out_0, vec2<i32>(id.xy), vec4f(max(h, 0.0), 0.0, 0.0, 1.0));
}
`;

const WGSL_COLORIZE = `
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let xy = vec2<i32>(id.xy);
    
    let h = textureLoad(in_0, xy, 0).r;
    
    // Nice Palette
    let deep = vec3f(0.01, 0.05, 0.2);
    let shallow = vec3f(0.0, 0.4, 0.6);
    let sand = vec3f(0.86, 0.8, 0.6);
    let grass = vec3f(0.1, 0.5, 0.2);
    let rock = vec3f(0.4, 0.35, 0.3);
    let snow = vec3f(0.95, 0.95, 1.0);
    
    var col = mix(deep, shallow, clamp(h*3.0, 0.0, 1.0));
    if(h > 0.3) { col = sand; }
    if(h > 0.35) { col = grass; }
    if(h > 0.6) { col = rock; }
    if(h > 0.8) { col = snow; }
    
    // Fake Normal mapping
    let hR = textureLoad(in_0, xy + vec2<i32>(1,0), 0).r;
    let hU = textureLoad(in_0, xy + vec2<i32>(0,1), 0).r;
    let ddx = (hR - h) * 40.0;
    let ddy = (hU - h) * 40.0;
    let n = normalize(vec3f(-ddx, -ddy, 1.0));
    let light = normalize(vec3f(-1.0, -1.0, 1.0));
    let diff = max(dot(n, light), 0.0);
    
    if(h > 0.3) { col *= (0.6 + 0.4*diff); } // Apply light only to land

    textureStore(out_0, xy, vec4f(col, 1.0));
}
`;

async function main() {
    const gpu = new GPUContext('gpuCanvas');
    await gpu.init();

    const heightMap = new Field(gpu, 'Height');
    const displayMap = new VisualField(gpu, 'Display', 'rgba8unorm');

    // 1. Terrain Generator
    // Note: We run this every frame to show animation, or once for static
    const genAlgo = new MapAlgorithm(gpu, WGSL_GEN, [], [heightMap]);
    
    // 2. Colorizer (Converts R32F -> RGBA8 for Rendering)
    const colorAlgo = new MapAlgorithm(gpu, WGSL_COLORIZE, [heightMap], [displayMap]);

    const terrainRenderer = new ViewportRenderer(gpu);
    const lineRenderer = new LineRenderer(gpu);
    const textRenderer = new TextRenderer(gpu);

    // Mock Data
    const rivers = new Float32Array([512, 512, 600, 600, 600, 600, 700, 550]);
    lineRenderer.updateData(rivers);

    function frame() {
        gpu.update();
        
        const cmd = gpu.device.createCommandEncoder();
        
        // Sim Steps
        genAlgo.run(cmd);   // Updates HeightMap
        colorAlgo.run(cmd); // Updates DisplayMap (Filterable)
        
        // Render Steps
        const pass = gpu.getRenderPass(cmd);
        
        // Draw Terrain (smoothly filtered)
        terrainRenderer.draw(pass, displayMap);
        
        // Draw Vectors
        lineRenderer.draw(pass);
        
        // Draw Text
        textRenderer.begin();
        textRenderer.addText(`Time: ${gpu.views.time[0].toFixed(2)}`, 10, 10);
        textRenderer.addText(`Zoom: ${gpu.camera.zoom.toFixed(2)}`, 10, 40);
        textRenderer.addText(`River Alpha`, 512, 512, 0.8);
        textRenderer.addText(`Target`, 700, 550, 0.8);
        textRenderer.draw(pass);
        
        pass.end();
        gpu.device.queue.submit([cmd.finish()]);
        
        requestAnimationFrame(frame);
    }
    requestAnimationFrame(frame);
}

main().catch(console.error);