import { GPUContext, Field, VisualField, StaticTexture } from './gpu-core.js';
import { MapAlgorithm } from './compute.js';
import { ViewportRenderer, LineRenderer } from './renderers.js';
import { TextRenderer } from './text-renderer.js';
import { GENERATORS } from './generators.js';
import { buildParamControls } from './gui-builder.js';
import { NOISE_FLAVORS } from './noise-lib.js';

function makeRgbaNoise(size, channels = 1, seed = 1337) {
    let s = seed >>> 0;
    const rnd = () => { s = (1664525 * s + 1013904223) >>> 0; return s / 4294967296; };
    const data = new Uint8Array(size * size * 4);
    for (let i = 0; i < size * size; i++) {
        const r = Math.floor(rnd() * 256);
        const g = channels >= 2 ? Math.floor(rnd() * 256) : r;
        const b = channels >= 3 ? Math.floor(rnd() * 256) : r;
        data[i * 4 + 0] = r;
        data[i * 4 + 1] = g;
        data[i * 4 + 2] = b;
        data[i * 4 + 3] = 255;
    }
    return data;
}

// --- Colormap presets and shader builder ---
const COLORMAPS = {
    grayscaleIso: `fn colorMap(h: f32) -> vec3f { let base = vec3f(h); let lines = smoothstep(0.0, 0.002, abs(fract(h*20.0) - 0.5)); let iso = mix(base, vec3f(1.0, 0.8, 0.4), lines*0.35); return iso; }`,
    heat: `fn colorMap(h: f32) -> vec3f { let t = clamp(h,0.0,1.0); return mix(vec3f(0.05,0.02,0.1), vec3f(1.0,0.9,0.2), t); }`,
    terrain: `fn colorMap(h: f32) -> vec3f { let sea = mix(vec3f(0.02,0.08,0.25), vec3f(0.0,0.4,0.6), smoothstep(0.0,0.3,h)); let land = mix(vec3f(0.12,0.5,0.2), vec3f(0.5,0.4,0.35), smoothstep(0.35,0.7,h)); let snow = mix(land, vec3f(0.95,0.95,1.0), smoothstep(0.7,0.95,h)); var col = snow; if(h < 0.3) { col = sea; } return col; }`
};

function buildColorizeShader(colormapBody) {
    return `
${colormapBody}
fn isoOverlay(h: f32, freq: f32, sharp: f32) -> f32 {
    // exp(-B*sin) comb to sharpen isolines
    let s = sin(h * freq * 6.28318530718);
    return exp(-sharp * abs(s));
}
@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    let xy = vec2<i32>(id.xy);
    let invRange = 1.0 / max(globals.maxHeight - globals.minHeight, 1e-6);
    let hRaw = textureLoad(in_0, xy, 0).r;
    var h = clamp((hRaw - globals.minHeight) * invRange, 0.0, 1.0);
    var col = colorMap(h);
    let iso = isoOverlay(h, 10.0, 8.0);
    col = mix(col, vec3f(1.0, 0.2, 0.1), iso * 0.25);
    let hR = clamp((textureLoad(in_0, xy + vec2<i32>(1,0), 0).r - globals.minHeight) * invRange, 0.0, 1.0);
    let hU = clamp((textureLoad(in_0, xy + vec2<i32>(0,1), 0).r - globals.minHeight) * invRange, 0.0, 1.0);
    let ddx = (hR - h) * 40.0;
    let ddy = (hU - h) * 40.0;
    let n = normalize(vec3f(-ddx, -ddy, 1.0));
    let light = normalize(vec3f(-1.0, -1.0, 1.0));
    let diff = max(dot(n, light), 0.0);
    col *= (0.55 + 0.45 * diff);
    textureStore(out_0, xy, vec4f(col, 1.0));
}
`; }

async function main() {
    const gpu = new GPUContext('gpuCanvas');
    await gpu.init();

    const heightMap = new Field(gpu, 'Height');
    const displayMap = new VisualField(gpu, 'Display', 'rgba8unorm');
    const accumBuffer = new Field(gpu, 'Accumulator');

    const clearShader = `@compute @workgroup_size(8,8)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x >= u32(globals.resX) || id.y >= u32(globals.resY)) { return; }
    textureStore(out_0, vec2<i32>(id.xy), vec4f(0.0, 0.0, 0.0, 0.0));
}`;
    const clearAlgo = new MapAlgorithm(gpu, clearShader, [], [accumBuffer]);

    function resetAccumulator() {
        const cmd = gpu.device.createCommandEncoder();
        clearAlgo.run(cmd);
        gpu.device.queue.submit([cmd.finish()]);
    }

    // Reference noise lookup textures (Shadertoy-style iChannel)
    const noise256 = new StaticTexture(gpu, 'Noise256', 256, 256, 'rgba8unorm');
    noise256.uploadRGBA8(makeRgbaNoise(256, 1, 1337));
    const noise1024 = new StaticTexture(gpu, 'Noise1024', 1024, 1024, 'rgba8unorm');
    noise1024.uploadRGBA8(makeRgbaNoise(1024, 2, 424242));

    // Setup UI
    const select = document.getElementById('generator-select');
    const colormapSelect = document.getElementById('colormap-select');
    const generateBtn = document.getElementById('generate-btn');
    const paramsContainer = document.createElement('div');
    paramsContainer.id = 'params-container';
    paramsContainer.style.marginTop = '10px';
    document.getElementById('sidebar').appendChild(paramsContainer);
    const generatorKeys = Object.keys(GENERATORS);
    generatorKeys.forEach(key => {
        const opt = document.createElement('option');
        opt.value = key;
        opt.textContent = GENERATORS[key].label;
        select.appendChild(opt);
    });

    const colormapKeys = Object.keys(COLORMAPS);
    colormapKeys.forEach(key => {
        const opt = document.createElement('option');
        opt.value = key;
        opt.textContent = key;
        colormapSelect.appendChild(opt);
    });

    const generatorParams = {};
    const generatorFlavors = {};
    function initParams(gen) {
        if (!generatorParams[gen.id]) {
            const obj = {}; gen.params?.forEach(p => obj[p.name] = p.default); generatorParams[gen.id] = obj;
        }
        return generatorParams[gen.id];
    }
    function initFlavor(gen) {
        if (!gen.noiseFlavor) return undefined;
        if (!generatorFlavors[gen.id]) generatorFlavors[gen.id] = gen.noiseFlavor;
        return generatorFlavors[gen.id];
    }

    let currentGenerator = GENERATORS[generatorKeys[0]];
    let currentParams = initParams(currentGenerator);
    let currentFlavor = initFlavor(currentGenerator);
    select.value = currentGenerator.id;
    let currentColormap = colormapKeys[0];
    colormapSelect.value = currentColormap;
    const extraInputsFor = (gen) => (gen.extraInputs ? gen.extraInputs({ noise256, noise1024 }) : []);
    let genAlgo = new MapAlgorithm(gpu, currentGenerator.wgsl(currentParams, currentFlavor), extraInputsFor(currentGenerator), [heightMap]);
    let colorAlgo = new MapAlgorithm(gpu, buildColorizeShader(COLORMAPS[currentColormap]), [heightMap], [displayMap]);
    gpu.views.heightRange[0] = 0.0;
    gpu.views.heightRange[1] = 1.0;
    let needsGeneration = false;
    let needsColorize = true;

    const terrainRenderer = new ViewportRenderer(gpu);
    const lineRenderer = new LineRenderer(gpu);
    const textRenderer = new TextRenderer(gpu);

    const rivers = new Float32Array([512, 512, 600, 600, 600, 600, 700, 550]);
    lineRenderer.updateData(rivers);

    function rebuildGeneratorShader() {
        const flavorKey = initFlavor(currentGenerator);
        const wgsl = currentGenerator.wgsl.length >= 2 ? currentGenerator.wgsl(currentParams, flavorKey) : currentGenerator.wgsl(currentParams);
        const inputs = extraInputsFor(currentGenerator);
        const outputs = [heightMap];
        const wasProgressive = isProgressive;
        isProgressive = currentGenerator.progressive && currentParams.progressive;
        if (isProgressive) {
            inputs.push(accumBuffer);
            outputs.push(accumBuffer);
            if (!wasProgressive || needsAccumReset) { resetAccumulator(); }
            progressiveFrame = 0;
            progressiveTotalFrames = Math.ceil(currentParams.samples / currentParams.samplesPerFrame);
            needsAccumReset = false;
        }
        genAlgo = new MapAlgorithm(gpu, wgsl, inputs, outputs);
    }

    function updateGenerator() {
        paramsContainer.innerHTML = '';
        if (currentGenerator.noiseFlavor) {
            const flavorRow = document.createElement('div');
            flavorRow.style.marginBottom = '10px';
            const label = document.createElement('label');
            label.textContent = 'Noise Flavor:';
            label.style.display = 'block';
            label.style.fontSize = '12px';
            flavorRow.appendChild(label);
            const selectEl = document.createElement('select');
            Object.entries(NOISE_FLAVORS).forEach(([key, meta]) => {
                const opt = document.createElement('option');
                opt.value = key;
                opt.textContent = meta.label;
                selectEl.appendChild(opt);
            });
            selectEl.value = initFlavor(currentGenerator);
            selectEl.addEventListener('change', () => {
                generatorFlavors[currentGenerator.id] = selectEl.value;
                currentFlavor = selectEl.value;
                rebuildGeneratorShader();
                requestGeneration();
            });
            selectEl.style.width = '100%';
            flavorRow.appendChild(selectEl);
            paramsContainer.appendChild(flavorRow);
        }
        if (currentGenerator.params) buildParamControls(paramsContainer, currentGenerator.params, (name, value) => {
            currentParams[name] = value;
            if (currentGenerator.progressive) needsAccumReset = true;
            rebuildGeneratorShader();
            requestGeneration();
        }, currentParams);
    }

    function requestGeneration() { needsGeneration = true; }

    generateBtn.addEventListener('click', () => {
        requestGeneration();
    });

    select.addEventListener('change', () => {
        const key = select.value;
        currentGenerator = GENERATORS[key];
        currentParams = initParams(currentGenerator);
        currentFlavor = initFlavor(currentGenerator);
        rebuildGeneratorShader();
        updateGenerator();
        requestGeneration();
    });

    colormapSelect.addEventListener('change', () => {
        currentColormap = colormapSelect.value;
        colorAlgo = new MapAlgorithm(gpu, buildColorizeShader(COLORMAPS[currentColormap]), [heightMap], [displayMap]);
        needsColorize = true;
    });

    // Initial generation
    updateGenerator();
    requestGeneration();

    let heightStatsPromise = null;
    let lastHeightStats = { min: 0.0, max: 1.0 };
    const enableAutoHeightStats = true;
    let lastGenTime = 0;
    let progressiveFrame = 0;
    let progressiveTotalFrames = 0;
    let isProgressive = false;
    let needsAccumReset = false;

    function autoStatsEnabled() {
        return enableAutoHeightStats && currentGenerator.id !== 'procErosionHeightmap_ref';
    }

    function kickHeightStats() {
        if (!autoStatsEnabled()) return;
        if (isProgressive && progressiveFrame < progressiveTotalFrames) return;
        if (heightStatsPromise || !heightMap.pendingReadback) return;
        heightStatsPromise = heightMap.collectReadback().then(data => {
            heightStatsPromise = null;
            if (!data) return;
            let min = data[0]; let max = data[0];
            for (let i = 1; i < data.length; i++) {
                const v = data[i];
                if (v < min) min = v;
                if (v > max) max = v;
            }
            if (!Number.isFinite(min) || !Number.isFinite(max)) return;
            if (max - min < 1e-6) { max = min + 1e-6; }
            gpu.views.heightRange[0] = min;
            gpu.views.heightRange[1] = max;
            lastHeightStats.min = min;
            lastHeightStats.max = max;
            console.log('Height stats', { min, max });
            needsColorize = true;
        }).catch(err => {
            console.error('Height readback failed', err);
            heightStatsPromise = null;
        });
    }

    function frame() {
        gpu.update();
        const cmd = gpu.device.createCommandEncoder();
        if (needsGeneration) {
            const t0 = performance.now();
            genAlgo.run(cmd);
            const t1 = performance.now();
            console.log(`Compute pass: ${(t1 - t0).toFixed(2)} ms`);
            lastGenTime = t1 - t0;

            if (isProgressive) {
                progressiveFrame++;
                const progress = progressiveFrame / progressiveTotalFrames;
                console.log(`Progressive: ${progressiveFrame}/${progressiveTotalFrames} (${(progress * 100).toFixed(1)}%)`);
                if (progressiveFrame >= progressiveTotalFrames) {
                    isProgressive = false;
                    console.log('Progressive accumulation complete');
                    if (autoStatsEnabled() && !heightMap.pendingReadback) { heightMap.encodeReadback(cmd); }
                    needsGeneration = false;
                } else {
                    needsGeneration = true;
                }
            } else {
                if (autoStatsEnabled() && !heightMap.pendingReadback) { heightMap.encodeReadback(cmd); }
                needsGeneration = false;
            }
            needsColorize = true;
        }
        if (needsColorize) {
            colorAlgo.run(cmd);
            needsColorize = false;
        }
        const pass = gpu.getRenderPass(cmd);
        terrainRenderer.draw(pass, displayMap);
        lineRenderer.draw(pass);
        textRenderer.begin();
        textRenderer.addText(`Time: ${gpu.views.time[0].toFixed(2)}`, 10, 10);
        textRenderer.addText(`Zoom: ${gpu.camera.zoom.toFixed(2)}`, 10, 40);
        textRenderer.addText(`Height min ${lastHeightStats.min.toFixed(3)} max ${lastHeightStats.max.toFixed(3)}`, 10, 70);
        textRenderer.addText(`River Alpha`, 512, 512, 0.8);
        textRenderer.addText(`Target`, 700, 550, 0.8);
        textRenderer.draw(pass);
        pass.end();
        gpu.device.queue.submit([cmd.finish()]);
        kickHeightStats();
        requestAnimationFrame(frame);
    }
    requestGeneration();
    requestAnimationFrame(frame);
}

main().catch(console.error);