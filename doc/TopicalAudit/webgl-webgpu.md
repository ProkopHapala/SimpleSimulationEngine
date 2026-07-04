---
type: TopicalAudit
title: WebGL & WebGPU
tags: [topic, javascript, webgl, webgpu, wgsl, gpu-compute, browser, canvas, three-js]
---

## Summary

Browser-based GPU rendering and compute. WebGPU (WGSL) for LandCraft terrain generation with compute shaders, render pipelines, and texture I/O. WebGL (Three.js/GLSL) for MHD plasma visualization and SDF ray-marching. WebGPU N-body simulation. Architecture: `GPUContext` (device/canvas/uniforms), `MapAlgorithm` (compute pipeline), `ViewportRenderer` (fullscreen texture), `LineRenderer` (debug lines).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| JS+WGSL | `js/LandCraft_web/gpu-core.js` | active | `GPUContext`: WebGPU device init, canvas configuration, uniform buffer (GlobalParams struct), mouse/camera input, resize handling. `SHADER_HEADER` — shared WGSL uniform struct. |
| JS+WGSL | `js/LandCraft_web/compute.js` | active | `MapAlgorithm`: auto-generates WGSL bindings for input/output textures, creates compute pipeline, dispatches 8×8 workgroups. Supports `Field` and `VisualField` with `getReadView()`/`getWriteView()`. |
| JS+WGSL | `js/LandCraft_web/renderers.js` | active | `ViewportRenderer`: fullscreen triangle-strip, samples texture with camera transform. `LineRenderer`: line-list rendering with world-to-screen transform. `TextRenderer`: bitmap text. |
| JS+WGSL | `js/LandCraft_web/generators.js` | active | WGSL terrain generators. See `noise-procedural.md`. |
| JS+WGSL | `js/LandCraft_web/noise-lib.js` | active | WGSL noise functions. See `noise-procedural.md`. |
| JS+WGSL | `js/LandCraft_web/main.js` | active | Main loop: `frame()` — GPU update, compute pass (generation), colorize pass, render pass (terrain + lines + text), readback for stats. Progressive accumulation support. |
| JS+WGSL | `js/LandCraft_web/LandCraft.html` | active | Entry point: imports modules, initializes GPU, UI, generators, renderers |
| JS | `js/NBody2D_WebGPU/index_opt.html` | active | WebGPU N-body simulation: compute kernels for force calculation, render pipeline for particle visualization |
| JS+WebGL | `js/mhd_demo/render.js` | active | Three.js renderer for MHD plasma. See `mhd-plasma.md`. |
| JS+WebGL | `js/GLSL_solid_modeling/GLSLscreen.js` | active | WebGL screen-space SDF ray-marching. See `solid-modeling-csg.md`. |
| JS+WebGL | `js/GLSL_solid_modeling/ListOfPrimitives.html` | active | Interactive primitive list with SDF ray-marching |
| JS+WebGL | `js/FlowField/FlowField.HTML` | active | WebGL flow field visualization |
| JS+WebGL | `js/mhd_demo/main.js` | active | MHD demo entry point |
| Doc | `docs/MolGUI_web.md` | doc | Molecular editor web version design |
| Doc | `docs/JS_ESModule_Refactor_Plan.md` | doc | ES module refactoring plan for JS codebase |
| Doc | `js/LandCraft_web/doc/FractalTerrain.md` | doc | Terrain generation algorithms. See `noise-procedural.md`. |

## Sub-topics

### WebGPU Compute Pipeline

`MapAlgorithm` (compute.js):
- Auto-generates WGSL binding declarations from input/output texture lists
- Bind group layout: group(0) = globals, group(1) = textures
- Input textures: `texture_2d<f32>` (read), output: `texture_storage_2d<format, write>`
- Dispatch: `ceil(width/8) × ceil(height/8)` workgroups of 8×8
- Post-dispatch: `f.swap()` for double-buffered fields

### WebGPU Render Pipeline

`ViewportRenderer` (renderers.js):
- Fullscreen triangle-strip (4 vertices, no VBO)
- Vertex shader: reconstructs world position from camera uniforms, maps to texture UV
- Fragment shader: `textureSampleLevel(tex, samp, uv, 0.0)` — avoids non-uniform control flow errors
- Bounds check: returns background color outside `[0,1]` UV

`LineRenderer`:
- Vertex buffer with `float32x2` positions
- `worldToScreen()` common vertex shader
- Line-list topology, fixed color

### Uniform Management

`GPUContext` (gpu-core.js):
- `GlobalParams` struct: time, dt, frame, resX/Y, mouseX/Y, camX/Y, zoom, minHeight/maxHeight (64 bytes)
- `ArrayBuffer` with typed array views for direct write
- `device.queue.writeBuffer()` per frame

### Camera & Input

- Orthographic: `visibleHeight = 1000 / zoom`, `visibleWidth = visibleHeight * aspect`
- Mouse-to-world: NDC → world coordinates
- Wheel zoom: `zoom *= (1 - deltaY * 0.001)`, clamped `[0.05, 20]`
- Drag: camera pan scaled by zoom

### WebGL (Three.js)

MHD plasma (`render.js`):
- Orthographic camera
- Three.js meshes for coil rings
- GLSL fragment shader for B-field visualization (elliptic integrals in shader)
- HSV coloring: hue = field direction, value = magnitude

SDF ray-marching (`GLSLscreen.js`):
- Screen-space quad
- Fragment shader with `map()` function from scene tree
- Ray-march loop with distance threshold

## Parity Status

- **WebGPU (LandCraft) ↔ WebGL (MHD/SDF)**: Different APIs, different use cases. No parity needed.
- **WebGPU N-body ↔ Python PyOpenCL N-body**: Both GPU particle simulation. Different algorithms. See `parallel-particle-cell.md`.
- **WebGPU compute ↔ Python OpenCL compute**: Similar architecture (kernel + texture/buffer + dispatch). Different APIs (WGSL vs OpenCL C).

## Open Issues

- `GPUContext` has hardcoded `mapSize = {1024, 1024}` — not configurable
- `ViewportRenderer` uses `triangle-strip` with 4 vertices — draws 2 triangles (6 verts), 2 wasted
- `LineRenderer` has fixed color `vec4f(0.0, 0.5, 1.0, 1.0)` — no per-line color
- No WebGPU fallback to WebGL — `navigator.gpu` check only, throws if unsupported
- `MapAlgorithm` dispatches fixed 8×8 workgroups — no tuning for different GPUs
- No error handling for shader compilation failures in `MapAlgorithm`
- `NBody2D_WebGPU/index_opt.html` not fully audited
- No WebGPU compute shader support for particle dynamics (only terrain generation)
- `JS_ESModule_Refactor_Plan.md` suggests ongoing migration — mixed module systems
