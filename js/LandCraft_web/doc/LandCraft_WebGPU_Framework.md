## USER

I want to make simulation game called LandCraft in raw WebGPU (no use of Three.js or anything) which is based on operation of woth terrain height maps, and hydraulics (flooding, path finding, river basins etc.). Your task now it to make some general framework for that. The goal is not to implement individual algorithms (hydraulic errosion, A*, diffusion limited agregation, fluid dynamics, outflow, flooding), I already have copied thse algorithms (e.g. from Shadertoy). Your job is to thing nabout some general framework or helpers, which allow us to esaily make interface for many such algorithms without code redudancy/duplicity, it shioud be somewhat streamlined as all these algorithms woull be kinda similar working with 2D grids/textures/eventually particle buffers and list of indexes on that map (1D arrays of 2D coordinates which define polylines like pahts, river streams, railroads, etc.). Do not start to program yet just think ith throug, the general architecture and desing decisions.

Also figure out how to make general renderer which will be just 2D top-view, no isometric or perspective 3D view. This should make it easy. But generally maps will be huge (e.g. 16k*16k pixels), much larger than we can fit on screen so we will always view just part of it, resp zoom in zoom out, so we need some rendering shaders for resampling (subsampling as well as reduction - minimaps)...

---

## Gemini-3-pro

This is a great engineering challenge. Writing raw WebGPU allows for maximum performance, but without a framework, the verbose API (BindGroups, Layouts, Pipelines) will suffocate your prototyping speed.

Since you want to focus on **Heightmaps (Fields)** and **Polylines/Particles (Agents)** on a massive **16k** scale, here is a breakdown of the architecture and design decisions to create a flexible `LandCraft` framework.

---

### 1. The Core Philosophy: "Fields & Agents"
The engine should abstract data into two primary types. You should almost never manually create a `GPUBuffer` or `GPUTexture` in your game logic code.

1.  **Fields (The Grid):** These represent continuous data (Terrain Height, Water Depth, Sediment, Flux). Since simulations involve differential equations (Time $t$ depends on Time $t-1$), these must automatically implement **Ping-Pong (Double Buffering)**.
2.  **Agents (The Vector):** These represent discrete data (Droplets, River nodes, Pathfinding A* open-sets). These are stored in **Storage Buffers**.

---

### 2. Data Architecture

#### A. The `Field` Class (Texture Wrapper)
Handling 16k textures is VRAM intensive. The framework needs a standardized wrapper.

*   **Internal Structure:** Holds *two* `GPUTexture` handles (`Read` and `Write`).
*   **Swap Mechanism:** A `swap()` method that flips the read/write handles after a compute pass.
*   **Format abstraction:** Support for `R32Float` (Height), `RG32Float` (Velocity), and `RGBA8` (Visuals).
*   **Virtual Addressing (Optional but recommended):** For 16k maps, you might not want to simulate the whole map every frame. The `Field` should know about "Active Regions" (offsets and sizes) to pass to shaders.

#### B. The `PolyBuffer` Class (Structured Buffer Wrapper)
For paths and river streams.
*   **Structure:** A massive `GPUBuffer` used as `storage`.
*   **Indirect Draw/Dispatch:** Since the number of active particles/river nodes changes, this class must manage an **IndirectBuffer**. This allows the GPU to tell itself how many workgroups to dispatch in the next pass without CPU readback (which is slow).

#### C. The `GlobalUniforms`
Don't bind time, resolution, and mouse coordinates to every single shader individually. Create a single "Frame State" Uniform Buffer bound to `Group(0)` in *every* shader.
*   *Contents:* `Time`, `DeltaTime`, `MapResolution`, `CameraBounds`, `CursorPos`.

---

### 3. The Compute Pipeline (The "Algorithm Container")

This is where you save time on redundancy. Instead of writing full WebGPU pipelines for every erosion step, create a **`ComputeKernel`** builder.

#### Design Concept: The "Shader Injector"
You want to copy-paste Shadertoy-like code. Shadertoy works because it hides the boilerplate. Your framework should do the same.

**The boilerplate WGSL generator:**
Your framework should automatically prepend the header needed to interface with your `Field` objects.

*Input (Your Algorithm Code):*
```wgsl
// Your specific logic
fn simulate_water(pos: vec2<i32>) -> f32 {
    let current_h = textureLoad(heightMap, pos, 0).r;
    // ... complex hydraulic code ...
    return new_water_level;
}
```

*The Framework Wraps it into:*
```wgsl
@group(0) @binding(0) var<uniform> globals : GlobalParams;
@group(1) @binding(0) var heightMap : texture_2d<f32>; // Auto-bound Read texture
@group(1) @binding(1) var outputMap : texture_storage_2d<r32float, write>; // Auto-bound Write texture

@compute @workgroup_size(8, 8)
fn main(@builtin(global_invocation_id) id : vec3<u32>) {
    // Framework handles bounds checking
    if (id.x >= globals.res.x || id.y >= globals.res.y) { return; }
    
    // Call user function
    let val = simulate_water(vec2<i32>(id.xy));
    
    textureStore(outputMap, vec2<i32>(id.xy), vec4<f32>(val));
}
```

**The `Kernel` Class:**
*   **Constructor:** `new Kernel(shaderSource, inputs: [Field, Field], outputs: [Field])`
*   **Method `run(commandEncoder)`:**
    1.  Automatically creates BindGroups based on the input/output fields.
    2.  Dispatches workgroups based on map size.
    3.  If an output is a `Field`, it triggers the `swap()` on that field automatically after execution.

---

### 4. The Renderer (2D Top-Down View)

Rendering a 16k x 16k texture requires smart sampling. You cannot simply draw a 16k quad; it will cause massive aliasing (shimmering) when zoomed out.

#### A. The "Virtual Camera"
Don't use a perspective matrix. Use a simple 2D View transform:
*   `Offset` (World coordinate at center of screen).
*   `Scale` (World units per screen pixel).

#### B. Hierarchical Reduction (Auto-Mipmapping)
Since your simulation updates the texture every frame, hardware mipmaps might be costly or unavailable for `storage` textures depending on the format.
**The Solution: The "Reduction" Compute Pass.**
Before rendering, run a fast compute pass that downsamples the simulation data into a smaller texture (e.g., 1024x1024) specifically for "zoomed out" viewing.
*   **Zoom Level < X:** Render the raw 16k texture (viewport covers a small area).
*   **Zoom Level > X:** Render the cached downsampled texture (viewport covers the whole map).

#### C. The Visualization Pipeline
Do not make the simulation textures RGBA. Keep them float (R32 / RG32).
Create a dedicated **Colorizer Shader** (Render Pipeline).
*   **Input:** Height Field, Water Field, Sediment Field.
*   **Logic:**
    *   Calculate Normals on the fly (using derivatives `dpdx`, `dpdy` of height).
    *   Mix colors based on height (Grass/Rock/Snow).
    *   Overlay water based on depth (Blue tint + Specular).
    *   Output: Final screen color.
*   **Benefit:** You can change the "look" of the game without recalculating the simulation.

---

### 5. Interaction Strategy (CPU <-> GPU)

Since algorithms like A* might need to run on the GPU, getting data back to the CPU (e.g., to tell a unit where to move) is slow (`mapAsync`).

*   **Design Decision:** Keep the logic on the GPU as much as possible.
*   **Visual Debugging:** Instead of `console.log`, write values to a debug texture and render it.
*   **Mouse Interaction:** Pass the mouse world coordinates to the simulation.
    *   *Brush:* "If distance(id.xy, mouse.xy) < radius, add height."
    *   *Path:* "If click, add point to PolyBuffer."

---

### Summary of Classes/Modules to Implement

1.  **`SimulationContext`**: Initializes Device, Queue, Canvas.
2.  **`ResourceManager`**:
    *   `createField(width, height, format)` -> `PingPongTexture`
    *   `createAgentBuffer(structSize, count)` -> `StorageBuffer`
3.  **`ComputeSystem`**:
    *   `createKernel(wgslCode, bindings)`
    *   `dispatch(kernel)`
4.  **`Renderer`**:
    *   `setCamera(x, y, zoom)`
    *   `draw(fields...)` -> Runs the "Colorizer" shader to the canvas context.
5.  **`Loop`**:
    *   Updates Uniforms.
    *   Encodes Simulation Kernels.
    *   Encodes Render Pass.
    *   Submits.

### The "Game Loop" Flow

1.  **Input:** Update `GlobalUniforms` (Mouse pos, DeltaTime).
2.  **Sim 1 (Rain):** Read `Noise`, Write to `WaterField`.
3.  **Sim 2 (Flow):** Read `Height`, `Water`, Write to `VelocityField` + `WaterField` (PingPong).
4.  **Sim 3 (Erosion):** Read `Velocity`, Write `Height`, `Sediment`.
5.  **Render Prep:** Generate Min/Max/Average of heightmap (downsample) for zoomed-out view if necessary.
6.  **Render:** Draw Full-Screen Quad.
    *   Fragment Shader reads `Height`, `Water` (sampled with linear filtering).
    *   Apply lighting based on zooming level.

This architecture decouples the "How to run WebGPU" from "How water flows," allowing you to drop in Shadertoy algorithms by just mapping inputs and outputs.

---

## USER

OK, these are good ideas. now write it all. First write the renderer and general utils.

You can use this to generate intitial height map texture, and then render it, in next phase we add terrain errosion, but lets make the general framework with minimal working prototype fiorst before we start to add features.


---

This is a complete, single-file WebGPU prototype. It implements the architecture we discussed: **Fields** (Ping-Pong textures), **Kernels** (Simulation steps), and a **Renderer** (Top-down view with camera controls).

I have ported the noise logic from the "Sirenian Dawn" reference (replacing the texture lookups with procedural noise so it runs without external assets) to generate the initial terrain.

### Architecture Overview

1.  **`SimulationContext`**: Manages the GPU Device and Global Uniforms (Time, Mouse, Camera).
2.  **`Field`**: A wrapper around **two** textures (Read/Write). It handles the "Ping-Pong" logic automatically.
3.  **`Kernel`**: A reusable compute shader wrapper. You pass it WGSL code, and it handles the pipeline creation and dispatching. It auto-injects the framework's header code.
4.  **`Renderer`**: A dedicated render pipeline that draws a specific `Field` to the screen, handling Pan and Zoom.

### How to use
Save this code as an `.html` file (e.g., `landcraft.html`) and open it in a browser with WebGPU support (Chrome/Edge).

### Key Design Decisions Explained

1.  **The "Header Injection"**:
    *   Look at `SHADER_HEADER`. I defined a struct `GlobalParams` with binding `Group(0) Binding(0)`.
    *   Every shader (Compute or Render) gets this prepended. This means you **never** have to manually bind Time, Resolution, or Mouse in your algorithms. You just use `globals.time` or `globals.resX`.

2.  **The `Field` Class (Ping-Pong)**:
    *   Simulation algorithms (like hydraulic erosion) usually read from position `(x,y)` to calculate the next state. You cannot read and write to the same texture safely in parallel.
    *   The `Field` class holds `texA` and `texB`.
    *   The `Kernel.dispatch()` method automatically binds the *current* read texture to the shader's Input slot, and the *current* write texture to the Output slot.
    *   At the end of dispatch, it calls `field.swap()`. This logic is completely hidden from the "Game Logic" section.

3.  **Rendering Massive Maps (16k)**:
    *   I included a `mapSize` configuration (defaulting to 4096 for broad compatibility, but you can bump it to 16384).
    *   The Renderer uses a **Vertex Shader trick**. It doesn't modify the mesh (which is just 2 triangles). Instead, it calculates the UV coordinates based on `globals.camX`, `globals.camY`, and `globals.zoom`.
    *   **Zooming:** When you zoom out, the UV coordinates span a larger range (e.g., 0.0 to 1.0). When you zoom in, they span a tiny range (e.g., 0.45 to 0.55).
    *   **Subsampling:** The Fragment Shader uses `textureSampleLevel(..., 0.0)`. WebGPU hardware handles the linear interpolation. For extreme zoom-outs on 16k maps, this might shimmer slightly (aliasing) because we aren't generating mipmaps every frame (which is slow). This is the standard trade-off for real-time terrain editing.

4.  **Hydraulics Foundation**:
    *   I included `WGSL_EROSION_PLACEHOLDER`. Currently, it's just a blur kernel.
    *   **Next Step for You:** Replace the code inside `WGSL_EROSION_PLACEHOLDER` with the actual hydraulic erosion math (calculating water height, carrying capacity, sediment deposition). The framework is ready to accept it:
        1.  Add a new `waterField` (Field).
        2.  Update the Kernel to accept `[heightField, waterField]` as inputs/outputs.
        3.  Paste your Shadertoy logic into the WGSL string.

---

## USER

Can you please separate this into reusable modules.
1) Viewport renderer
2) line-rendeder (render the 1D structures like roads, rivers etc. 
3) particle rendered (render wolkeds on the teraon and along roads and rivers)
4) 2D map algorith (general class for any compute shader operating on 2D grid to 2D grid), like fluid dynamics, partical differential equations, hydraulic erosion, path finding etc.
5) main.js

try to make these low-level classes general and extendable so we can use them also in other projects. This will be just one application using this genral frameworks

---

## USER

look we have the proejct LandCraft described here
@LandCraft_WebGPU_Framework.md 

we started to implement it in
@LandCraft_main.html 
@main.js @renderers.js @compute.js @gpu-core.js 

but gemini left there many errors and unfinished things. Can you correct it? First make a map, how everything works so we have some overview, the general indeas and wehere is shich function what is doing what, like documentation one line buletpoint of purpose of each function how it fitrs into the biggir picture/whole framework context (not exhastive map of input output argument)

---

## GPT5-Codex

High-level module map
----------------------
1. `doc/LandCraft_WebGPU_Framework.md` – specification generated by Gemini describing the intended “Fields & Agents” architecture, ping-pong textures, kernel builder, renderer, and loop. Serves as blueprint for implementing the JS modules. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/doc/LandCraft_WebGPU_Framework.md#1-790
2. [gpu-core.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:0:0-0:0) – foundational runtime: global uniform header, GPU device/canvas setup, per-frame uniform updates, input handling, and the [Field](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:110:0-136:1) ping-pong texture wrapper. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#1-137
3. [compute.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:0:0-0:0) – minimal “kernel” abstraction ([MapAlgorithm](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:2:0-65:1)) that wraps a WGSL source string into a reusable compute pipeline with automatic binding of globals + field textures and ping-pong swaps. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js#1-66
4. [renderers.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/renderers.js:0:0-0:0) – currently empty; intended to house `ViewportRenderer`, `LineRenderer`, `ParticleRenderer` declared in [main.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js:0:0-0:0). Needs real implementations. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/renderers.js#1-1
5. [main.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js:0:0-0:0) – wiring script: instantiates GPU context, data fields, compute algorithms, placeholder renderers, dummy line/particle buffers, and the per-frame loop. Contains inline WGSL strings for terrain generation and fragment colorizer. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#1-118
6. `LandCraft_main.html` – (not re-read here) hosts canvas/UI and appears to contain another copy of the framework prototype (probable legacy). Needs reconciliation with JS modules.

Function/class inventory (per file)
-----------------------------------
### [gpu-core.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:0:0-0:0)
- `SHADER_HEADER`: WGSL snippet defining `GlobalParams` uniform struct bound at group(0)/binding(0), to be prepended to every shader so kernels/renderers share camera/time inputs. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#1-15
- [GPUContext(canvasId)](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:16:0-107:1) constructor: stores canvas handle, initializes device-related defaults (map size, camera/mouse state), allocates shared uniform buffer backing array views, and seeds timing counters. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#17-43
- [GPUContext.init()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:44:4-56:5): requests adapter/device, configures WebGPU canvas, creates uniform buffer, wires up input listeners. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#45-57
- [GPUContext.update()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:58:4-76:5): refreshes time/delta/frame counters, updates uniforms (resolution, mouse, camera, zoom), and writes them into GPU buffer each frame. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#59-77
- [GPUContext.getRenderPass(encoder)](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:78:4-86:5): convenience wrapper that sets up a color pass clearing the swapchain texture. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#79-87
- [GPUContext._setupInput()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:88:4-106:5): installs wheel zoom, mouse drag panning, and cursor tracking (currently in screen-relative coords). @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#89-107
- [Field(gpu, format)](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:110:0-136:1): wraps two GPU textures for ping-pong double buffering; stores map dimensions from context. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#110-127
- [Field.swap()/getReadView()/getWriteView()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:128:4-132:5): manages and exposes current read/write texture views to kernels. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js#129-136

### [compute.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:0:0-0:0)
- [MapAlgorithm(gpu, wgslCode, inputs, outputs)](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:2:0-65:1): generic compute kernel builder. Creates bind group layout (globals + texture inputs + storage outputs), compiles provided WGSL with shared header, builds compute pipeline. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js#3-39
- [MapAlgorithm.run(commandEncoder)](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:40:4-64:5): recreates bind group each frame using current field views, dispatches workgroups over map size (fixed 8×8 groups), ends pass, then swaps each output field to keep ping-pong consistent. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js#41-65

### [renderers.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/renderers.js:0:0-0:0)
- File is empty placeholder—`ViewportRenderer`, `LineRenderer`, `ParticleRenderer` referenced elsewhere are missing, explaining runtime errors. Needs implementations for fullscreen quad, polyline draw, and particle instancing. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/renderers.js#1-1

### [main.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js:0:0-0:0)
- WGSL string `GEN_TERRAIN`: compute shader (with custom helper functions `hash`, `noise`, `fbm`) generating island-height field; writes to storage texture bound manually at `@group(1) @binding(1)` (note mismatch with framework auto-binding). @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#7-31
- WGSL string `TERRAIN_COLOR_FS`: fragment shader sampling a height texture + sampler to colorize terrain by height bands. Depends on VS outputs defined inside yet-to-be-written renderer. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#33-57
- [main()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js:58:0-115:1):
  - Instantiates [GPUContext](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:16:0-107:1), awaits initialization. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#59-63
  - Creates one [Field](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:110:0-136:1) (`heightField`). @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#64-66
  - Builds [MapAlgorithm](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:2:0-65:1) using `GEN_TERRAIN` (no inputs, one output). @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#67-69
  - Instantiates placeholder renderers (`ViewportRenderer`, `LineRenderer`, `ParticleRenderer`) that currently do not exist. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#70-74
  - Runs generation kernel once. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#75-79
  - Seeds dummy river polyline and worker particle buffers, passes to renderers. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#80-93
  - Defines [frame()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js:93:4-113:5) loop:
    - Calls [gpu.update()](cci:1://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/gpu-core.js:58:4-76:5).
    - Creates command encoder.
    - [Placeholder] where simulation passes would be dispatched.
    - Begins render pass, calls renderer draw methods, submits, schedules next frame. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#94-116
  - Starts loop and handles errors. @/home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js#117-118

Notable gaps / next steps
-------------------------
1. [renderers.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/renderers.js:0:0-0:0) must implement the three renderer classes plus supporting vertex shader(s)/bind groups. Without this, [main.js](cci:7://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/main.js:0:0-0:0) throws import errors.
2. `GEN_TERRAIN` shader bypasses the framework’s automatic binding conventions (uses group(1)/binding(1) manually). Either adapt [MapAlgorithm](cci:2://file:///home/prokophapala/git/SimpleSimulationEngine/js/LandCraft_web/compute.js:2:0-65:1) to auto-inject bindings or rewrite the WGSL to use the standardized layout.
3. Need additional kernels (erosion, water, reduction) as outlined in the doc, along with camera-aware uniform population (convert mouse coords to world space).
4. `LandCraft_main.html` hosts a parallel prototype (duplicated logic). Decide which entry point to keep and reconcile assets.

Let me know if you want me to switch to Code mode to start fixing the missing renderers and reconcile the WGSL bindings.