# pySymGLSL Framework Documentation

## 1. Overview and Goal

pySymGLSL is a lightweight framework for building GPU-accelerated simulations purely with OpenGL and GLSL fragment shaders. Inspired by Shadertoy, it treats simulation state as 2D textures and composes multi-pass render graphs to advance the state each frame.

Goals:

- Focus on algorithm logic inside GLSL, not OpenGL boilerplate
- Define simulations declaratively via JSON pipelines
- Provide a small, debuggable Python core (`GLSL_Simulation`) and a minimal PyQt5 GUI (`GLSL_GUI`)

Typical use cases include grid/PDE-style simulations (fluids, diffusion, waves), image processing, or any ping‑pong texture compute.

## 2. Core Design Principles

- __Generality__: The engine is agnostic to specific simulations. Shader file names, uniform names, and texture names are free-form and provided in the pipeline JSON.

- __Baking__: Pipelines are parsed once, then “baked” into callables with pre-resolved uniforms, samplers, and draw state. The runtime loop is minimal.

- __Low Runtime Cost__: Per frame: update dynamic uniforms, bind textures/samplers, render full-screen quad to target FBO, repeat for each pass.

- __Shadertoy Compatibility__: Defaults mimic Shadertoy conventions: `iResolution`, `iFrame`, and `iChannelN` (N = 0..7). The default vertex shader provides `v_texcoord` in [0,1].

- __Simplicity__: One class (`GLSL_Simulation`) owns all ModernGL resources. GUI is kept thin and optional.

## 3. Quick Start

- __Run GUI__ (loads a default pipeline):

```bash
python -m pySymGLSL.GLSL_GUI
```

- __Load a pipeline__ in the GUI: Click “Load” and choose a JSON (e.g. `python/pySymGLSL/pipelines/fluid.json`).

- __Display texture__: Use the combo box “Display texture” to select which named texture to show.

- __Play / step__: Toggle “Play” to advance continuously. Updating parameters will re-run a frame when “auto” is enabled.

- __Programmatic use__ (headless OK):

```python
from pathlib import Path
from pySymGLSL import GLSL_Simulation

sim = GLSL_Simulation(sim_size=(512, 512))
base_dir = Path(__file__).parent / "python/pySymGLSL/shaders"

# Example pipeline: [fragment_path, output_tex, inputs_map, dynamic_uniforms]
pipeline = [
    ["common/gauss.glslf", "view0", {}, ["pos0"]],
]

baked, tex_names = sim.build_pipeline(pipeline, base_dir)
# Set dynamic uniforms each frame
values = {
    "iResolution": (512.0, 512.0, 1.0),
    "iFrame": 0,
    "pos0": (0.0, 0.0, 0.2, 1.0)
}
sim.run_graph(baked, values)
```

## 4. JSON Pipeline Schema

Pipelines live in `python/pySymGLSL/pipelines/`. Examples: `fluid.json`, `gauss.json`.

Top-level fields:

- `parameters`: map name → [type, default(s), step]
  - Type is an informational string (e.g., `"vec4"`, `"float"`). The GUI ignores this field; the number of components is inferred from the length of the defaults list.
  - Defaults: single number or list (for vectors). Example: `[0.5, 0.5, 0.0, 1.0]`.
  - Step: spin box step size.

- `Pipeline`: array of passes. Each pass is a 4‑tuple:
  1) `fragment_path` (relative to `python/pySymGLSL/shaders/`),
  2) `output_texture_name` (created automatically if missing),
  3) `inputs_map`: e.g. `{ "iChannel0": "fieldA", "iChannel1": "fieldB" }`,
  4) `dynamic_uniforms`: list of uniform names the engine should update from the per-frame dictionary.

Example (`fluid.json`):

```json
{
  "parameters": {
    "driver": ["vec4", [0.5, 0.5, 0.0, 1.0], 0.1]
  },
  "Pipeline": [
    ["fluid/solveFluid.glslf", "fieldA", {"iChannel0": "fieldC"}, ["driver"]],
    ["fluid/solveFluid.glslf", "fieldB", {"iChannel0": "fieldA"}, ["driver"]],
    ["fluid/solveFluid.glslf", "fieldC", {"iChannel0": "fieldB"}, ["driver"]],
    ["fluid/view.glslf",      "view0",  {"iChannel0": "fieldC"}, []]
  ]
}
```

Notes:

- All textures referenced in the pipeline (as outputs or in `iChannelN`) are auto‑allocated with size `sim_size`, components=4, dtype from `GLSL_Simulation(dtype=...)`.
- Textures are initialized once via `initialize_textures()` (default clear color 0,0,0,0) to avoid undefined feedback behavior.
- `iResolution` and `iFrame` are automatically appended to the `dynamic_uniforms` list if missing, so you can omit them in JSON and still receive updates.

## 5. Shader Conventions

- __Vertex shader__: If none is provided, a default full-screen quad VS is injected. It writes `v_texcoord` in [0,1].
- __Fragment uniforms__ commonly used:
  - `uniform vec3 iResolution;` // (width, height, 1.0)
  - `uniform int  iFrame;`
  - `uniform sampler2D iChannel0;` .. `iChannel7;`
  - Any custom uniforms listed in the pass `dynamic_uniforms` (e.g. `driver`, `pos0`).

- __Accessing textures__: Use Shadertoy-style sampling with `v_texcoord` or compute `uv = gl_FragCoord.xy / iResolution.xy`.

- __Outputs__: Write next state to `gl_FragColor` (RGBA). Multi-channel textures are supported; the meaning is simulation-defined.

See examples:

- `shaders/fluid/solveFluid.glslf` (multi-pass stable fluid with vorticity)
- `shaders/fluid/view.glslf` (visualization)
- `shaders/common/gauss.glslf` (simple analytic field)

## 6. Key Components (Developers)

- __`GLSL_Simulation`__ (`python/pySymGLSL/GLSL_Simulation.py`)
  - Manages ModernGL context, programs, textures, samplers, and FBOs
  - Methods:
    - `load_program(name, vertex_path=None, fragment_path)`: compile/link program; default VS used if none given. Comments are stripped from fragment source to avoid false includes.
    - `build_pipeline(pipeline, base_dir)`: ensure programs are loaded; auto-allocate textures/FBOs; return `(baked_graph, texture_names)`
    - `bake_pass(program_name, output_name, input_names, dynamic_uniforms)`: pre-resolve uniform setters, bind texture units (`iChannelN` or `u_texture_N`), create a pass VAO, and return an executable callable
    - `bake_graph(graph_definition)`: list of callables from a list of pass descriptors
    - `initialize_textures(v=(0,0,0,0))`: clear all FBOs to a value (called inside `build_pipeline`)
    - `run_graph(baked_graph, dynamic_values)`: execute all pass callables; increments `iFrame`
    - `release()`: free GL resources (programs, textures, samplers, FBOs, VAOs)

- __`GLRenderWidget`__ (`python/pySymGLSL/GLSL_GUI.py`)
  - QGLWidget that owns a `GLSL_Simulation` and an on-screen display program sampling a texture
  - In `paintGL()`: ensures dynamic uniforms (adds `iResolution`/`iFrame`), runs the baked graph, then displays the selected texture to the default framebuffer

- __`MainWindow`__ (`python/pySymGLSL/GLSL_GUI.py`)
  - Lightweight PyQt5 GUI for loading/saving pipelines, editing parameters, selecting the display texture, and running the simulation
  - Uses helpers from `BaseGUI.py` (widget factory + JSON utilities)

- __`BaseGUI`__ (`python/pySymGLSL/BaseGUI.py`)
  - Reusable widget builders: `button`, `checkBox`, `comboBox`, `spinBox`, `textEdit`, `spin_row`
  - JSON helpers: `strip_json_comments()`, `extract_json_block()`
  - `populate_params_from_json(params_dict)`: creates scalar/vector spin boxes from the `parameters` schema

## 7. Lifecycle

- __Load Pipeline__ (GUI `load_pipeline`):
  1. Read JSON, strip comments
  2. Compute `base_dir = <pipeline_dir>/../shaders`
  3. `GLSL_Simulation.build_pipeline(Pipeline, base_dir)`
  4. Update UI parameter widgets from `parameters`

- __Initialize GL__ (`GLRenderWidget.initializeGL`):
  - Create ModernGL context
  - Create on‑screen display program and VAO (separate from simulation passes)

- __Per Frame__ (`GLRenderWidget.paintGL`):
  1. Ensure `iResolution` and `iFrame` are present in `dynamic_values`
  2. `sim.run_graph(baked_graph, dynamic_values)`
  3. Bind chosen texture and draw to screen using display VAO

## 8. Debugging Guide

- __Shader compile errors__: Raised as `RuntimeError` in `load_program()` with ModernGL error log and fragment path.
- __Missing uniforms__: `bake_pass()` validates uniform presence and raises `KeyError` if a dynamic uniform is not found in the program.
- __Uniform arity__: A mismatch (e.g. sending scalar where vec3 expected) raises `ValueError` in the baked setter.
- __No output drawn__: Check that you selected an existing texture in the “Display texture” combo. Use the first item as a sanity check.
- __Feedback artifacts__: Ensure UVs in `[0,1]` (e.g. `gl_FragCoord.xy / iResolution.xy`), correct sampler usage, and that initial textures are properly cleared.
- __Performance__: Keep passes minimal; only list uniforms that change per frame in `dynamic_uniforms`.

## 9. Examples

- __Gaussian field__ (`pipelines/gauss.json`):

```json
{
  "parameters": {
    "pos0": ["vec4", [0.0, 0.0, 0.2, 1.0], 0.1]
  },
  "Pipeline": [
    ["common/gauss.glslf", "view0", {}, ["pos0"]]
  ]
}
```

- __Simple fluid__ (`pipelines/fluid.json`): triplicated advection/force solve feeding into a view pass. See `shaders/fluid/solveFluid.glslf` and `shaders/fluid/view.glslf`.

## 10. Integration Notes (with GLCL2)

- GLCL2 combines OpenCL kernels (compute) and OpenGL shaders (rendering); pySymGLSL focuses solely on OpenGL/GLSL multi-pass compute.
- Both systems use a “baked pipeline” idea and Shadertoy-like shader conventions. It’s straightforward to prototype visual GLSL pipelines in pySymGLSL and later port parts to GLCL2 render stages if needed.
- For shared GUI patterns, `BaseGUI.py` provides reusable widget helpers and JSON utilities that mirror the approach used in GLCL2 (see `GLCL2/doc/GLCL_manifest.md` for broader architecture).

## 11. Repository Layout

```
python/pySymGLSL/
├── BaseGUI.py
├── GLSL_GUI.py
├── GLSL_Simulation.py
├── __init__.py
├── doc/
│   └── pySymGLSL_manifest.md   ← this file
├── pipelines/
│   ├── fluid.json
│   ├── fluid_old.json (legacy format)
│   └── gauss.json
└── shaders/
    ├── common/
    │   ├── copy.glslf
    │   └── gauss.glslf
    ├── fluid/
    │   ├── solveFluid.glslf
    │   └── view.glslf
    └── shader_toy_ref/ (reference examples)
```

## 12. Best Practices

- Use `gl_FragCoord.xy / iResolution.xy` for UVs unless you need custom mapping
- Keep `dynamic_uniforms` minimal; update from GUI or your Python loop
- Prefer `iChannelN` names for texture samplers to match defaults
- Initialize state explicitly in the first frames (see `iFrame` usage in fluid example)

