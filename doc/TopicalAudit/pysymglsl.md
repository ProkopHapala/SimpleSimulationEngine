---
type: TopicalAudit
title: pySymGLSL
tags: [topic, python, glsl, moderngl, gpgpu, render-pass, framebuffer, shader, simulation]
---

## Summary

Lightweight Python GPGPU simulation toolkit using GLSL fragment shaders via `moderngl`. Implements a "baked render graph" approach: users explicitly define textures, framebuffers, and render passes as simple tuple-based declarations. The system bakes all name lookups and object resolutions once at setup, then executes a pre-compiled list of `moderngl` handles in the hot loop with zero string lookups. Designed for grid-based physics simulations (fluid dynamics, particle systems) using ping-pong framebuffer techniques.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| Python | `python/pySymGLSL/GLSL_Simulation.py` | active | Core `GLSL_Simulation` class. Manages moderngl context, programs, textures, framebuffers. Bakes render passes into callable functions. Default vertex shader + fullscreen quad. |
| Python | `python/pySymGLSL/GLSL_GUI.py` | active | GUI integration. Playback control, viewport resize, shader graph rebuild from disk, uniform updates. Extends `BaseGUI`. |
| Doc | `doc/Markdown/pySymGLSL.md` | doc | Design discussion (960 lines): iterative design process from dictionary-based to baked graph approach. Final API: tuple-based `[shader, output, [inputs], [uniforms]]`. |

## Sub-topics

### Core Architecture

`GLSL_Simulation` class provides:
- **Context management**: `moderngl.create_context()`, fullscreen quad VAO
- **Resource dictionaries** (setup only): `programs`, `textures`, `framebuffers` — used during baking, not in hot loop
- **`create_texture(name, components, dtype)`**: Creates texture + corresponding framebuffer (FBO name = texture name)
- **`bake_pass(shader, output, inputs, uniforms)`**: Resolves names to `moderngl` object handles, returns baked callable
- **`run_graph(baked_passes, dynamic_uniforms)`**: Iterates over pre-baked passes — no lookups, no string parsing

### Baked Render Pass

A baked pass holds direct `moderngl` object references:
- `program` (compiled GLSL)
- `output_fbo` (framebuffer to write to)
- `inputs` (list of `(texture, texture_unit)` tuples)
- `dynamic_uniforms` (dict of uniform name → `moderngl.Uniform` object)

Execution: `fbo.use()` → `program.use()` → bind textures → set uniforms → `quad_fs.render()`

### Render Graph Definition

User defines passes as simple tuples:
```python
["shaders/advect.glsl", "velocity_b", ["velocity_a", "velocity_a"], ["u_deltatime"]]
```
Format: `[shader_path, output_texture, [input_textures], [dynamic_uniform_names]]`

User explicitly manages ping-pong — no automatic buffer swapping. Both A→B and B→A sequences defined separately.

### GUI Integration

`GLSL_GUI.py`:
- Playback start/stop control
- Viewport resize handling
- Shader pipeline parsing from text file (line-by-line)
- Recompilation of shaders from disk (hot-reload)
- Uniform value updates via GUI controls

### Design History

From `pySymGLSL.md` (960-line design doc):
1. Initial: dictionary-driven, automatic ping-pong → too slow, too magical
2. Refined: `RenderPass` dataclass with automatic ping-pong → user lacked control
3. Final: explicit user-managed ping-pong, tuple-based graph, baked execution → current implementation

Key design decisions:
- User is 100% responsible for texture naming and ping-pong
- No same texture as both input and output (OpenGL limitation)
- FBO name = texture name (simplification)
- Tuple format `[shader, output, [inputs], [uniforms]]` instead of dicts

## Parity Status

- **pySymGLSL ↔ C++ OpenCL simulations**: Different GPU computing approaches (GLSL fragment shader vs OpenCL kernels). Both use ping-pong buffer technique. No formal parity — different use cases (Python rapid prototyping vs C++ production).
- **pySymGLSL ↔ js/LandCraft_web WebGL**: Both use GLSL for GPGPU. pySymGLSL uses moderngl (desktop OpenGL), JS uses WebGL (browser). Different APIs but same concept.

## Open Issues

- No error handling for shader compilation failures in baked passes (moderngl raises but no custom error messages)
- No support for compute shaders (design decision, but limits flexibility for non-grid data)
- No texture format validation — user must ensure input/output texture formats are compatible
- No automatic resource cleanup (textures, framebuffers, programs not explicitly freed)
- `GLSL_GUI.py` shader pipeline parsing is fragile — format not well documented
- No tests for pySymGLSL functionality
- Design doc (`pySymGLSL.md`) is 960-line LLM chat transcript — useful for design rationale but not API reference
- No examples directory with ready-to-run simulations
