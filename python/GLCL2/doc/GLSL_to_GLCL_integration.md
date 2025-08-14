# GLCL2 + pySymGLSL Integration Plan: Advanced OpenGL Features and Compute

This document proposes a concrete plan to integrate advanced OpenGL features (textures, framebuffers, instancing) into GLCL2, clarify compute shader support across GLCL2 and pySymGLSL, and align both frameworks around a shared pipeline model.


# Findings

- __pySymGLSL (moderngl)__
  - Uses `moderngl` and ships a small but capable pipeline in `python/pySymGLSL/GLSL_Simulation.py`.
  - Supports textures + samplers and per-target FBOs (ping-pong), see `GLSL_Simulation.setup_texture()` and `build_pipeline()` allocating textures and `self.ctx.framebuffer(...)`.
  - Render passes are fragment-shader based, executed by drawing a fullscreen quad via a VAO. Default vertex shader is injected if none is provided.
  - No compute shader abstraction present in `GLSL_Simulation.py`. Moderngl itself supports compute shaders, but the current class implements only VS/FS-based passes.
  - Default OpenGL context requirement is 3.3 core (see `moderngl.create_standalone_context(require=330)`).

- __GLCL2 (PyOpenGL + OpenCL)__
  - `python/GLCL2/OGLsystem.py` implements basic shader compilation (VS+FS), a minimal `GLobject` for VAO/VBO, and an `OGLSystem` registry for shader programs. No texture/fbo management yet.
  - `python/GLCL2/ModularGL.py` defines higher-level abstractions (`BufferObject`, `RenderObject`, `InstancedRenderObject`, `RenderSystem`). It references symbols not present in `OGLsystem.py` (e.g., `InstancedData`, `Mesh`, `create_sphere_mesh`), indicating incomplete or stale linkage. Instancing scaffolding exists but depends on missing pieces.
  - `python/GLCL2/GLCLBrowser.py` coordinates OpenCL kernels and OpenGL shaders, but contains no FBO/texture management and no OpenGL compute path.
  - OpenCL compute is present via `python/GLCL2/OCLsystem.py`.

- __GLCL NBody example__
  - `python/GLCL/NBody_glsl.py` uses OpenGL 4.3 compute shaders (SSBOs + `glDispatchCompute`) and a simple render pass of GL_POINTS. No geometry shader is used here. Context requested is 4.3 core.


# Goals

- __Unify pipeline concepts__ across pySymGLSL and GLCL2: common language for passes, resources (textures, FBOs, SSBOs), and scheduling.
- __Add textures + FBO management to GLCL2__ comparable to pySymGLSL’s usability (ping-pong, named targets, sampler binding, size management).
- __Stabilize instancing in GLCL2__: provide working `InstancedData`/`Mesh` layer or a minimal instancing helper that doesn’t rely on missing code.
- __Introduce OpenGL compute passes__
  - In pySymGLSL: optional compute-pass nodes leveraging moderngl’s compute support.
  - In GLCL2: OpenGL compute as a first-class pipeline node, alongside existing OpenCL kernels.
- __Synchronize resource sharing and barriers__ across compute and render passes (e.g., SSBO updates -> draw, textures written by FS -> next pass reads, CL/GL interop).


# Architecture Proposal

- __Common Resource Model__
  - Textures: 2D RGBA (default) with named registry and samplers. Add size management and auto-(re)allocation on resize.
  - Framebuffers: one FBO per named texture target; support MRT later.
  - SSBO/UBO: named SSBO registry with binding index assignment helper; utility to upload numpy arrays; typed views are responsibility of the shader.
  - Programs: program registry by name with stages (VS/FS/CS). Uniform caching optional.

- __Pass Abstractions__
  - RenderPass: (program=VS/FS, output=FBO name, inputs=[texture names], dynamic uniforms, draw=fullscreen quad or mesh).
  - ComputePass_GL: (program=CS, dispatch=(gx, gy, gz), buffers={binding: ssbo_name}, dynamic uniforms, optional memory barriers after dispatch).
  - ComputePass_CL: existing OpenCL kernel launch description in `GLCL2/OCLsystem.py` (no change in v1).

- __Scheduler__
  - Pipeline: sequential list of passes; each pass declares reads/writes. After each pass insert necessary barriers:
    - After CS writing SSBO -> `glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT)`
    - After FS writing textures -> implicit framebuffer completion; reading in next pass is fine. For safety between FS->CS reading images, add `GL_TEXTURE_FETCH_BARRIER_BIT`.
    - For CL/GL interop: ensure ownership and blocking transitions (glFinish or explicit glFlush+fences; CL-GL shared buffer acquire/release if used).

- __Versioning__
  - GLCL2: require OpenGL 4.3 when any OpenGL compute pass is present. Otherwise continue to support 3.3 for pure VS/FS.
  - pySymGLSL: keep 3.3 for VS/FS-only pipelines; enable optional compute nodes that require 4.3 when used.


# Detailed Changes by Package

## pySymGLSL

- __Add compute support__ to `python/pySymGLSL/GLSL_Simulation.py` (non-breaking):
  - `load_compute(name, path) -> moderngl.ComputeShader` and registry `self.compute: Dict[str, ComputeShader]`.
  - `bake_compute_pass(program_name, dispatch=(gx, gy, gz), buffers: Dict[int, moderngl.Buffer], dynamic_uniforms: List[str]) -> callable`.
  - Uniform update model mirrors `bake_pass()`; buffer binding via `buffer.bind_to_storage_buffer(binding)`. After dispatch, call `self.ctx.memory_barrier()`.
  - Extend `build_pipeline()` to accept pass-type tag e.g., tuples `("cs/my_comp.comp", None, {binding: ssbo_name}, [uniforms], "compute", (gx,gy,gz))`. Maintain backward compatibility by defaulting to render pass when `output` is a texture name.
  - Minimal buffer helpers: `create_ssbo(name, np_array)` -> `self.ctx.buffer(np_array.tobytes()); buffer.bind_to_storage_buffer(binding)` at pass-exec time.

- __Why here?__ Keeps pySymGLSL self-contained and allows users to experiment with compute pipelines without switching to PyOpenGL.

## GLCL2

- __Textures and FBOs__
  - New module `Textures.py` (or extend `OGLsystem.py`) providing:
    - `create_texture_2d(name, size, fmt=GL_RGBA32F, filter=GL_LINEAR, wrap=GL_CLAMP_TO_EDGE)` and `get_texture(name)`
    - `create_fbo(name, color_tex_names=[...])` with a default one-attachment path; registry `self.textures`, `self.fbos`, `self.samplers`.
    - Simple fullscreen quad VAO setup utility for FS passes.

- __Instancing__
  - Implement missing pieces referenced in `ModularGL.py`: either
    - Option A: add minimal `Mesh`, `InstancedData`, `create_sphere_mesh` into `OGLsystem.py` (keeping scope tight), or
    - Option B: refactor `ModularGL.py` to avoid external deps: implement an internal `InstancedData` that calls `glVertexAttribDivisor` and draws via `glDrawArraysInstanced`/`glDrawElementsInstanced`.
  - Recommend Option B to avoid duplication and reduce coupling.

- __OpenGL Compute Pass__
  - Add `GLComputeProgram` helper: compile CS, set uniforms, bind SSBOs, dispatch, and issue barriers.
  - Extend `GLCLBrowser` to recognize a pipeline step like:
    - `{"type":"compute_gl", "program":"my_cs.glsl", "dispatch":[gx,gy,gz], "ssbos": {0:"positions", 1:"velocities"}, "uniforms": {"dt": 0.01}}`
  - Keep existing OpenCL steps as-is; allow mixing order with required sync points.

- __Scheduler and Sync__
  - For GL->GL: insert `glMemoryBarrier` between CS writes and subsequent reads.
  - For CL<->GL: if sharing buffers, adopt CL-GL interop acquire/release flow; otherwise, explicit copies and `glFinish()` before CL and `clFinish()` before GL.


# Pipeline Unification Sketch

- __Common JSON idea__ (backward-compatible, illustrative):
  - Render pass (VS/FS):
    ```json
    ["prog/my_fs.glsl", "outTex", {"iChannel0":"prevTex"}, ["iResolution", "iFrame"]]
    ```
  - Compute pass (GL):
    ```json
    ["prog/my_cs.glsl", null, {"ssbo": {"0":"positions", "1":"velocities"}}, ["dt", "count"], "compute", [128, 1, 1]]
    ```
  - Compute pass (CL): existing GLCL2 kernel descriptor.

- __Rule__: if `output` is a texture/FBO name -> render pass; if `type == "compute"` -> compute pass with optional `dispatch`.


# Migration and Compatibility

- __pySymGLSL__ remains 3.3 FS-first. Compute is optional and gated by context version. Existing pipelines unaffected.
- __GLCL2__ gains FBO+texture management and optional GL compute passes. Existing OpenCL workflows remain intact.
- __NBody_glsl__ can be re-expressed as pipeline nodes: one `compute_gl` pass + one render pass (GL_POINTS), validating the abstraction.


# Work Plan (Phased)

- __Phase 1: Baseline + Docs (this file)__
  - Validate findings. Decide Option B for instancing.

- __Phase 2: GLCL2 textures/FBO__
  - Implement texture/FBO registry and fullscreen-quad utilities.
  - Add minimal render-pass executor (bind inputs, set uniforms, draw quad to target FBO).

- __Phase 3: Instancing stabilization__
  - Implement local `InstancedData` in `ModularGL.py` with `glVertexAttribDivisor` and `glDraw*Instanced`.
  - Provide a simple example in `tests_bash` or a small Python demo.

- __Phase 4: GL compute in GLCL2__
  - Add `GLComputeProgram`, SSBO helpers, and a `compute_gl` pipeline node in `GLCLBrowser`.
  - Convert `GLCL/NBody_glsl.py` into a pipeline config as a validation step.

- __Phase 5: pySymGLSL compute__
  - Add `load_compute` and `bake_compute_pass` in `GLSL_Simulation.py` mirroring `bake_pass` ergonomics.
  - Provide an example compute+render pipeline.

- __Phase 6: Interop and sync__
  - Define explicit sync rules and helper calls for GL<->CL steps (no-ops unless mixing).


# Risks and Mitigations

- __Moderngl vs PyOpenGL split__: Keep each framework self-contained. No cross-context resource sharing; instead, converge on common pipeline semantics.
- __Version requirements__: Guard compute usage behind GL 4.3 checks; fail loudly with clear messages if unavailable.
- __Stale references in GLCL2__: Remove or replace missing imports in `ModularGL.py` to avoid broken builds; keep minimal surface area.


# Immediate Action Items

- __GLCL2__
  - Add texture/FBO registry and fullscreen quad utilities in a new small module or in `OGLsystem.py`.
  - Implement `InstancedData` locally in `ModularGL.py` and remove missing imports.
  - Add `GLComputeProgram` and a `compute_gl` pipeline step in `GLCLBrowser` (no CL changes yet).

- __pySymGLSL__
  - Extend `GLSL_Simulation.py` with optional compute pass support using moderngl compute shaders.

- __Validation__
  - Recreate `NBody_glsl` via the new pipeline abstraction (one compute pass, one render pass) to verify design.


# Open Questions

- Should GLCL2 migrate to moderngl to reduce boilerplate and align with pySymGLSL, or keep PyOpenGL for tighter control and easy CL-GL interop?
- Do we want MRT and multi-pass down the line (e.g., deferred shading)? If so, extend FBO creation to multiple color attachments early.
- How should we express per-pass barriers in JSON (implicit vs explicit)?


# Appendix: File References

- `python/pySymGLSL/GLSL_Simulation.py` – textures/FBO/pipeline (FS), no compute yet.
- `python/GLCL2/OGLsystem.py` – VS/FS compile, `GLobject` basic VAO/VBO.
- `python/GLCL2/ModularGL.py` – high-level render + instancing scaffolding; references missing types.
- `python/GLCL/NBody_glsl.py` – OpenGL 4.3 compute + GL_POINTS rendering example.
