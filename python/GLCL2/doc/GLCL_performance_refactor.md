# GLCL2 Browser: Simplification and Performance Refactor Plan

This document audits the per-frame/per-step paths in GLCL2 and proposes concrete simplification and optimization actions. It also lists potential pitfalls to avoid when implementing them.

Scope reviewed:
- `python/GLCL2/GLCLBrowser.py`
- `python/GLCL2/GLCLGUI.py`
- `python/GLCL2/OGLsystem.py`
- `python/GLCL2/OCLsystem.py`
- `python/BaseGUI.py` (for UI refactors)


## Critical per-frame paths

- `GLCLBrowser.update_simulation()`
  - Loops kernels: `for kernel_call in self.baked_kernel_calls: kernel_call.execute()`
  - Buffer sync CPU path: for each `buf` in `self.buffers_to_sync`: `np.empty(..)` → `OCLSystem.fromGPU(..)` → `GLCLWidget.update_buffer_data(..)`
  - `self.glcl_widget.update()` schedules repaint.

- `GLCLWidget.paintGL()`
  - `_execute_fs_pipeline()` for full-screen passes.
  - Iterates `render_pipeline_info` to draw GL objects (program lookup, uniform sets, draw).
  - Optional display overlay `_draw_display_texture()`.


## High-impact, low-risk optimizations

- __Preallocate host buffers for CL→GL sync__
  - Today: `np.empty(self.buffer_shapes[buf_name], dtype=np.float32)` per frame per buffer in `GLCLBrowser.update_simulation()`.
  - Action: allocate once in `apply_simulation_config()` (or after `_precompute_buffer_sync_list()`), store in a dict `self.host_sync_buffers[name]`. Reallocate only when shape changes.
  - Risk: must re-init when parameters alter buffer sizes; hook into `reset_simulation()`/rebake.

- __Replace `glBufferData` with `glBufferSubData` when size unchanged__
  - `OGLsystem.GLobject.upload_vbo()` always calls `glBufferData` (realloc path). Reallocation can stall the driver.
  - Action: track current VBO size in `GLobject`; call `glBufferSubData` if `nbytes` unchanged, else `glBufferData`.
  - Risk: none if size is tracked correctly.

- __Cache uniform locations and program handles__
  - `GLCLWidget._execute_fs_pipeline()` performs `glGetUniformLocation` for `iChannelN`, `iFrame`, `iResolution`, `driver` every pass, every frame.
  - `GLCLWidget.update_matrices()` queries locations for `projection`, `view` every draw call.
  - Actions:
    - Build a cache map per program: `{uniform_name: loc}` at compile time (`OGLSystem.compile_all_shaders()`), or lazily the first time, then reuse.
    - For FS passes, precompute a small struct per pass: `program`, `out_fbo_id`, `input_units`, `uniform_locs`.
  - Risk: invalidate caches on shader rebuilds (`GLCLWidget.rebuild_gl_resources()`) and when FS config changes (`set_fs_config()`).

- __Avoid redundant state changes in FS pipeline__
  - Today each pass toggles depth/scissor and calls `glViewport` even if unchanged.
  - Actions:
    - Track last depth/scissor state and viewport; only change when needed.
    - Group state changes around the entire FS pipeline (depth/scissor off once; restore once).
  - Risk: Ensure default framebuffer viewport is restored correctly for QOpenGLWidget.

- __Reduce Python name/dict lookups in hot loops__
  - Actions:
    - Pre-resolve: in `bake_kernels()`, store plain tuples `(kernel, global_size, local_size, args)` in a list to iterate (keep `BakedKernelCall` as a lightweight container but access attributes into locals before loop to reduce attribute lookups).
    - For FS pipeline, pre-resolve texture names → texture IDs, program name → program int. At frame time only bind IDs and draw.
  - Risk: refresh these caches on reload/rebuild.

- __Batch OpenCL queue waits__
  - Today `OCLSystem.execute_kernel()` calls `.wait()` per kernel; this serializes and incurs Python↔driver transitions.
  - Actions:
    - Enqueue all kernels without waiting; call `self.queue.finish()` once per frame (unless debugging). Optionally use returned events and wait on the last one only.
  - Risk: in debug modes it’s often useful to fail fast; make behavior configurable (respect `bDebugCL`).


## Medium-impact optimizations (consider after above)

- __OpenCL↔OpenGL interop to avoid CPU copies__
  - Today: CL→CPU→GL path for VBO updates.
  - Actions:
    - Create GL VBOs first (`GLobject`), then wrap them as CL buffers via `cl.GLBuffer` (requires shared context). In frame loop, `clEnqueueAcquireGLObjects` → run kernels → `clEnqueueReleaseGLObjects`. GL draws directly; eliminate `fromGPU()` and `upload_vbo()` per-frame.
  - Risk: platform-specific context creation and synchronization; must rework `OCLSystem.select_device()` to use an interoperable context and manage acquire/release ordering.

- __Persistent-mapped buffers or PBOs for uploads__
  - If interop is not feasible, use persistent mapping on GL side or PBOs to overlap uploads. Benefit is smaller than full interop but still reduces stalls.
  - Risk: complexity vs. gain; driver-dependent.

- __Minimize `glUseProgram`/`glBindVertexArray` churn__
  - Sort draw calls by program and VAO to reduce state flips (if multiple passes with different buffers exist).
  - Risk: changes to render order may affect blending or later effects; only reorder when safe.


## Simplifications and code organization

- __Refactor `_build_ui()` to BaseGUI patterns__
  - Today: manual group/box setup with a mix of `QGroupBox` + `QVBoxLayout` + `BaseGUI.button()`/`comboBox()`.
  - Actions:
    - Add small helpers in `BaseGUI` (future): `vgroup(title, build_fn)` that creates a `QGroupBox` + `QVBoxLayout` and calls `build_fn(layout)`; returns group for `control_layout.addWidget(group)`. This removes repeated boilerplate in `GLCLBrowser._build_ui()`.
    - Use existing helpers consistently (already applied for buttons and combo box). Potentially add `form_group(title)` for parameter forms.
  - Risk: None, but keep compatibility with `params_layout` being a `QFormLayout` used by `populate_params_from_dict()`.

- __Move `BakedKernelCall` to `OCLsystem.py`__
  - Rationale: conceptually belongs with CL execution; avoids scattering CL abstractions in `GLCLBrowser`.
  - Actions: relocate class and import in `GLCLBrowser`; no behavior change.
  - Risk: circular imports if `OCLsystem.py` references `GLCLBrowser` (it does not); keep it independent.

- __Centralize expression resolution__
  - `_resolve_expression()` currently in `GLCLBrowser`. Make it a small module-level util in the same file or move to a tiny `utils` section, to reuse for textures and buffer sizes.
  - Risk: very low; keep behavior identical.


## Detailed targets by file

- `python/GLCL2/GLCLBrowser.py`
  - __`update_simulation()`__
    - Preallocate `host_data` per buffer; reuse arrays.
    - Optionally switch to `self.ocl_system.queue.finish()` once per frame after enqueues.
  - __`bake_kernels()`__
    - Store lightweight tuples and iterate without extra method indirection; or keep class but hoist attributes to locals before loop.
  - __`_precompute_buffer_sync_list()`__
    - OK; ensure list order stable to allow prealloc arrays to match iteration order.
  - __`_build_ui()`__
    - Replace repetitive group setup with `BaseGUI` helper(s) when added; current usage of `button()`/`comboBox()` is already good.

- `python/GLCL2/GLCLGUI.py`
  - __`update_matrices()`__
    - Cache uniform locations per program (projection/view) in a dict; avoid `glGetUniformLocation` every draw.
  - __`_execute_fs_pipeline()`__
    - Pre-bake pass descriptors: program id, out FBO id, input texture ids, uniform locations, pass size (w,h), flags indicating which uniforms exist.
    - Keep a small state cache for current viewport and depth/scissor flags.
  - __`_draw_display_texture()`__
    - Minor: cache uniform location `uTex` after program creation.

- `python/GLCL2/OGLsystem.py`
  - __`GLobject.upload_vbo()`__
    - Track `self._size_bytes`; use `glBufferSubData` when size unchanged.
  - __Shader management__
    - Optional: cache uniform locations when programs are compiled; expose `get_uniform_location(program, name)` that uses a dict lookup instead of GL query.
  - __Texture binding__
    - Optional: bind by texture id if caller can pre-resolve IDs; retain current API for clarity but add a fast path.

- `python/GLCL2/OCLsystem.py`
  - __Queue usage__
    - Provide `execute_many(kernels)` which enqueues a list of `(kernel, global, local, args)` and calls `queue.finish()` once (configurable).
  - __House `BakedKernelCall`__
    - Move the class here; optionally replace by a simple `dataclass` or `namedtuple` for lower overhead.


## Potential pitfalls and what can break

- __Uniform location caching lifetime__
  - After `GLCLWidget.rebuild_gl_resources()` or any shader recompile, cached locations are invalid. Ensure caches are dropped on rebuild and re-populated.

- __Viewport/state restoration with QOpenGLWidget__
  - QOpenGLWidget uses an internal FBO; do not leave a custom FBO bound. Always restore previous framebuffer (already handled in `render_fs_to_fbo`). Be careful when caching viewport; restore to the value from `resizeGL()` or the captured `GL_VIEWPORT`.

- __OpenCL/GL interop context creation__
  - Creating a CL context that shares with the current GL context is platform dependent (WGL/GLX/CGL/EGL). `select_device()` will need an alternate path and error handling. Interop introduces `acquire/release` ordering requirements and error modes.

- __Host buffer preallocation on parameter changes__
  - If a parameter (e.g., `particle_count`) changes buffer sizes, arrays must be reallocated. Hook into `reset_simulation()` and rebake to rebuild `host_sync_buffers` and GL VBOs.

- __Batching kernel waits__
  - Debugging becomes harder if all kernels are enqueued and only a final wait is performed. Keep a debug mode that waits per-kernel as today (guard by `bDebugCL`).

- __GL state caching invalidation__
  - If other code changes GL state out-of-band, a simplistic cache can desync. Encapsulate state changes inside `OGLSystem` and let it own the cache; avoid modifying state elsewhere.


## Refactor plan (incremental, safe steps)

1) Preallocate host sync arrays and switch `GLobject.upload_vbo` to use `glBufferSubData` when possible. Add a tiny uniform cache for `projection`/`view`. Measure.

2) Cache FS pipeline metadata (programs, FBO ids, uniform locations, pass sizes). Remove per-frame `glGetUniformLocation` calls and redundant state flips. Measure.

3) Batch OpenCL kernel waits to a single `queue.finish()` per frame (when not in `--debug-cl` mode). Measure.

4) Optional: introduce GL state caching in `OGLSystem` for `glUseProgram`, `glBindFramebuffer`, `glViewport`.

5) Consider CL/GL interop to eliminate CPU copies for render-bound buffers. Prototype on one buffer to validate reliability on the target system(s).


## _build_ui() consolidation (no code changes yet)

- Introduce in `BaseGUI`:
  - `def vgroup(self, title, builder):` creates `QGroupBox/QVBoxLayout`, calls `builder(layout)`, returns group.
  - `def form_group(self, title):` returns a `QGroupBox` with a `QFormLayout` assigned and bound to `self.params_layout`.
- Re-express `GLCLBrowser._build_ui()` via these helpers to remove boilerplate while keeping current control structure.


## Moving `BakedKernelCall`

- Move class to `OCLsystem.py` to keep CL abstractions localized. `GLCLBrowser` will import it from there; behavior unchanged.
- Alternatively, drop the class altogether and store plain tuples in `self.baked_kernel_calls`:
  - `(kernel_obj, global, local, args)`; in `update_simulation()` do a simple for-loop invoking kernels directly. This removes one Python method call per kernel per frame.


## Measurement guidance

- Add frame-time logging around:
  - Kernel enqueues vs. total CL time (`queue.finish()` scoped).
  - Buffer sync path (CL→CPU→GL upload).
  - FS pipeline time and per-pass time.
- Use a simple macro toggle to avoid logging in release runs.


## Implemented: Dynamic fullscreen uniform handling

- Default uniforms are handled explicitly:
  - `iResolution` is set once per program in `GLCLWidget.apply_fs_static_uniforms()` based on the current FS target size.
  - `iFrame` is updated per frame in `GLCLWidget._execute_fs_pipeline()`.
- All other FS uniforms are applied dynamically from the user script configuration:
  - Declared uniforms are read from `config["opengl_shaders"][shader_name][2]` and matched against `config["parameters"]` by name and type.
  - Supported scalar/vector types: `int`, `float`, `vec2`, `vec3`, `vec4`.
  - Sampler aliases (`iChannelN`, `texN`, `sN`) are bound from the pass inputs; they are not expected in parameters.
  - Hardcoded special-cases (e.g., `driver`) were removed; behavior now fully agnostic to uniform names.
- Uniform locations are cached per program to avoid repeated `glGetUniformLocation` calls.
- Uniforms are (re)applied on GL init and whenever GUI parameters change, minimizing per-frame work.

## Summary

Focus first on memory allocation and GL/CL query overhead in hot loops (preallocate arrays, cache uniforms/IDs, avoid redundant state changes). These are safe, contained changes with immediate payoff. Interop is the long-term win for eliminating CPU copies; do it once the simpler steps are stable.
