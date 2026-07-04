---
type: TopicalAudit
title: PyOpenCL Workflows
tags: [topic, python, opencl, pyopencl, buffer-management, kernel-compilation, glcl, cl-gl-sharing]
---

## Summary

Python OpenCL infrastructure spanning two wrapper generations: `OpenCLBase` (pyMolecular, with auto kernel header parsing and arg generation) and `OCLSystem` (GLCL2, with kernel caching and `BakedKernelCall`). GLCL Browser provides GUI for loading simulation scripts, editing parameters, and visualizing results via OpenGL. CL↔GL interop for zero-copy rendering. Eulerian fluid simulation uses direct PyOpenCL without wrapper class.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| Python | `python/pyMolecular/OCL/OpenCLBase.py` | active | `OpenCLBase`: context/queue init, `load_program()`, `extract_kernel_headers()` (parse `__kernel` signatures from source), `create_buffer()`, `toGPU()`/`fromGPU()`, `generate_kernel_args()` (auto-bind args from header + `kernel_params` dict), `parse_cl_lib()` (extract `//>>>function` / `//>>>macro` sections). ~540 lines. |
| Python | `python/pyMolecular/OCL/clUtils.py` | active | `get_nvidia_device()`, `roundup_global_size()`, `get_cl_info()` (device info, local memory characteristics), `try_load_clFFT()` (gpyfft), `make_inds_pbc()` (PBC index patterns) |
| Python | `python/GLCL2/OCLsystem.py` | active | `OCLSystem`: context/queue, multi-program loading (`load_program(name, path)` with kernel caching via `prg.all_kernels()`), buffer dict with auto-reallocation, `toGPU()`/`fromGPU()`, `execute_kernel()` / `enqueue_kernel()` (batched), `BakedKernelCall` (pre-baked args + work sizes). ~194 lines. |
| Python | `python/GLCL2/OGLsystem.py` | active | `OGLSystem`: OpenGL context, shader compilation, texture/FBO management for CL↔GL interop |
| Python | `python/GLCL2/GLCLBrowser.py` | active | `GLCLBrowser(BaseGUI)`: script loading, parameter widgets, simulation control, CL→GL buffer sync, display texture selection. ~678 lines. |
| Python | `python/GLCL2/GLCLGUI.py` | active | `GLCLWidget`: QOpenGLWidget with CL↔GL sharing, rendering pipeline |
| Python | `python/GLCL2/ModularGL.py` | active | Modular OpenGL rendering utilities |
| Python | `python/GLCL/NBody_glcl.py` | active | Simple N-body: inline OpenCL source, PyQt5 OpenGL widget, CL buffer (no GL sharing — fallback to `cl.Buffer` copy). See `parallel-particle-cell.md`. |
| Python | `python/GLCL/NBody_glsl.py` | active | N-body with GLSL compute (alternative to OpenCL) |
| Python | `python/EulerianImpacFluid/EulerianImpacFluid.py` | active | Eulerian fluid: direct PyOpenCL (no wrapper class), buffer management, kernel dispatch. See `fluid-dynamics.md`. |
| Python | `python/pyMolecular/OCL/MolecularDynamics.py` | active | `MolecularDynamics(OpenCLBase)`: multi-system MD. See `parallel-particle-cell.md`. |
| Doc | `docs/OpenCLBase.md` | doc | OpenCL base class documentation |
| Doc | `docs/MolGUI_web.md` | doc | Molecular editor web version |

## Sub-topics

### OpenCLBase (pyMolecular)

Key features:
- **Kernel header parsing**: `extract_kernel_headers(source)` parses `__kernel void name(args)` from source, stores in `kernelheaders` dict
- **Auto arg generation**: `generate_kernel_args(kernel_name)` matches header parameters to `kernel_params` dict and `buffer_dict` by name — eliminates manual `clSetKernelArg` calls
- **CL library parsing**: `parse_cl_lib(path)` extracts `//>>>function NAME` and `//>>>macro NAME` sections for code injection/preprocessing
- **Device selection**: `select_device(preferred_vendor='nvidia')` with fallback to `cl.create_some_context()`

### OCLSystem (GLCL2)

Key features:
- **Multi-program support**: `load_program(name, path)` — multiple `.cl` files compiled and cached by name
- **Kernel caching**: `prg.all_kernels()` → global `kernels` dict by function name
- **Buffer auto-reallocation**: `create_buffer(name, size, flags)` — releases and recreates if size changed
- **Batched execution**: `enqueue_kernel()` (no `.wait()`) for batched per-frame finish
- **BakedKernelCall**: Pre-baked kernel call with frozen args + work sizes for hot loops

### CL↔GL Interop

- `GLCLBrowser`: CL computes simulation → CL buffer downloaded to host → uploaded to GL VBO for rendering
- `NBody_glcl.py`: Attempted GL sharing but falls back to `cl.Buffer` copy (`use_gl_sharing = False`)
- `GLCLWidget`: QOpenGLWidget with CL context sharing
- `OGLSystem`: Texture/FBO management for display

### Script-Driven Simulation (GLCLBrowser)

- Load Python script → script defines config (buffers, kernels, parameters, render pipeline)
- `baked_kernel_calls` list: pre-baked kernel calls executed per frame
- `buffers_to_sync`: CL→GL buffer copy list per frame
- Parameter widgets: auto-generated from config, editable at runtime
- Display texture selection: multiple output textures from CL kernels

## Parity Status

- **`OpenCLBase` ↔ `OCLSystem`**: Two parallel implementations with overlapping functionality. `OpenCLBase` has auto-arg generation (kernel header parsing). `OCLSystem` has multi-program support and `BakedKernelCall`. No convergence — both maintained independently.
- **`NBody_glcl.py` (GLCL) ↔ `GLCLBrowser` (GLCL2)**: GLCL is standalone demo, GLCL2 is framework. Different approaches to CL↔GL.
- **Eulerian fluid (direct PyOpenCL) ↔ wrapper classes**: Uses neither `OpenCLBase` nor `OCLSystem` — direct `cl.Buffer` and `cl.Program` calls.

## Open Issues

- Two parallel OpenCL wrapper hierarchies (`OpenCLBase` vs `OCLSystem`) — DRY violation, should converge
- `NBody_glcl.py` CL↔GL sharing not working — falls back to host copy
- `OpenCLBase::generate_kernel_args()` not fully audited — complex name-matching logic
- `OCLSystem::create_buffer()` silently allows `size <= 0` (sets buffer to `None`)
- No error handling for kernel execution failures in `BakedKernelCall` beyond try/catch + print
- `clUtils::get_nvidia_device()` has f-string bug: `"No {what} device found"` (missing `f` prefix)
- `GLCLBrowser` is ~678 lines — large monolithic class
- No CL↔GL interop in `OpenCLBase` — only `OCLSystem`/`OGLSystem` path supports it
- Eulerian fluid uses direct PyOpenCL — inconsistent with wrapper pattern
