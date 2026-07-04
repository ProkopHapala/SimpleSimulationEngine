---
type: TopicalAudit
title: OpenCL & GPU Computing
tags: [topic, cross-language, gpu, opencl, parallel]
---

## Summary

OpenCL infrastructure in `cpp/common/OCL/` providing GPU compute wrappers. Python OpenCL via `pyOpenCL` in `python/GLCL2/` (GLCL framework), `python/pyRay/` (ray-marching CSG), `python/terrain_ocl/` (terrain erosion). Design docs for OpenCLBase and API optimization. Key concern: high-precision physics on single-precision GPU hardware.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/OCL/` | active | OpenCL wrapper: context, queue, kernel, buffer management |
| C++ | `cpp/common/OCL/OpenCL_API_opt.md` | active | API optimization notes (doc) |
| C++ | `cpp/sketches_OCL/` | experimental | Standalone OpenCL simulation demos |
| C++ | `cpp/apps_OCL/MolecularEditorOCL/` | experimental | GPU-accelerated molecular editor |
| C++ | `cpp/common/dynamics/TrussDynamics_d.cpp` | active | `run_Cholesky_omp_simd` — CPU SIMD parallel, not GPU but related parallel strategy |
| Python | `python/GLCL2/` | active | GLCL: GLSL+OpenCL integration framework, persistent buffers, kernel caching |
| Python | `python/pyRay/` | active | Real-time CSG ray-marching on GPU via pyOpenCL |
| Python | `python/terrain_ocl/` | active | GPU terrain erosion simulation |
| Python | `python/GLCL2/doc/GLCL_manifest.md` | active | GLCL design manifest |
| Python | `python/GLCL2/doc/GLCL_performance_refactor.md` | active | Performance refactoring notes |
| Doc | `docs/OpenCLBase.md` | active | OpenCL base infrastructure documentation |
| Doc | `doc/High_Precision_Physics_on_singlepoint_GPU_hardware.md` | active | Single vs double precision on GPU |
| Doc | `doc/Parallel_Particle_To_Cell_accumulation.md` | active | Parallel particle-to-cell algorithm |

## Parity Status

- **C++ OCL ↔ Python GLCL**: Different abstraction levels. C++ provides low-level wrapper, Python GLCL adds GLSL integration and kernel caching. No direct numerical parity test.
- **terrain_ocl**: Python OpenCL kernels should have CPU reference in `cpp/common/maps/TerrainHydraulics.h` — parity not verified.

## Open Issues

- Single-precision GPU limitations for physics simulations — documented but no systematic solution
- `cpp/sketches_OCL/` has `NOTES.md` and `MySimulationCodes.md` — status unclear
- OpenCL dev notes: `cpp/common/OCL/OpenCL_dev_notes.md` (low quality)
- GLCL-to-GLCL integration: `python/GLCL2/doc/GLSL_to_GLCL_integration.md` (low quality)
- ctypes bindings workflow: `doc/CodingRules/workflows/python/python_ctypes_workflow.md`
- See skill: `doc/AGENTs/skills/port-to-opencl/SKILL.md`, `doc/AGENTs/skills/gpu-optimize/SKILL.md`, `doc/AGENTs/skills/gpu-debug/SKILL.md`
