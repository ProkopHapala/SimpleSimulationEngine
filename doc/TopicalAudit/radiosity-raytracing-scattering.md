---
type: TopicalAudit
title: Radiosity, Ray Tracing & Radiation Scattering
tags: [topic, cross-language, radiosity, ray-tracing, scattering, radiation, gpu]
---

## Summary

Multiple approaches to radiative heat transport, radiosity, and particle scattering (X-ray, neutron). C++ implementations: `TriangleRayTracer.h` (occlusion), `Radiosity.h` (matrix-based radiosity), `Scatterer.h` (matrix-based scattering), `Scatterer2.h` (channel-network scattering). Python: 2D/3D radiosity solvers in `python/Radiosity/`, OpenCL X-ray attenuation in `python/pyScatter/`. Shared OpenCL kernel `radiosity.cl` for GPU occlusion. Primary application: spacecraft design — radiation from engines, weapons, stars; structural thermal loading.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/TriangleRayTracer.h` | active | Triangle sampling to `SurfElement`s, sequential occlusion test via `rayInTriangle()`. No spatial acceleration. Base class for Radiosity and Scatterer. |
| C++ | `cpp/common/dynamics/Radiosity.h` | active | Dense coupling matrix `M[i][j]` with Lambertian cosine/r² + occlusion. `step_Direct` iterative update. `processTriangles` TODO. Toy/demo quality. |
| C++ | `cpp/common/dynamics/Scatterer.h` | unfinished | Matrix-based scattering with direction bins, `ScatterElem`. Many TODOs, unclear normalization. Largely abandoned. |
| C++ | `cpp/common/dynamics/Scatterer2.h` | unfinished | Channel-network approach: `Channel`/`HalfChannel` between `ScatterElem2`s. `step_Direct` has scattering + transfer phases. **Bug**: uses undefined `n` (should be `scattelems.size()`), `backChans` not populated. |
| C++ | `cpp/sketches_SDL/3D/test_Radiosity.cpp` | active | Demo: load mesh, sample elements, build coupling, iterate, visualize panels + coupling lines |
| C++ | `cpp/sketches_SDL/3D/test_Scatterer.cpp` | active | Demo for Scattering2: two triangles, flux sources, channels, visualization |
| C++ | `cpp/sketches_SDL/3D/test_RayScattererMMC.cpp` | experimental | Ray-based Monte Carlo scatterer test |
| OpenCL | `cpp/common_resources/cl/radiosity.cl` | active | GPU occlusion kernel: brute-force + broad-phase (group-level) + narrow-phase design. Local memory tiling, barrier sync. Used by Python `TriangleOcclusionRaytracer.py`. |
| Python | `python/Radiosity/Radiosity.py` | active | 2D radiosity: line elements, view factors, solve (I-W)B=P, T⁴ formulation. No occlusion. |
| Python | `python/Radiosity/Radiosity3D.py` | active | 3D radiosity: hex-rasterize mesh faces to elements, 3D view factors, solve, visualize. No occlusion. |
| Python | `python/Radiosity/TriangleOcclusionRaytracer.py` | active | PyOpenCL occlusion matrix between element COGs using `radiosity.cl`. Separate from Radiosity3D — integration TODO. |
| Python | `python/Radiosity/PolygonRasterization.py` | active | Hex rasterization of 2D polygons: clipping, centroid/area, boundary merging. Used by Radiosity3D. |
| Python | `python/Radiosity/OctahedralSphereMaping.py` | active | Octahedral sphere parameterization for direction sampling |
| Python | `python/pyScatter/xray_sim.py` | active | OpenCL X-ray attenuation: Beer-Lambert through clustered tubes/triangles. Deterministic, no MC scatter. |
| Python | `python/pyScatter/cl/SpacecraftSimulation.cl` | active | OpenCL kernel for spacecraft radiation simulation |
| Python | `python/pyScatter/cl/cl/Scatter.cl` | unfinished | Placeholder kernel — no real scattering logic, needs full rewrite |
| Python | `python/pyScatter/AttenuationTest.py` | active | Attenuation test driver |
| Python | `python/pyScatter/SimulationPipeline.py` | active | Pipeline orchestration for pyScatter |
| Python | `python/pyScatter/SpacecraftGeometry.py` | active | Spacecraft geometry generation for scattering tests |
| Doc | `doc/RayTracing.md` | active | CPU vs OpenCL implementation comparison for radiosity coupling |
| Doc | `docs/RadiosityRaytracing/ScatteringSolverComparison.md` | active | Detailed comparison: `Scatterer.h` (matrix) vs `Scatterer2.h` (channel) |
| Doc | `docs/RadiosityRaytracing/RadiosityAndScatteringOverviewDiscussion.md` | active | LLM chat: full review of all scattering/radiosity files, status, gaps |
| Doc | `docs/RadiosityRaytracing/TriangleOcclusionRaytracer.md` | active | Design plan: GPU occlusion, broad-phase/narrow-phase hierarchy |
| Doc | `docs/RadiosityRaytracing/SphereSampling.md` | active | Icosahedral diamond mapping for S² sampling (used for direction discretization) |
| Doc | `python/Radiosity/Radiosity.md` | active | 2D radiosity tutorial/design |
| Doc | `python/Radiosity/Radiosity3D.md` | active | 3D radiosity design note |
| Doc | `python/Radiosity/OctahedralSphereMaping.md` | active | Octahedral mapping design |
| Doc | `python/Radiosity/PolygonRasterization.md` | active | Hex rasterization design |
| Doc | `python/pyScatter/MonteCarloScattere.md` | active | MC scattering design dialogue (ideas only, no code) |

## Parity Status

- **C++ `Radiosity.h` ↔ Python `Radiosity3D.py`**: Same algorithm (coupling matrix + iterative solve). C++ has occlusion via `TriangleRayTracer`, Python does not (separate `TriangleOcclusionRaytracer.py` not integrated).
- **C++ `Scatterer.h` ↔ `Scatterer2.h`**: Different approaches to same problem. `Scatterer2.h` is newer but buggy. `Scatterer.h` is abandoned.
- **OpenCL `radiosity.cl` ↔ C++ `TriangleRayTracer.h`**: GPU vs CPU occlusion. GPU uses local memory tiling + broad/narrow phase. CPU is brute-force sequential. No automated parity test.
- **Python `xray_sim.py` ↔ C++**: Python-only, no C++ counterpart for X-ray attenuation.

## Open Issues

- `Scatterer2.h` `step_Direct` uses undefined `n` — should be `scattelems.size()`
- `Scatterer2.h` `backChans` allocated but never populated — backward channel mapping missing
- `Scatterer.h` largely unfinished, many TODOs — candidate for deprecation
- `Scatter.cl` is a non-functional placeholder
- Python `Radiosity3D.py` lacks occlusion — `TriangleOcclusionRaytracer.py` exists but is separate
- `Radiosity.h` `processTriangles` is commented out (TODO)
- `Xray.cl` assumes `cluster_size ≤ 64` (CACHE_CAPACITY), no caps on tubes
- No cross-language parity tests for any radiosity/scattering result
- Spacecraft application: radiation from engines/stars/weapons uses radiosity+scattering tools with different source terms (see `encyclopedia/space_warfare/07_Physical_model_and_equations.md`)
