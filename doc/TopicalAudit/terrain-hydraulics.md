---
type: TopicalAudit
title: Terrain & Hydraulic Erosion
tags: [topic, cross-language, terrain, erosion, hydraulic, gpu]
---

## Summary

Hydraulic erosion simulation for procedural terrain generation. C++ implementation in `TerrainHydraulics.h` with simplex noise terrain. Python OpenCL GPU version in `terrain_ocl/`. Extensive JavaScript web version in `LandCraft_web/` with ShaderToy-inspired erosion shaders. Multiple documentation files covering 1D hydraulics, basin filling, terrain cubic rendering.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/maps/TerrainHydraulics.h` | active | CPU hydraulic erosion on terrain heightmap |
| C++ | `cpp/common/maps/` (TerrainSimplex) | active | Simplex noise terrain generation |
| C++ | `cpp/sketches_SDL/2D/test_TerrainHydraulics.cpp` | active | SDL2 demo of hydraulic erosion |
| Python | `python/terrain_ocl/` | active | GPU OpenCL terrain erosion |
| Python | `python/terrain_ocl copy/` | deprecated | Backup copy — exclude |
| JS | `js/LandCraft_web/` | active | WebGL terrain with erosion shaders, ShaderToy inspiration |
| JS | `js/LandCraft_web_bak/` | deprecated | Backup — exclude |
| Doc | `docs/LandCraft/TerrainHydraulics.md` | active | Hydraulic erosion algorithm doc |
| Doc | `docs/LandCraft/hydraulics1D.md` | active | 1D hydraulic model |
| Doc | `docs/LandCraft/BasinFilling.md` | active | Basin filling algorithm (48KB) |
| Doc | `docs/LandCraft/TerrainCubic.md` | active | Cubic terrain rendering |
| Doc | `python/terrain_ocl/terrain_ocl.md` | active | Python OCL terrain doc (LLM chat) |
| ShaderToy | `js/LandCraft_web/ShaderToy_inspiration/` | reference | Downloaded erosion/noise/terrain shaders — not primary source |

## Parity Status

- **C++ `TerrainHydraulics.h` ↔ Python `terrain_ocl/`**: CPU vs GPU implementation of same erosion algorithm. No automated parity test found.
- **C++ ↔ JS `LandCraft_web/`**: Web version is independent implementation, likely diverged algorithm.

## Open Issues

- `terrain_ocl copy/` directory is a backup — should be excluded/cleaned
- `LandCraft_web_bak/` is a backup — should be excluded
- Accumulative Signed Distance Field erosion: `python/terrain_ocl/Accumulative_Signed_Distance_Field_Errosion.md` (low quality)
- Watershed transform: `python/terrain_ocl/Watershed_Transform.md` (low quality)
- No cross-language parity test for erosion results
- See `doc/Markdown/cpp/common/maps/TerrainSimplex.md` for auto-generated API doc
