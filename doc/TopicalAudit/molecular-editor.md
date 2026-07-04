---
type: TopicalAudit
title: Molecular Editor
tags: [topic, cpp, python, opencl, molecular, gui, editor, mmff, forcefield]
---

## Summary

Interactive molecular editor and simulator with 3D visualization, force field evaluation, and molecular dynamics. Three C++ variants exist: `MolecularEditor` (CPU), `MolecularEditor2` (unused/commented out in build), and `MolecularEditorOCL` (OpenCL-accelerated). A Python OpenCL wrapper (`MolecularDynamics.py`) provides headless MD simulation. The editor supports atom/bond manipulation, MMFF force field, geometry optimization (FIRE/DynamicOpt), and real-time dynamics.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/apps/MolecularEditor/` | active | CPU molecular editor app. Uses `MolecularEngine` OBJECT library, `DynamicOpt` for optimization. Built via `cpp/apps/CMakeLists.txt`. |
| C++ | `cpp/apps_OCL/MolecularEditorOCL/MolecularEditorOCL_main.cpp` | active | OpenCL-accelerated editor. `AppMolecularEditorOCL` extends `AppSDL2OGL_3D`. Manages MMFF force field, OpenCL buffers, interactive editing, CPU/GPU stepping. |
| C++ | `cpp/apps_OCL/MolecularEditorOCL/MolecularEditorOCL_scanner.cpp` | active | Scanner variant for force field scanning. |
| C++ | `cpp/apps/MolecularEditor2/` | inactive | Commented out in `cpp/apps/CMakeLists.txt`. Second generation editor, not built. |
| Python | `python/pyMolecular/OCL/MolecularDynamics.py` | active | OpenCL MD wrapper. Extends `OpenCLBase`. Manages buffers, kernels, force scanning, MD steps. Supports MMFF and raw atom modes. |
| Python | `python/pyMolecular/OCL/MMFF.py` | active | MMFF force field parameter assignment: bond parameters (UFF/simple), neighbors, electron pairs, pi orbitals. |
| Doc | `docs/MolGUI_web.md` | doc | Design document for web-based molecular GUI. |
| Doc | `cpp/apps/MolecularEditor/doc/GlobalOptimization.md` | doc | Documentation on global optimization strategies. |

## Sub-topics

### C++ Molecular Editor (CPU)

- App class derived from `AppSDL2OGL_3D` with SDL2/OpenGL rendering
- Uses `MolecularEngine` OBJECT library (shared between CPU and OCL variants)
- `DynamicOpt` for geometry optimization (FIRE algorithm)
- Interactive atom picking, bond creation/breaking, molecule loading from files

### MolecularEditorOCL (OpenCL-accelerated)

`AppMolecularEditorOCL` in `MolecularEditorOCL_main.cpp`:
- Extends `AppSDL2OGL_3D` for 3D camera and input
- OpenCL buffers for atom positions, forces, neighbor lists
- MMFF force field evaluation on GPU
- CPU fallback stepping (`cpuStep`) for debugging
- Real-time visualization of forces, bonds, molecular geometry
- Scanner variant for systematic force field parameter scanning
- CMake links against `MolecularEngine`, `DynamicOpt`, `SDL2OGL` object libraries + OpenCL/OpenGL/SDL2

### Python OpenCL Molecular Dynamics

`MolecularDynamics.py` extends `OpenCLBase`:
- Buffer management: positions, velocities, forces, neighbor grids
- Kernel execution for force calculation and integration
- Two initialization modes: MMFF molecular system or raw atom positions
- Grid force field with B-spline interpolation
- Texture-based grid force field sampling
- Methods: `scanForces()`, `runMD()`, `initGridFF()`, `uploadAtoms()`, `downloadForces()`

### MMFF Force Field Parameters

`MMFF.py` handles:
- `MMFFsp3_loc` representation conversion
- Bond parameter assignment via UFF and simple models
- Neighbor detection, electron pair assignment, pi orbital handling
- Force constant calculations and reallocation of FF arrays

## Parity Status

- **C++ MolecularEditor ↔ MolecularEditorOCL**: Both use `MolecularEngine` and `DynamicOpt` object libraries. OCL variant offloads force evaluation to GPU. Core molecular data structures shared.
- **C++ MolecularEditorOCL ↔ Python MolecularDynamics.py**: Both use OpenCL for force evaluation. Python version is headless (no GUI), C++ version is interactive. Both support MMFF force field. Buffer layouts and kernel interfaces should match but no formal parity test exists.
- **MMFF.py ↔ C++ MMFF.h**: Python `MMFF.py` assigns parameters that are consumed by OpenCL kernels. C++ `MMFF.h` defines the force field data structures. Parameter format compatibility is implicit, not verified by tests.

## Open Issues

- `MolecularEditor2` is commented out in build — status unclear, may be incomplete or abandoned
- No formal parity tests between C++ and Python OpenCL MD implementations
- `MolecularEditorOCL_scanner.cpp` purpose and usage not well documented
- Web-based molecular GUI (`MolGUI_web.md`) is design-only, no implementation found
- OpenCL kernel source files referenced via symlinks in CMake (`cl/` directory) — fragile path setup
- No automated tests for molecular editor functionality
