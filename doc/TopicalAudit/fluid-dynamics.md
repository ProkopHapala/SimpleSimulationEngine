---
type: TopicalAudit
title: Fluid Dynamics
tags: [topic, cross-language, fluid, eulerian, lagrangian, opencl, python, potential-flow]
---

## Summary

Fluid dynamics implementations span Eulerian compressible flow (5-equation multi-material model with level set), potential flow (vortex lattice method, Biot-Savart), and Lagrangian particle methods. The Eulerian impact fluid is the most complete implementation with OpenCL kernels. Potential flow is used for aerodynamics. See also `continuum-mechanics-impact.md` for the Eulerian fluid details and `aerodynamics-hydrodynamics.md` for potential flow.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| Python+OpenCL | `python/EulerianImpacFluid/EulerianImpacFluid.py` | active | Eulerian compressible multi-material flow. 5-equation model (ρ, ρu, ρv, E, α₁). PyOpenCL wrapper: buffer management, kernel dispatch, step(), initialization. See `continuum-mechanics-impact.md` for details. |
| OpenCL | `python/EulerianImpacFluid/EulerianImpacFluid.cl` | active | Kernels: EOS (stiffened gas), Rusanov flux, update, redistance (level set reinitialization). Positivity floors. See `continuum-mechanics-impact.md`. |
| C++ | `cpp/common/math/PotentialFlow.h` | active | Biot-Savart law, vortex line, horseshoe vortex, vortex lattice method. Inline functions for aerodynamic lift calculation. See `aerodynamics-hydrodynamics.md`. |
| C++ | `cpp/sketches_SDL/3D/test_VortexLattice.cpp` | active | CFD demo: vortex lattice method for wing lift calculation with visualization |
| C++ | `cpp/sketches_SDL/3D/test_Electromagnetic.cpp` | active | 2D grid-based electromagnetic solver with Poisson solver (iterative + FFT). `solvePoisson()`, `solvePoissonFFT()`, `dampDens()`. |
| JS | `docs/BiotSavart/AeroPotentialFlow_simulation_javascript.md` | doc | Design doc for JavaScript potential flow simulation |
| JS | `docs/BiotSavart/MHD_Plasma_Nozzle_webgl.md` | doc | MHD plasma nozzle — see `mhd-plasma.md` |
| C++ | `cpp/common/dynamics/MechEuler2D.h` | active | 2D Eulerian multi-material fluid: per-cell molar amounts, temperature, pressure, velocity. Advection + cohesion to limit interface diffusion. |
| C++ | `cpp/common/dynamics/MechGrid2D.h` | active | Multi-material rectangular grid: pressure-to-force conversion and momentum update. Simpler than MechEuler2D (no advection/thermodynamics). |
| C++ | `cpp/common/dynamics/MechMesh2D.h` | active | Visco-elastic simulation on arbitrary triangular meshes. Per-triangle pressure/volume/internal energy. Builder class and edge-swap optimization. |
| C++ | `cpp/common/dynamics/MechPIC2D.h` | active | 2D Particle-In-Cell: mass/velocity on particles, pressure on grid. P2G deposit, EOS, G2P force interpolation, Lagrangian advection. |
| C++ | `cpp/common/dynamics/MechPIC2D_Temperature.h` | active | PIC with temperature: thermodynamic state on both particles and cells. Heat conduction and convection. |
| C++ | `cpp/common/dynamics/CompressiveParticles.h` | active | Compressible sphere particles with adiabatic EOS, overlap forces, pressure-driven radius dynamics. For granular flow and compaction. |
| C++ | `cpp/common/dynamics/Shock1D.h` | active | 1D shock waves with planar/cylindrical/spherical symmetry. Layered materials, adiabatic + isochoric EOS. For implosion, shaped charges, explosions. |
| C++ | `cpp/common/dynamics/SuperSonic2D.h` | active | Oblique shock relations (theta-beta-Mach) with polynomial lookup tables for fast supersonic flow deflection solving. |
| Doc | `doc/Continum_Mechanics_Solver.md` | doc | Continuum mechanics solver design document |

## Sub-topics

### Eulerian Compressible Flow

See `continuum-mechanics-impact.md` for full details:
- 5-equation all-speed model: mass, momentum, energy, volume fraction
- Stiffened gas EOS: `p = (γ-1)(E - ½ρv²) - γ·p∞`
- Rusanov (local Lax-Friedrichs) flux with positivity floors
- Level set for material interface tracking + redistancing
- OpenCL implementation with `update`, `flux`, `eos`, `redistance` kernels

### Potential Flow (Vortex Methods)

See `aerodynamics-hydrodynamics.md` for full details:
- Biot-Savart law for vortex-induced velocity
- Vortex lattice method: horseshoe vortices on wing panels
- Vortex line integration
- Used for aerodynamic lift/drag calculation

### Grid-based Poisson Solver

In `test_Electromagnetic.cpp`:
- Iterative: `solvePoisson(nstep, eps)` — Jacobi-type relaxation with optimal `cR` parameter
- FFT-based: `solvePoissonFFT()` — 2D FFT of density, frequency-domain Green's function (incomplete)
- `potential2force()` — gradient of potential to force field
- `dampDens(bmix)` — density field damping for stability

## Parity Status

- **Eulerian fluid Python↔OpenCL**: Direct PyOpenCL binding. See `continuum-mechanics-impact.md`.
- **Potential flow C++ ↔ JS**: JS version is design document only, not implemented.
- **Poisson iterative ↔ Poisson FFT**: Both in `test_Electromagnetic.cpp` but FFT version incomplete (Green's function commented out).

## Open Issues

- FFT Poisson solver incomplete — Green's function multiplication commented out
- No SPH (Smoothed Particle Hydrodynamics) implementation
- No incompressible Navier-Stokes solver
- No GPU implementation of potential flow / vortex lattice
- `test_Electromagnetic.cpp` is a sketch — not production code
- No coupling between Eulerian fluid and potential flow methods
- See `continuum-mechanics-impact.md` and `aerodynamics-hydrodynamics.md` for additional open issues

## Related Audits

- **`continuum-mechanics-impact.md`** — Detailed deep-dive into the Eulerian compressible multi-material solver (5-equation model, level set, OpenCL kernels). Also contains a comprehensive **Algorithm Review: Heterogeneous Material Simulation Methods** covering level set, VOF, CLSVOF, Ghost Fluid, PIC, FLIP, MPM, SPH, LBM, phase field, and more with pros/cons and implementation status.
- **`parallel-particle-cell.md`** — Grid infrastructure (`CubeGridRuler`), cell lists, neighbor search. Relevant for PIC/SPH methods that use the same grid infrastructure for particle-to-grid deposition.
- **`aerodynamics-hydrodynamics.md`** — Potential flow / vortex methods for aerodynamics. Complementary regime (incompressible, irrotational) to the compressible Eulerian solver.
- **`soft-body-truss-dynamics.md`** — Lagrangian elastic dynamics. Different paradigm: structural mechanics vs. fluid mechanics.
- **`collision-detection.md`** — Broad-phase collision detection (hash grid, sweep-and-prune). Relevant for particle-based fluid methods.
