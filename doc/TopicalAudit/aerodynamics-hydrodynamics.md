---
type: TopicalAudit
title: Aerodynamics & Hydrodynamics Simulation
tags: [topic, cross-language, aerodynamics, hydrodynamics, fluid, potential-flow, vortex]
---

## Summary

Aerodynamic modeling for aircraft and spacecraft: panel-based aero surfaces with polar stall model (`AeroSurf.h`), rigid-body aircraft with propellers (`AeroCraft.h`), vortex lattice method via Biot-Savart law (`PotentialFlow.h`), and combat aerodynamics model (`AirCombatModel.h`). Potential flow and vortex particle methods documented for web implementations. Applications: `AeroCombat` (3D flight sim), `SailWar` (2D sail ships), `SpaceTactics` (spacecraft in atmosphere).

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/dynamics/AeroSurf.h` | active | `AeroSurface` with `AeroPolar` model: CD0, dCD, dCL, stall transition via `trashold_cub`. `applyForce()` computes lift/drag from angle of attack. Supports both polar model and coefficient model. |
| C++ | `cpp/common/dynamics/AeroSurf2D.h` | active | 2D aerodynamic surface with compact polar model (CD0, dCD, dCL, sStall, wStall). Fewer parameters than 3D version. |
| C++ | `cpp/common/dynamics/AeroCraftControl.h` | active | `SurfControl` wraps a surface with angle limits; `AeroCraftControler` coordinates ailerons/rudder/elevator for altitude hold and roll stabilization. |
| C++ | `cpp/common/dynamics/AeroCraft.h` / `.cpp` | active | `AeroCraft : RigidBody` with panels + propellers. `Propeler` class: thrust from momentum theory (quadratic solve for Δv). `applyAeroForces()` iterates panels + propellers. File-based aircraft definition. |
| C++ | `cpp/common/math/PotentialFlow.h` | active | Biot-Savart law: `dBiotSawart()`, `ILineFinite()`, `ILineSemiInf()`, `ILineSemiInfDecay()` (Lorenzian vortex core). Horseshoe vortex: `horseshoe()`, `horseshoeDecay()`. Semi-infinite vortex sheet: `ISemiInfSheet()`. Source dipole: `sourceDipol()`. |
| C++ | `cpp/sketches_SDL/3D/test_VortexLattice.cpp` | active | Vortex lattice method demo: lift-line wing, horseshoe vortices with decay, velocity field visualization, numerical integration verification |
| C++ | `cpp/common/CombatModels/AirCombatModel.h` | active | Analytical aerodynamic performance: `maxSpeed_simple()`, `trunRate_simple()`, `climbRate_simple()`, `climbRate_CLmax()`. Drag/lift/thrust/gravity. Slip-away tactics. Energy management. |
| C++ | `cpp/apps/AeroCombat/` | active | 3D subsonic aircraft combat simulator using `AeroCraft` |
| C++ | `cpp/apps/SailWar/` | active | 2D sail-ship simulator (wind aerodynamics on sails) |
| JS | `docs/BiotSavart/AeroPotentialFlow_simulation_javascript.md` | active | Web-based potential flow simulation design doc |
| JS | `docs/BiotSavart/VortexParticleMethod.md` | active | Vortex particle method design (25KB, detailed) |
| JS | `docs/BiotSavart/MHD_Plasma_Nozzle_webgl.md` | active | MHD plasma nozzle simulation (webgl) |
| Doc | `docs/BiotSavart/PotentialFlow.md` | active | Potential flow theory and implementation notes |
| Doc | `docs/AeroCombat/AeroCombat.md` | active | AeroCombat app documentation |
| Doc | `docs/AeroCombat/AeroSurf.md` | active | Aerodynamic surface model documentation |
| Doc | `docs/LandTactics/AirCombatModel.md` | active | Air combat model: aerodynamics, turn performance, tactical maneuvering |
| Doc | `doc/Markdown/cpp/common/dynamics/AeroSurf.h.md` | active | Auto-generated API doc |
| Doc | `doc/Markdown/cpp/common/dynamics/AeroCraft.h.md` | active | Auto-generated API doc |

## Sub-topics

### AeroSurface Polar Model

Stall transition using smoothed threshold function:
- Pre-stall: `CD = CD0 + dCD*|sin(α)|*|sin(α)|`, `CL = dCL*sin(α)`
- Post-stall: `CD = CD0 + dCDS*|sin(α)|*|sin(α)|`, `CL = dCLS*cos(α)*sin(α)`
- Transition via `trashold_cub(|sin(α)|, sStall, sStall+wStall)` — cubic blend

### Vortex Lattice Method (PotentialFlow.h)

- Biot-Savart law for vortex filament induced velocity
- Horseshoe vortex: bound vortex + two trailing semi-infinite vortices
- Lorenzian decay core (`ILineSemiInfDecay`) for regularization
- Semi-infinite vortex sheet integration for lifting surface theory

### Propeller Model

Momentum theory: `dm = ρ*S*(v0+vstatic)`, solve quadratic `0.5*dm*Δv² + dm*v0*Δv - P = 0` for Δv, then `thrust = dm*Δv*efficiency - dm*v0*CD`

## Parity Status

- **C++ `AeroSurf.h` polar model ↔ `AirCombatModel.h` analytical**: Different abstraction levels. `AeroSurf` is per-panel force computation, `AirCombatModel` is whole-aircraft performance envelope. No direct parity test.
- **C++ `PotentialFlow.h` ↔ JS VortexParticleMethod**: C++ implements vortex lattice (horseshoe), JS design doc describes vortex particle method. Different discretization of same physics.
- **C++ `AeroCraft` ↔ `AeroCombat` app**: App uses `AeroCraft` directly. Integration tested through gameplay.

## Open Issues

- `AeroSurface.applyForceSimple()` is commented out — replaced by polar model version
- `AeroPolar` struct duplicated in both `AeroSurf.h` and as member variables in `AeroSurface` — DRY violation
- Vortex lattice method (`test_VortexLattice.cpp`) is a standalone demo, not integrated into `AeroCraft`
- No hydrodynamics (water) implementation found — only aerodynamics
- `SailWar` aerodynamics not audited in detail
- MHD plasma nozzle (`MHD_Plasma_Nozzle_webgl.md`) is web-only, no C++ counterpart
- Potential flow docs in `docs/BiotSavart/` are design docs / LLM chats, not implementation references

## Related Audits

- **`fluid-dynamics.md`** — Overview of all fluid dynamics implementations. References this audit for potential flow / vortex methods.
- **`continuum-mechanics-impact.md`** — Compressible Eulerian multi-material solver (level set, 5-equation model). Different regime: compressible flow with shocks vs. incompressible potential flow here. Also contains a comprehensive **Algorithm Review: Heterogeneous Material Simulation Methods** covering all multi-material algorithms.
- **`parallel-particle-cell.md`** — Grid infrastructure and neighbor search. Vortex particle methods (designed in JS docs) would use this infrastructure for O(N) neighbor finding.
- **`soft-body-truss-dynamics.md`** — Rigid body dynamics (`RigidBody` base class used by `AeroCraft`).
