---
type: TopicalAudit
title: Orbital Mechanics & Trajectory Optimization
tags: [topic, cross-language, orbital, trajectory, optimization, spline, n-body, ode]
---

## Summary

Orbital mechanics with Keplerian elements, N-body gravitational simulation, and trajectory optimization via spline-based variational methods. C++ `SpaceBodies.h` defines `OrbitalElements`, `Orbit`, `SpaceBody` with Keplerian propagation and astorb database loading. `SpaceWorld.h` implements N-body integration via RKF45. Trajectory optimization uses cubic B-spline control points optimized by FIRE (`DynamicOpt.h`) with `TrajectoryVariation.h` objective. `SplineManager.h` for non-uniform hermite spline evaluation. Application: `OrbitalWar`/`SpaceTactics` game. Python Jupyter notebooks for combat calculations.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/Orbital/SpaceBodies.h` | active | `OrbitalElements` (semi-major, eccentricity, inclination, arg_periapsis, node_longitude, mean_anomaly). `Orbit` with vis-viva equation, true anomaly solver, `pointAtEpoch()`, `toPoints()`. `SpaceBody` with trajectory arrays (`trjPos`, `trjThrust`), hermite spline interpolation. `SpaceBodyIntegrator` for single-body RKF45. `fromString_astorb()` for asteroid database loading. |
| C++ | `cpp/common/Orbital/SpaceWorld.h` | active | N-body simulation: `getDerivODE()` computes gravitational forces + ship thrust. `predictTrjs()` integrates all bodies. `ODEintegrator_RKF45` with adaptive timestep. Dual representation: Keplerian orbits for planets + N-body for ships. `load_astorb()` for real asteroid data. |
| C++ | `cpp/common/dynamics/ODEintegrator.h` | active | `ODEintegrator` base (Euler) + `ODEintegrator_RKF45` (Runge-Kutta-Fehlberg 4(5) with adaptive step). Supports both function pointer and virtual object derivatives. |
| C++ | `cpp/common/math/spline_hermite.h` | active | Hermite spline basis functions: `basis()`, `dbasis()`, `ddbasis()`, `curve_point()` |
| C++ | `cpp/common/math/SplineManager.h` | active | Non-uniform hermite spline manager: control points + optional explicit derivatives. `eval()`, `evalUniform()`, `insertPoint()`, `removePoint()`. Binary search for time index. |
| C++ | `cpp/common/math/CubicBSpline.h` | active | Cubic B-spline basis: `val()`, `dval()`, `ddval()`, `basis()`, `dbasis()`, `ddbasis()` |
| C++ | `cpp/common/Orbital/TrajectoryVariation.h` | active | `getTrajectoryVariation()` — objective function for trajectory optimization: minimizes thrust integral along spline control points |
| C++ | `cpp/common/dynamics/DynamicOpt.h` | active | FIRE (Fast Inertial Relaxation Engine) optimizer: `bindArrays()`, `move_FIRE()`. Used for trajectory optimization. |
| C++ | `cpp/apps/OrbitalWar/test_OptContinuousThrust.cpp` | active | Trajectory optimization demo: 20 B-spline control points, FIRE optimizer, `getTrajectoryVariation` objective. Visualizes control points + forces. |
| C++ | `cpp/apps/OrbitalWar/spaceTactics.cpp` | active | Space combat game: N-body physics, trajectory splines, time scrubbing, weapon tests |
| C++ | `cpp/common/Orbital/spaceCombat.h` | active | Combat assembly: guns, projectiles, hit detection, damage evaluation |
| C++ | `cpp/common/Orbital/asteroidEngineering.h` | active | Asteroid composition (mineral/rock types), orbital maneuver cost (`manuever_planeChange`) |
| C++ | `cpp/common/dynamics/appliedPhysics.h` | active | `centralGravityForce()`, physical constants (G, AU) |
| Python | `projects/SpaceCombat/ch1_basics.ipynb` | active | Orbital mechanics basics notebook |
| Python | `projects/SpaceCombat/Manueverability.ipynb` | active | Spacecraft maneuverability analysis |
| Python | `projects/SpaceCombat/Weapons.ipynb` | active | Weapon performance calculations |
| Python | `projects/SpaceCombat/MISO_paper_draft.ipynb` | active | MISO (Multi-Input Single-Output) paper draft for space combat |
| Python | `projects/SpaceCombat/Industry.ipynb` | active | Space industry/economy model |
| Python | `projects/SpaceCombat/py/` | active | Python support modules |
| Doc | `docs/SpaceTactics/SpaceTactics.md` | active | Comprehensive code map: N-body, weapons, propulsion, shields |
| Doc | `docs/SpaceTactics/SpaceBodies.md` | active | Celestial body definitions |
| Doc | `docs/SpaceTactics/SpaceWorld.md` | active | World model documentation |
| Doc | `docs/SpaceTactics/SplineManager.md` | active | Spline trajectory management |
| Doc | `docs/SpaceTactics/appliedPhysics.md` | active | Applied physics models |
| Doc | `docs/SpaceTactics/asteroidEngineering.md` | active | Asteroid mining/engineering |
| Doc | `docs/SpaceTactics/spaceCombat.md` | active | Combat model documentation |

## Sub-topics

### Keplerian Orbital Mechanics

- `OrbitalElements` → `Orbit` conversion via Euler angles (`fromEuler_orb`)
- True anomaly from mean anomaly: iterative Kepler equation solver (`true_anomaly()`)
- Vis-viva equation: `v = sqrt(GM*(2/r - 1/a))`
- Orbital point at epoch: `pointAtEpoch()` → 3D position from true anomaly + orbit rotation
- Epoch shifting: `shift_epoch()` updates mean anomaly
- Real asteroid loading from astorb database format

### N-Body Gravitational Simulation

- `SpaceWorld::getDerivODE()` — gravitational forces between all bodies + ship thrust
- RKF45 integrator with adaptive timestep (SAFETY=0.2, PGROW=1.2, PSHRINK=0.7)
- Trajectory storage: `trjPos[]` arrays per body, hermite spline interpolation
- Dual representation: planets on Keplerian orbits, ships on N-body trajectories
- `predictTrjs()` — integrate forward, store positions at each timestep

### Trajectory Optimization

- **B-spline control points** represent trajectory
- **FIRE optimizer** (`DynamicOpt.h`) minimizes objective
- **`getTrajectoryVariation()`** — thrust integral objective (minimize total thrust)
- `test_OptContinuousThrust.cpp` — demo with 20 control points, real-time visualization
- Commented-out spline editor code suggests interactive trajectory editing was planned

## Parity Status

- **C++ N-body ↔ Python notebooks**: Notebooks perform analytical calculations (delta-v, maneuverability) but don't simulate N-body dynamics. No direct parity.
- **Keplerian orbits ↔ N-body integration**: Dual representation in `SpaceWorld`. Planets use analytical Keplerian, ships use numerical integration. No consistency check found.
- **`test_OptContinuousThrust` ↔ `SpaceTactics` game**: Optimization demo is standalone; game uses spline trajectories but optimization not integrated into gameplay.

## Open Issues

- `SplineManager.insertPoint()` is a placeholder (`return -1`)
- `SpaceBodyIntegrator` thrust interpolation commented out (`//f.add()` in `getDerivODE`)
- `test_OptContinuousThrust.cpp` has commented-out interactive editing (mouse handling, event handling)
- No automated test for Keplerian ↔ N-body consistency
- `appliedPhysics.h` location inconsistent — included from `cpp/common/dynamics/` but referenced as `cpp/common/Orbital/`
- Trajectory optimization uses FIRE but no comparison with other optimizers (CG, L-BFGS)
- No multi-body trajectory optimization (only single spacecraft)
- `projects/SpaceCombat/` notebooks are standalone calculations, not integrated with C++ simulation
- See `doc/Markdown/cpp/common/Orbital/` for auto-generated API docs
