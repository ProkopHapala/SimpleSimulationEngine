# Spacecraft Simulator TODO

**Generated from**: `docs/SpaceCrafting/Codebase_Reference_Map.md` and codebase analysis.
**Last updated**: 2025-07-18

---

## Priority Legend

- **P0** — Blocks all other work. Fix first.
- **P1** — Core feature gap. Needed for minimum viable simulation.
- **P2** — Important for realism / gameplay. Can be deferred.
- **P3** — Enhancement / polish.

---

## P0: Bond Breaking (Structural Damage)

### Problem
`evalBondTension()` in `TrussDynamics_d.cpp:1788` has a TODO: "break the bond if strain > maxStrain". `maxStrain` is populated per bond from material properties (`Spull/Kpull`, `Spush/Kpush`) in `exportSim()` but never checked.

### Tasks
- [ ] Implement bond breaking: when `strain[i] > maxStrain[i].x` (pull) or `< maxStrain[i].y` (push), remove bond from neighbor lists (`neighs`, `neighBs`, `neighB2s`) and mark as broken.
- [ ] Handle solver topology change: Cholesky decomposition (`LDLT_L`, `LDLT_D`) is invalidated when bonds break. Either:
  - Switch to iterative solver (Jacobi/CG) when damage occurs, or
  - Re-decompose Cholesky periodically (expensive for large meshes).
  - **Reference**: Python `pyTruss/truss.cl` has 12 OpenCL kernels for iterative solvers (Jacobi, GS, VBD, Momentum) that handle topology changes naturally. Consider porting to C++ OpenCL path.
- [ ] Test: Create a simple truss, apply extreme load, verify bonds break and structure separates.
- [ ] Verify `maxStrain` values are sensible for real materials (steel, carbon fiber).

### Files
- `cpp/common/dynamics/TrussDynamics_d.h` — `maxStrain`, `strain`, `bonds`, `neighs`
- `cpp/common/dynamics/TrussDynamics_d.cpp` — `evalBondTension()` (line 1775), `linsolve()` (line 1030)
- `cpp/common/Orbital/SpaceCraft2Mesh2.h` — `exportSim()` populates `maxStrain` (line 137)
- `python/pyTruss/truss.cl` — OpenCL iterative solver kernels (reference for topology-robust solvers)
- `python/pyTruss/truss_solver_ocl_new.py` — GPU solver with VBD/Jacobi (reference)

---

## P0: Self-Collision Detection

### Problem
After bonds break (ship disintegration), truss fragments can pass through each other. `Buckets` spatial hashing (`pointBBs`, `edgeBBs`, `faceBBs`) exists but only `evalTrussCollisionImpulses_bonds()` is implemented, which handles only bonded point pairs.

### Existing Python Implementation (Reference)
- **`python/particleCollisions2D/constrains.cl`** — OpenCL kernels for unified bond+collision neighbor system. `updateJacobi_neighs` handles both bonds and collisions in one Jacobi pass (collision pairs flagged with `l0 < 0`). `updateCollisonNeighbors` does broad-phase with spatial groups + local memory. This approach could be ported to C++ `TrussDynamics_d`.
- **`python/particleCollisions2D/constrains_cl.py`** — `CLConstrains` class wrapping the OpenCL kernels.

### Tasks
- [ ] Implement broad-phase collision: use `pointBBs`/`edgeBBs` to find candidate colliding pairs between disconnected fragments. Consider porting `updateCollisonNeighbors` from Python.
- [ ] Implement narrow-phase: point-edge and edge-edge collision response (impulse-based, similar to `evalTrussCollisionImpulses_bonds`).
- [ ] Consider unified bond+collision neighbor approach from `constrains.cl` — same `neighs`/`kngs`/`l0s` arrays, `l0 < 0` sentinel for collisions.
- [ ] Handle point-face collision using `faceBBs` (for plate/shield penetration).
- [ ] Test: Break a truss into two pieces, push them toward each other, verify collision response.

### Files
- `cpp/common/dynamics/TrussDynamics_d.h` — `pointBBs`, `edgeBBs`, `faceBBs`, `collision_damping`
- `cpp/common/dynamics/TrussDynamics_d.cpp` — `evalTrussCollisionImpulses_bonds()` (line 1740)
- `python/particleCollisions2D/constrains.cl` — reference implementation (OpenCL)
- `python/particleCollisions2D/constrains_cl.py` — Python wrapper

---

## P1: Projectile Damage Model

### Problem
No code computes which truss edges/nodes are hit by a projectile. `spaceCombat.h` has `whippleShieldType` and `ProjectedTarget` but these are high-level abstractions not connected to the truss mesh.

### Tasks
- [ ] Implement ray-truss intersection: cast ray (projectile trajectory) through `TrussDynamics` points/edges/faces. Use `TriangleRayTracer` for face hits.
- [ ] Probabilistic hit model: for a projectile cloud (spread of parallel rays), compute probability of hitting each girder/edge based on cross-sectional area and occlusion.
- [ ] Apply damage: break bonds near hit point, transfer kinetic energy as impulse to nearby nodes.
- [ ] Nuclear/radiative damage: use radiosity framework to distribute energy from a flash point across surface elements.
- [ ] Test: Fire a projectile at a simple truss, verify it breaks bonds near impact.

### Files
- `cpp/common/Orbital/spaceCombat.h` — `ProjectileType`, `SpaceGunType`, `whippleShieldType`, `ProjectedTarget`
- `cpp/common/dynamics/TriangleRayTracer.h` — `getOcclusion()`, ray-triangle test
- `cpp/common/dynamics/Radiosity.h` — surface element coupling for radiative damage

---

## P1: Pipe / Line Attachment System

### Problem
`Pipe` class exists in `SpaceCraftComponents.h:648` with `maxFlow`, `Path`, and component pointers `a`/`b`, but:
- No mesh generation (pipes don't appear in `BuildCraft_blocks`).
- No simulation integration (no flow dynamics, no thermal coupling).
- `component_kind()` returns `Rope` (placeholder).

### Tasks
- [ ] Define pipe types: fuel pipe, coolant pipe, power cable, heat pipe. Each with different flow models.
- [ ] Mesh generation: pipe follows a `Path` along girders (like slider rails). Generate thin edges representing the pipe/cable.
- [ ] Attachment points: define input/output ports on components (tanks, thrusters, radiators). Pipe connects ports.
- [ ] Flow simulation: simple 1D flow network (graph with nodes = junctions, edges = pipes). Solve for pressure/flow rate.
- [ ] Thermal coupling: heat pipes transfer thermal energy between connected components.
- [ ] Power grid: cables transfer electrical power from generators (thrusters, reactors) to consumers (guns, sliders).

### Files
- `cpp/common/Orbital/SpaceCraftComponents.h` — `Pipe` class (line 648)
- `cpp/common/Orbital/SpaceCraft.h` — `std::vector<Pipe*> pipes` (line 50)
- `cpp/common/Orbital/EditSpaceCraft.h` — add `l_Pipe` Lua function

---

## P1: Slider Path Initialization (C++)

### Problem
C++ `Slider::updatePath()` extracts path vertices from a `StructuralComponent` side. `findNearestPoint()` is commented out. `StructuralComponent::updateSlidersPaths()` iterates attached nodes but path initialization may be incomplete.

### Tasks
- [ ] Uncomment and test `Slider::findNearestPoint()` — needed to snap slider to closest path position.
- [ ] Verify `updateSlidersPaths()` correctly initializes all sliders on a ring with multiple girders.
- [ ] Port JS `pAttach` + `weldDist` mechanism to C++ for more robust slider-to-girder attachment.
- [ ] Port JS cubic B-spline interpolation (`bsplineInterpolate`) to C++ `Path::getPos()`.
- [ ] Test: Create a ring with 4 sliders on girders, verify paths are correct and sliders move smoothly.

### Files
- `cpp/common/Orbital/SpaceCraftComponents.h` — `Slider::updatePath()` (line 804), `findNearestPoint()` (line 815, commented)
- `js/common_js/SplineCubic.js` — `bsplineInterpolate()` (reference for porting)
- `js/spacecraft_editor/js/SpaceCraft2Mesh.js` — `pAttach` mechanism (line 237)

---

## P1: Actuator / Control System

### Problem
`Slider::icontrol` maps to a `control_speed[]` array in `applySliders2sim()`, but there is no higher-level control system. No user input mapping, no PID controllers, no actuator sequencing.

### Tasks
- [ ] Define `SpaceCraftControl` class: array of DOF indices, each with target value, current value, and control mode (position/velocity/force).
- [ ] User input mapping: keyboard/joystick → DOF targets.
- [ ] PID controller per DOF: compute `control_speed` from target vs. current position.
- [ ] Actuator limits: enforce `forceMax`, `powerMax`, `maxSpeed` from `Slider` class.
- [ ] Multi-slider coordination: group sliders (e.g., 4 sliders on a ring) to share a single DOF.
- [ ] Test: Drive a ring with 4 sliders via keyboard, verify coordinated movement.

### Files
- `cpp/common/Orbital/SpaceCraftComponents.h` — `Slider::icontrol`, `forceMax`, `powerMax`, `maxSpeed`
- `cpp/common/Orbital/SpaceCraft2Mesh2.h` — `applySliders2sim()` (line 261)

---

## P1: Port Parabolic Mesh Generators to C++

### Problem
`ParabolaSheet`, `ParabolaSlab_wrap`, `ParametricParabolaPatch` exist only in JS (`MeshesUV.js`, `MeshGenerators.js`). C++ `Mesh::Builder2` lacks these. The nozzle/damper prototype only works in JS.

### Tasks
- [ ] Port `ParabolaUVfunc` to C++ (simple: `z = f(r)` with `R0`, `R`, `L` parameters).
- [ ] Port `ParabolaSheet` to `Mesh::Builder2` (single-layer truss following parabolic surface).
- [ ] Port `ParabolaSlab_wrap` to `Mesh::Builder2` (double-layer with bracing between layers).
- [ ] Port `ParametricParabolaPatch` to `Mesh::Builder2` (annular patch with independent inner/outer tessellation).
- [ ] Port `TubeSheetBond` (two tubes with bond-length analysis and welding).
- [ ] Add C++ test in `constructionBlockTests` equivalent.
- [ ] Test: Generate nozzle mesh in C++, verify same vertex/edge count as JS.

### Files
- `js/spacecraft_editor/js/MeshesUV.js` — `ParabolaUVfunc`, `ParabolaSheet`, `ParabolaSlab_wrap`
- `js/spacecraft_editor/js/MeshGenerators.js` — `ParametricParabolaPatch`, `TubeSheetBond`
- `cpp/common/geometry/MeshBuilder2.h` — target for porting

---

## P2: Integrate MHD with Spacecraft Thruster Model

### Problem
MHD plasma nozzle simulation (`js/mhd_demo/physics.js`, `doc/python/MHD/inductance_core.py`) is standalone. `Thruster` class in `SpaceCraftComponents.h` has no magnetic field modeling.

### Tasks
- [ ] Define thruster types: chemical (simple impulse), nuclear pulse (Orion-style pusher plate + damper), magnetic nozzle (MHD plasma).
- [ ] For magnetic nozzle: link `Thruster` to MHD simulation parameters (coil geometry, plasma properties).
- [ ] Compute thrust from MHD: integrate magnetic pressure on cage coils → net axial force.
- [ ] Couple MHD force to `TrussDynamics` as external force on nozzle attachment points.
- [ ] For nuclear pulse: model pusher plate impulse + damper compression (slider travel).
- [ ] Test: Attach magnetic nozzle to simple truss, verify thrust force propagates through structure.

### Files
- `cpp/common/Orbital/SpaceCraftComponents.h` — `Thruster` class
- `js/mhd_demo/physics.js` — `stepSimulation()`, `coilField()`, `sampleFieldAtPoint()`
- `doc/python/MHD/inductance_core.py` — `build_inductance_matrix()`, `calc_mutual_inductance()`
- `doc/python/MHD/doc/MHD_Plasma_Nozzle_Simulation.md` — design document

---

## P2: Radiosity Integration with Spacecraft Simulation

### Problem
`Radiosity` class is included in `EditSpaceCraft.h` but not wired into the simulation loop. Surface elements from spacecraft mesh are not fed to the radiosity solver during dynamics.

### Existing Python Implementation (Reference)
- **`python/Radiosity/Radiosity3D.py`** — Full 3D radiosity solver: face rasterization (`rasterize_face_to_elems`), view factors (`compute_view_factors_3d`), linear solve (`solve_radiosity_system_3d`), temperature calculation (`calculate_temperatures_3d`). Can be used as reference for C++ integration.
- **`python/Radiosity/TriangleOcclusionRaytracer.py`** — PyOpenCL wrapper for GPU occlusion matrix using shared `radiosity.cl` kernel.
- **`python/pyScatter/SimulationPipeline.py`** — Full pipeline: geometry → rasterize → occlusion → radiosity → temperatures. Uses `SpacecraftSimulation.cl` kernel.
- **`python/pyScatter/SpacecraftSimulation.cl`** — OpenCL kernels: `compute_occlusion` (ray-triangle), `compute_attenuation` (X-ray optical depth).

### Tasks
- [ ] Extract surface triangles from `TrussDynamics` faces (or from `SpaceCraft` plates/shields).
- [ ] Run `Radiosity::processTriangles()` to build coupling matrix at startup.
- [ ] During simulation: update element positions from `TrussDynamics` points, recompute coupling (or approximate with rigid-body transform).
- [ ] Solve radiosity for thermal equilibrium: `M * vals = sources` (using existing `LinSolver`).
- [ ] Feed thermal results back: temperature affects material properties (stiffness, strength).
- [ ] GPU acceleration: use `radiosity.cl` kernels for occlusion computation. Python `TriangleOcclusionRaytracer.py` already wraps this kernel — use as reference.
- [ ] Validate C++ results against Python `Radiosity3D.py` solver on same geometry.
- [ ] Test: Heat one side of a spacecraft, verify temperature distribution across surfaces.

### Files
- `cpp/common/dynamics/Radiosity.h` — `makeCouplingMatrix()`, `processTriangles()`, `step_Direct()`
- `cpp/common/dynamics/TriangleRayTracer.h` — `getOcclusion()`, `addTriangle()`
- `cpp/common_resources/cl/radiosity.cl` — `occlusion_matrix`, `coupling_matrix` kernels
- `python/Radiosity/Radiosity3D.py` — 3D radiosity solver (reference for validation)
- `python/Radiosity/TriangleOcclusionRaytracer.py` — PyOpenCL occlusion wrapper
- `python/pyScatter/SimulationPipeline.py` — full radiosity pipeline (reference)
- `python/pyScatter/SpacecraftSimulation.cl` — OpenCL occlusion + attenuation kernels
- `cpp/common/Orbital/EditSpaceCraft.h` — `radiositySolver` instance (line 27)

---

## P2: JS → C++ Building System Port (Evaluation)

### Problem
JS spacecraft editor has more advanced mesh generators and UI, but no dynamics. C++ has dynamics but older mesh generation. User asked to evaluate porting JS building system back to C++.

### Pros
- C++ has OpenGL rendering (`SpaceCraftDraw.h`) for visual debugging.
- C++ has `TrussDynamics` for immediate simulation feedback.
- C++ Lua scripting is mature for design automation.
- Better performance for large meshes.

### Cons
- JS editor is actively developed with modern UI (web-based, collaborative potential).
- Porting all JS generators back to C++ is significant effort.
- JS MeshBuilder is already a faithful port of C++ `Mesh::Builder2` — most operations exist in both.

### Recommendation
- **Don't wholesale port JS editor to C++.** Instead:
  1. Port specific missing generators (parabolic dishes, `TubeSheetBond`) to C++ `Mesh::Builder2`.
  2. Keep JS editor as the primary interactive design tool.
  3. Use mesh export CLI (`spaceCraft_mesh_export_cli.md`) to transfer designs from JS to C++ simulation.
  4. Add C++ visual test harness for mesh generators (simple OpenGL viewer with parameter sliders).

### Tasks
- [ ] Port `ParabolaSheet`, `ParabolaSlab_wrap`, `ParametricParabolaPatch` (see P1 above).
- [ ] Port `TubeSheetBond` with bond-length analysis.
- [ ] Port JS `pAttach` slider attachment to C++.
- [ ] Port cubic B-spline interpolation to C++ `Path`.
- [ ] Add C++ mesh generator test GUI (parameter sliders + OpenGL preview).
- [ ] Verify mesh export CLI can transfer JS designs to C++ simulation format.

---

## P2: Collision Avoidance in Telescopic Damper

### Problem
Telescopic truss recoil damper (hexagonal outer + triangular inner tube) requires collision avoidance between inner and outer tube during compression. No collision avoidance code exists for this specific case.

### Tasks
- [ ] Define collision shapes for tube cross-sections (hexagonal, triangular).
- [ ] Implement per-segment collision check: for each slider position, verify inner tube vertices don't penetrate outer tube walls.
- [ ] Add collision response: if penetration detected, apply impulse to push inner tube back.
- [ ] Test: Compress telescopic damper, verify no interpenetration.

### Files
- `js/spacecraft_editor/js/constructionBlockTests.js` — "TubeSheetBond Hex Ring" test (line 340)
- `docs/SpaceCrafting/SpaceCraftConstructionProblems.md` — §1.1-1.5 damper design
- `cpp/common/dynamics/TrussDynamics_d.h` — collision infrastructure

---

## P3: JS Editor Pending Tasks

From `SpaceCraftConstructionProblems.md`:

- [ ] **[J2] Visualizers for index/angle/SDF rail selection strategies** — show which vertices are selected by each method for debugging.
- [ ] **[J7] Advanced SDF-based selection regions** — n-fold angular quadrants, axial slabs for more complex selection patterns.

---

## P3: SpaceCraft2Truss Integration Gaps

### Problem
`SpaceCraft2Truss.h` has placeholders for shields, radiators, thrusters (lines 88-92) — these components are not converted to truss elements.

### Tasks
- [ ] Convert `Shield` plates to truss faces (triangulated quads with appropriate stiffness).
- [ ] Convert `Radiator` plates to truss faces with thermal properties.
- [ ] Convert `Thruster` to truss node with external force application point.
- [ ] Convert `Tank` to truss with mass distribution (partially done — `toTruss` has optional tank handling).
- [ ] Test: Build a spacecraft with shields + radiators + thrusters, verify all appear in truss.

### Files
- `cpp/common/Orbital/SpaceCraft2Truss.h` — `toTruss()` (line 50-92)

---

## P2: Port Python OpenCL Truss Solvers to C++

### Problem
Python `pyTruss/truss.cl` (849 lines, 12 kernels) has the most complete OpenCL truss solver suite: VBD (chunk + serial), Jacobi (sparse, local, fly, diff), Gauss-Seidel (colored, fly, diff), momentum mixing, precomputed RHS. C++ `spaceCraftDynamicsOCL.cpp` exists but is less documented and may not have all solver variants. Graph coloring (Vivace/Luby's MIS) is implemented in Python but not in C++.

### Tasks
- [ ] Audit C++ `spaceCraftDynamicsOCL.cpp` / `spaceCraftDynamicsOCL.h` for existing OpenCL kernels.
- [ ] Compare kernel set: which Python kernels are missing in C++?
- [ ] Port missing kernels (likely `jacobi_diff`, `GS_diff`, `GS_fly`, `precompute_dRHS`, `pd_momentum_mix`).
- [ ] Port graph coloring (`color_graph` from `pyTruss/truss.py`) to C++ for parallel GS.
- [ ] Port Chebyshev acceleration parameters from `pyTruss/IterativeLinearSolvers.py`.
- [ ] Validate: Run same truss problem on both Python and C++ OpenCL, compare results.

### Files
- `python/pyTruss/truss.cl` — 12 OpenCL kernels (849 lines, reference)
- `python/pyTruss/truss_solver_ocl_new.py` — GPU solver wrapper (reference)
- `python/pyTruss/truss.py` — `color_graph()` (Vivace graph coloring)
- `python/pyTruss/IterativeLinearSolvers.py` — Chebyshev/momentum acceleration
- `python/pyTruss/sparse.py` — sparse operations, `gauss_seidel_iteration_colored()`
- `cpp/common/dynamics/TrussDynamics_d.h` — C++ solver enum (14 methods)
- `cpp/apps/OrbitalWar/spaceCraftDynamicsOCL.cpp` — C++ OpenCL dynamics (audit target)
