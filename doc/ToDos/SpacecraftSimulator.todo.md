# Spacecraft Simulator TODO

**Generated from**: `docs/SpaceCrafting/Codebase_Reference_Map.md` and codebase analysis.
**Last updated**: 2025-07-05

---

## Priority Legend

- **P0** — Blocks all other work. Fix first.
- **P1** — Core feature gap. Needed for minimum viable simulation.
- **P2** — Important for realism / gameplay. Can be deferred.
- **P3** — Enhancement / polish.

---

## Completed: SpaceCraftEditorNew Crash & Visualization Fixes (2025-07-05)

### Problem
`SpaceCraftEditorNew` crashed on startup with `ship_ICF_marksman_2.lua` due to 4 separate `exit(0)` calls triggered by unconnected vertices and uninitialized edge types. Additionally, thick lines and broken 3D text rendering made visualization unusable.

### Fixes Applied

**Crash fixes (4 `exit(0)` calls eliminated):**
- [x] `checkAllPointsConnected(true, true)` → `checkAllPointsConnected(false, true)` in `SpaceCraft2Mesh2.h` — 12 unconnected bound nodes no longer cause exit.
- [x] `o->edge_type = 0` in `l_Node` (`EditSpaceCraft.h`) — uninitialized `edge_type=-1` caused `exportSim` to exit on invalid stick material index.
- [x] `exit(0)` on bad edge type in `exportSim` (both double and float versions) → `continue` (skip bad edges with warning).
- [x] `checkMasses()` → `checkMasses(1e-9, false)` in `SpaceCraft2Mesh2.h` — zero-mass unconnected points no longer cause exit.
- [x] Assign `mass=1.0` to unconnected points in both `exportSim` (double) and `exportSim_f` (float) — prevents NaN in Cholesky decomposition (D[j]=0 → division by zero).
- [x] Bound nodes (`l_BoundNode`, `l_Rope2`) registered in `components` vector for `find_mesh_element` lookup.
- [x] `SpaceCraft::clear()` fixed to delete only via `components` (avoids double-free).

**Visualization fixes:**
- [x] `glLineWidth(1.0)` added at start of `SpaceCraftDynamicsApp::draw()` — prevents stale thick line width from previous frame affecting truss rendering.
- [x] `glLineWidth(1.0)` reset added after picking visualization block in `SpaceCraftEditorNew.cpp` — was set to 5.0 but never reset.
- [x] `glEnable(GL_DEPTH_TEST)` in `SpaceCraftDynamicsApp::draw()` (was `glDisable`) — matches old editor, enables proper 3D occlusion for truss and text rendering.
- [x] `glDisable/glEnable(GL_DEPTH_TEST)` around vertex numbering overlay — 3D text labels render on top of geometry.

### Files Modified
- `cpp/common/Orbital/SpaceCraft2Mesh2.h` — crash fixes (exit→warn, mass fix, checkMasses bExit=false)
- `cpp/common/Orbital/EditSpaceCraft.h` — edge_type init, components registration
- `cpp/common/Orbital/SpaceCraft.h` — clear() double-free fix, find_mesh_element cleanup
- `cpp/apps/OrbitalWar/SpaceCraftDynamicsApp.h` — GL_DEPTH_TEST enable, glLineWidth reset
- `cpp/apps/OrbitalWar/SpaceCraftEditorNew.cpp` — glLineWidth resets, depth test for vertex labels

### Remaining
- [ ] 12 unconnected vertices from bound nodes still exist (non-fatal, warnings only). Root cause: `make_anchor_point` creates vertices before girders are processed, so no edges connect them. Fix would require reordering mesh building or deferring anchor connections.
- [x] Memory leaks fixed — see "Memory Leak Fixes" section below.

---

## Completed: Memory Leak Fixes (2025-07-05)

### Problem
`SpaceCraftEditorNew` had no destructors or deallocation logic. ASan reported **187802 bytes leaked in 591 allocations**. Major sources: `TrussDynamics_d`/`TrussDynamics_f` raw arrays (48MB), Lua binding objects not registered in `components`, global `lua_State` never closed, GUI panels/panels not deleted, `SpaceCraft`/`SpaceCraftSimulator` not destructed.

### Root Causes
1. **No `delete app`** in `main()` — entire app object tree leaked.
2. **No destructors** on classes owning raw `new[]` arrays: `TrussDynamics_d`, `TrussDynamics_f`, `SparseMatrix`, `SparseMatrix2`, `Buckets`, `CGsolver`.
3. **Missing `components.push_back`** in 10 Lua binding functions (`l_Rope2`, `l_Weld`, `l_Radiator`, `l_Shield`, `l_Tank`, `l_Thruster`, `l_Gun`, `l_Rock`, `l_Balloon`, `l_Slider`) — objects allocated with `new` but never registered for deletion by `SpaceCraft::clear()`.
4. **Double `initSpaceCraftingLua()`** — called in both `SpaceCraftDynamicsApp` and `SpaceCraftEditorNew` constructors; second call overwrote `theLua` pointer, leaking the first `lua_State`.
5. **No destructors** on `PTree`, `Dict`, `MultiPanel`, `BoundGUI`, `Path`, `SpaceCraft`, `SpaceCraftSimulator`, `SpaceCraftDynamicsApp`.
6. **`workshop` is a global** (not heap-allocated) — `~SpaceCraft` must not `delete` the borrowed `workshop` pointer.

### Fixes Applied

**Destructors added:**
- `~PTree()` in `Tree.h` — recursively deletes branch pointers.
- `~Dict()` in `containers.h` — deletes all owned `T*` in `vec`.
- `~MultiPanel()` in `GUI.h` — deletes all `GUIPanel*` in `subs`.
- `~BoundGUI()` in `GUI.h` — `delete[] drivers`.
- `~Path()` in `SpaceCraftComponents.h` — frees `ps` if not shared (`sharedFrom==0`).
- `~SpaceCraft()` in `SpaceCraft.h` — calls `clear()` (workshop is borrowed, not deleted).
- `~SpaceCraftSimulator()` in `spaceCraftSimulator.h` — calls `sim.dealloc()` + `sim_f.dealloc()`.
- `~SpaceCraftDynamicsApp()` in `SpaceCraftDynamicsApp.h` — `delete simulator`.
- `~SpaceCraftEditorNew()` in `SpaceCraftEditorNew.cpp` — `lua_close(theLua)` + `delete theSpaceCraft`.

**`dealloc()` methods added (for raw-array classes):**
- `TrussDynamics_d::dealloc()` — frees 25+ raw arrays + nested `SparseMatrix`/`SparseMatrix2`/`Buckets`/`CGsolver`.
- `TrussDynamics_f::dealloc()` — same for float version.
- `SparseMatrix::dealloc()` — frees `nng`, `inds`, `vals`.
- `SparseMatrix2::dealloc()` — frees `nngs`, `i0s`, `vals`, `inds`.
- `CGsolver::dealloc()` — frees `invD`, `res`, `d`, `Ad`, `z`, `alpha`, `beta`, `err`, `err2`.
- `Buckets::dealloc()` — frees `cellNs`, `cellI0s`, `cell2obj`, `obj2cell`.

**Ownership fixes (missing `components.push_back`):**
- `EditSpaceCraft.h` — Added `components.push_back(o)` to `l_Rope2`, `l_Weld`, `l_Radiator`, `l_Shield`, `l_Tank`, `l_Thruster`, `l_Gun`, `l_Rock`, `l_Balloon`, `l_Slider`.
- `SpaceCraft.h` — Added `components.push_back` for Slider nodes in `make_Ring2`.

**Other fixes:**
- `delete app` after `app->loop()` in `SpaceCraftEditorNew.cpp:main()`.
- Guard against double `initSpaceCraftingLua()` — close existing `theLua` before creating new state.
- `virtual` destructor on `SpaceCraftSimulator` (deleted through base pointer in `~SpaceCraftDynamicsApp`).

### Files Modified
- `cpp/apps/OrbitalWar/SpaceCraftEditorNew.cpp` — `delete app`, `~SpaceCraftEditorNew` destructor
- `cpp/apps/OrbitalWar/SpaceCraftDynamicsApp.h` — `~SpaceCraftDynamicsApp` destructor
- `cpp/apps/OrbitalWar/spaceCraftSimulator.h` — `~SpaceCraftSimulator` destructor
- `cpp/common/Orbital/SpaceCraft.h` — `~SpaceCraft` destructor, `components.push_back` in `make_Ring2`
- `cpp/common/Orbital/SpaceCraftComponents.h` — `~Path` destructor
- `cpp/common/Orbital/EditSpaceCraft.h` — `components.push_back` in 10 Lua bindings, double-init guard
- `cpp/common/dataStructures/Tree.h` — `~PTree` destructor
- `cpp/common/dataStructures/Buckets.h` — `dealloc()` method
- `cpp/common/utils/containers.h` — `~Dict` destructor
- `cpp/common/dynamics/TrussDynamics_d.h` — `dealloc()` method
- `cpp/common/dynamics/TrussDynamics_f.h` — `dealloc()` method
- `cpp/common/math/SparseMatrix.h` — `dealloc()` method
- `cpp/common/math/SparseMatrix2.h` — `dealloc()` method
- `cpp/common/math/CG.h` — `dealloc()` method
- `cpp/common_SDL/SDL2OGL/GUI.h` — `~MultiPanel`, `~BoundGUI` destructors

### Result
ASan leak report: **187802B / 591 allocs → 141748B / 196 allocs**. Remaining leaks are exclusively from `iris_dri.so` (Mesa OpenGL driver internals, 9216B direct + ~132KB indirect) — not application code.

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

### Completed
- [x] Port JS cubic B-spline interpolation (`bsplineInterpolate`) to C++ — done in `cpp/common/dynamics/KinematicSolver.h`.
- [x] Kinematic solver (`KinematicSystem`) with LM optimization, slider constraints, `sweepSolve` and `sweepSolveRange`.
- [x] Slider paths derived from actual girder vertex indices (same approach as `Girder::sideToPath()`).
- [x] Kinematic animation demo: triangular girder + parabolic dish sliding along 3 axial rails, exported as `.obj+` multi-frame.
- [x] Rail path visualization (static lines) + slider connection lines in animation frames.

### Tasks (remaining)
- [ ] Uncomment and test `Slider::findNearestPoint()` — needed to snap slider to closest path position.
- [ ] Verify `updateSlidersPaths()` correctly initializes all sliders on a ring with multiple girders.
- [ ] Port JS `pAttach` + `weldDist` mechanism to C++ for more robust slider-to-girder attachment.
- [ ] Integrate `KinematicSolver` with `SpaceCraftComponents.h` `Slider` class (replace ad-hoc path setup).
- [ ] Test: Create a ring with 4 sliders on girders, verify paths are correct and sliders move smoothly.

### Files
- `cpp/common/Orbital/SpaceCraftComponents.h` — `Slider::updatePath()` (line 804), `findNearestPoint()` (line 815, commented)
- `cpp/common/dynamics/KinematicSolver.h` — `KinematicSystem`, `bsplineInterpolate()`, `sweepSolve()`, `sweepSolveRange()`
- `cpp/apps/OrbitalWar/constructionBlockApp.cpp` — `-KinematicTest` demo (damper: dish slides along girder)
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

### Completed
- [x] Port `ParabolaUVfunc` to C++ — done in `cpp/common_SDL/SDL2OGL3/DrawUV.h`.
- [x] Port `ParabolaSheet` to C++ `Builder2` — done in `DrawUV.h`.
- [x] Port `ParametricParabolaPatch` to C++ `Builder2` — done in `DrawUV.h`.
- [x] C++ test in `constructionBlockApp.cpp` — `-KinematicTest` generates parabolic dish + triangular girder.

### Tasks (remaining)
- [ ] Port `ParabolaSlab_wrap` to `Mesh::Builder2` (double-layer with bracing between layers).
- [ ] Port `TubeSheetBond` (two tubes with bond-length analysis and welding).
- [ ] Test: Generate nozzle mesh in C++, verify same vertex/edge count as JS.

### Files
- `cpp/common_SDL/SDL2OGL3/DrawUV.h` — `ParabolaUVfunc`, `ParabolaSheet`, `ParametricParabolaPatch` (ported)
- `js/spacecraft_editor/js/MeshesUV.js` — `ParabolaSheet`, `ParabolaSlab_wrap` (reference for remaining)
- `js/spacecraft_editor/js/MeshGenerators.js` — `TubeSheetBond` (reference for remaining)
- `cpp/common/geometry/MeshBuilder2.h` — target for remaining ports

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

### Completed
- [x] Port `ParabolaSheet`, `ParametricParabolaPatch` to C++ (see P1 above).
- [x] Port cubic B-spline interpolation to C++ (`KinematicSolver.h`).
- [x] `dir2tree` refactored to shared utility `cpp/common/utils/dir2tree.h` — eliminates code duplication across MeshViewer, SpaceCraftGUI, and other apps. Accepts optional file filter callback.
- [x] `MeshViewer` (`cpp/apps/MeshViewer/MeshViewer.h`) established as the primary mesh/animation visualizer:
  - OBJ face parser handles quads/ngons (fan-triangulation), no fake wireframe diagonals.
  - Backface culling toggle (`showBackfaces`), on by default (renders both sides).
  - Wireframe only draws explicit edges from OBJ file (no derived triangle edges).
  - Multi-frame `.obj+` animation playback: `[`/`]` step frames, `Space` play/pause, frame slider with live update (`bCmdOnSlider`).
  - Camera no longer resets during frame navigation — stays fixed while stepping through animation.
  - File browser with `dir2tree` + mesh file filter (`isMeshFile`).
  - Display toggles: solid, wireframe, points, vertex/edge/face numbers, normals.
  - Used to visualize kinematic animation results (damper demo, slider connection lines, rail paths).

### Tasks (remaining)
- [ ] Port `ParabolaSlab_wrap` and `TubeSheetBond` with bond-length analysis.
- [ ] Port JS `pAttach` slider attachment to C++.
- [ ] Add C++ mesh generator test GUI (parameter sliders + OpenGL preview).
- [ ] Verify mesh export CLI can transfer JS designs to C++ simulation format.

---

## P2: Collision Avoidance in Telescopic Damper

### Problem
Telescopic truss recoil damper (hexagonal outer + triangular inner tube) requires collision avoidance between inner and outer tube during compression. No collision avoidance code exists for this specific case.

### Completed
- [x] Kinematic demo of damper mechanism: parabolic dish (pusher plate) slides along triangular girder (spine) using `KinematicSolver`. 3 slider constraints on axial rails derived from girder vertices. Animation exported as `.obj+` and visualized in `MeshViewer`.

### Tasks (remaining)
- [ ] Define collision shapes for tube cross-sections (hexagonal, triangular).
- [ ] Implement per-segment collision check: for each slider position, verify inner tube vertices don't penetrate outer tube walls.
- [ ] Add collision response: if penetration detected, apply impulse to push inner tube back.
- [ ] Test: Compress telescopic damper, verify no interpenetration.

### Files
- `js/spacecraft_editor/js/constructionBlockTests.js` — "TubeSheetBond Hex Ring" test (line 340)
- `docs/SpaceCrafting/SpaceCraftConstructionProblems.md` — §1.1-1.5 damper design
- `cpp/common/dynamics/TrussDynamics_d.h` — collision infrastructure
- `cpp/common/dynamics/KinematicSolver.h` — kinematic solver for slider constraints
- `cpp/apps/OrbitalWar/constructionBlockApp.cpp` — `-KinematicTest` damper demo

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

---

## P2: Integrate Burn1D Nuclear Propulsion Solvers

### Problem
1D nuclear/combustion solvers in `doc/python/Burn1D/` are functional but standalone. They are directly relevant to nuclear-powered spacecraft (Orion pusher plate, nuclear thermal rockets, fusion drives) but not yet connected to the spacecraft simulator.

### Solvers Available
- **`implosion_solver.py`** — Spherical implosion (Numba-accelerated staggered grid, DT fusion, neutron diffusion, material layers). Applicable to: Orion pusher plate impulse, fusion drive compression.
- **`neutron_diffusion_solver.py`** — Multi-group neutron diffusion with burnup/depletion, k_eff. Applicable to: space reactor criticality, radiation shielding design.
- **`euler_tube_solver.py`** — 1D Euler with Rusanov flux, variable geometry. Applicable to: nozzle flow, propellant lines.
- **`pulsejet_solver.py`** — Lagrangian pulsejet with 5-species combustion. Applicable to: pulse propulsion cycle.
- **`wave_tube_solver.py`** — Lagrangian acoustic wave. Applicable: pressure waves in propellant lines.

### Tasks
- [ ] Validate each solver produces correct physics (run old versions for parity, compare plots)
- [ ] Couple implosion solver impulse output to `TrussDynamics` as external force on pusher plate
- [ ] Use neutron diffusion solver for radiation shielding design (couple with `Radiosity` for surface heating)
- [ ] Port critical 1D solvers to C++ or OpenCL for real-time simulation integration
- [ ] Define spacecraft nuclear propulsion component model (reactor, shield, pusher plate, nozzle)
- [ ] Test: Simulate Orion-style nuclear pulse with pusher plate + damper + spacecraft truss

### Files
- `doc/python/Burn1D/implosion_solver.py` — `ImplosionSim` class (Numba staggered grid)
- `doc/python/Burn1D/neutron_diffusion_solver.py` — `NeutronDiffusionSim` class (multi-group diffusion)
- `doc/python/Burn1D/euler_tube_solver.py` — `CompressibleSolver` class (Rusanov flux)
- `doc/python/Burn1D/pulsejet_solver.py` — `PulsejetSolver` class (Lagrangian combustion)
- `doc/python/Burn1D/wave_tube_solver.py` — `WaveTubeSolver` class (acoustic wave)
- `docs/SpaceCrafting/Codebase_Reference_Map.md` — §14 Nuclear Propulsion Solvers
