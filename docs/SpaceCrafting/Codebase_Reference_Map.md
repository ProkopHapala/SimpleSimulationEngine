# Spacecraft Simulator Codebase Reference Map

**Purpose**: Assemble all code references relevant to the spacecraft simulator, organized by subsystem. This is a navigation document for developers.

---

## 1. Spacecraft Component Model

### C++ (Reference Implementation)
- **`cpp/common/Orbital/SpaceCraft.h`** — Central `SpaceCraft` class: holds vectors of nodes, ropes, girders, rings, guns, sliders, radiators, shields, tanks, pipes, thrusters, balloons, rocks, welds. `make_Rope`, `make_Girder`, `make_Ring`, `make_Ring2` for Lua-driven construction. `build_order` for sequencing.
- **`cpp/common/Orbital/SpaceCraftComponents.h`** — Component hierarchy: `ShipComponent` → `StructuralComponent` (girders, rings, ropes) → `Node`, `Weld`, `Modul` (tanks, balloons, rocks), `Pipe`, `Plate` (radiators, shields), `Thruster`, `Rotor`, `Slider`. `Slider` inherits `Node`, has `Path`, motor physics (`move()`), `icontrol` index for user control. `Path` class with `cur`, `fromCur()`, `fold()`, `findNearestPoint()`.
- **`cpp/common/Orbital/EditSpaceCraft.h`** — Lua scripting API: `l_Material`, `l_Node`, `l_Rope`, `l_Girder`, `l_Ring`, `l_Slider`, `l_Weld`, `l_Radiator`, `l_Shield`, `l_Tank`, `l_Thruster`, `l_Gun`, `l_Rock`, `l_Balloon`. Also includes `Radiosity radiositySolver` and `TriangleRayTracer`.

### JavaScript (Active Development)
- **`js/spacecraft_editor/js/SpaceCraft.js`** — JS port of component model: `Node`, `Girder`, `Rope`, `Plate`, `Path`, `Slider`, `Ring`, `SpaceCraft` classes. `Path.interpolate()` uses `bsplineInterpolate()` (cubic B-spline). `Slider` has `pAttach` for explicit attachment point, `weldDist` for welding tolerance.
- **`js/spacecraft_editor/js/SpaceCraftEngine.js`** — Engine: `rebuildMesh()` calls `BuildCraft_blocks_js`, `processCommands()` for incremental updates, `globalRegistry` for object tracking.
- **`js/spacecraft_editor/js/SpaceCraft2Mesh.js`** — `BuildCraft_blocks_js`: block-based mesh generation (nodes→cubes, girders→bridged trusses, ropes, plates, rings→wheels, sliders→path attachment). `buildPathOnEdge` with SDF and stride methods. `BuildCraft_aux_js` for debug visualization of slider paths.

---

## 2. Mesh Generation

### C++
- **`cpp/common/geometry/MeshBuilder2.h`** — `Mesh::Builder2` class: flat arrays for verts/edges/tris/chunks/strips. `addCube`, `bridgeFacingPolygons`, `rope`, `wheel`, `ParametricQuadPatch`, `vert_weld`, `vert_collapse`, SDF selection. Dynamic mode with soft removal.
- **`cpp/common/Orbital/SpaceCraft2Mesh2.h`** — `exportSim()` converts `Mesh::Builder2` → `TrussDynamics_d/f`. `BuildCraft_truss` and `BuildCraft_blocks` iterate components. `applySliders2sim()` connects slider motor to `EdgeVertBond` in simulation. `sliders2edgeverts()` sets up edge-vert constraints.
- **`cpp/common/Orbital/SpaceCraft2Mesh_blocks.h`** — Block-based conversion: `node_to_mesh`, `girder_to_mesh`, `ring_to_mesh`, `slider_to_mesh`, `rope_to_mesh`. `BuildCraft_blocks` orchestrates. `exportSimToFile` for external simulation.
- **`cpp/common/Orbital/SpaceCraft2Truss.h`** — `toTruss()` converts SpaceCraft → `Truss` (nodes, ropes, girders, rings, tanks). **TODOs**: shields, radiators, thrusters not yet integrated (lines 88-92).

### JavaScript
- **`js/common_js/MeshBuilder.js`** — Full JS port of `Mesh::Builder2` (2017 lines): `vert`, `edge`, `wheel`, `addCube`, `bridgeFacingPolygons`, `rope`, `ParametricQuadPatch`, SDF selection, `vert_weld`, `vert_collapse`, `mapBondLengthsFromVertexSelections`, duplicate detection. Selection banks for editing.
- **`js/spacecraft_editor/js/MeshGenerators.js`** — UV-based mesh generators: `TubeSheet`, `SlabTube`, `ParabolaSheet`, `ParabolaSlab`, `ParametricParabolaPatch`, `TorusSheet`, `QuadSheet`, `QuadSlab`, `TubeSheetBond`.
- **`js/spacecraft_editor/js/MeshesUV.js`** — UV surface functions: `ConeUVfunc`, `UV_slab`, `ParabolaUVfunc`, `slabEdges`. Parametric mapping from UV to 3D.
- **`js/spacecraft_editor/js/MeshGenTestGUI.js`** — Interactive GUI for testing mesh generators with parameter sliders.
- **`js/spacecraft_editor/js/constructionBlockTests.js`** — Test harness including "TubeSheetBond Hex Ring" test with **parabolic pusher-plate / plasma nozzle dish** mounted on telescopic damper (lines 340-436). This is the nozzle/damper prototype.
- **`js/common_js/SplineCubic.js`** — `bsplineInterpolate()`: cubic B-spline interpolation for slider paths.

---

## 3. Truss Dynamics Simulation

### C++
- **`cpp/common/dynamics/Truss.h`** — `Truss` class: `TrussNode`, `TrussEdge`, `TrussFace`. Methods for adding points, edges, faces, generating girders/rings/cylinders, managing "blocks" of points/edges.
- **`cpp/common/dynamics/TrussDynamics_d.h`** — Double-precision physics engine: points, forces, velocities, bonds, neighbor lists. **Projective Dynamics**: `make_PD_Matrix`, `make_PDmat_sparse`, `rhs_ProjectiveDynamics`, `dotPD`. **Linear solvers**: CG, Cholesky (LDLT), Jacobi, Gauss-Seidel, Momentum, FIRE. **Collision**: `Buckets` for pointBBs/edgeBBs/faceBBs, `evalTrussCollisionImpulses_bonds`. **EdgeVertBonds** for slider constraints. **Strain**: `evalBondTension()` computes strain per bond — **TODO: break bond if strain > maxStrain** (line 1788). `LinSolveMethod` enum with 14 methods.
- **`cpp/common/dynamics/TrussDynamics_d.cpp`** — Implementation (2042 lines): `prepare_LinearSystem` does Cholesky decomposition (sparse + dense variants). `linsolve()` dispatches to selected solver. `evalTrussCollisionImpulses_bonds` does momentum-conserving impulses. `maxStrain` stored per bond but **breaking not implemented**.
- **`cpp/common/dynamics/TrussDynamics_f.h`** — Float-precision variant (356 lines), similar functionality.
- **`cpp/common/dynamics/Cholesky.h`** — Cholesky LDLT decomposition utilities.
- **`cpp/common/dynamics/SoftBody.h`** — SoftBody with bond force clamping (`fmax`).

### Python — Truss Dynamics & GPU Solvers (`python/pyTruss/`)

This is a major implementation of truss dynamics with OpenCL GPU solvers, mirroring the C++ `TrussDynamics_d` architecture.

#### OpenCL Kernels (`truss.cl`, 849 lines)
- **`jacobi_iteration_sparse`** — One Jacobi iteration: sparse neighbor-based, float4 positions+mass.
- **`gauss_seidel_iteration_colored`** — GS iteration with graph coloring for parallelism.
- **`jacobi_iteration_sparse_local`** — Jacobi with local memory optimization.
- **`vbd_vertex_chunk`** — Vertex Block Descent: per-workgroup vertex solve with local neighbor data. Pre-packed chunks (32 vertices, ≤128 neighbors).
- **`vbd_vertex_serial`** — VBD serial kernel: single-threaded iteration over all vertices.
- **`jacobi_iteration_diff_serial`** — Jacobi on displacement variables (dp = p - p0).
- **`pd_momentum_mix_serial`** — Momentum mixing stage (Heavy Ball / Chebyshev acceleration).
- **`jacobi_fly`** — Jacobi with on-the-fly force computation (no precomputed RHS). Iterates in-register.
- **`GS_fly`** — Gauss-Seidel with on-the-fly force, per-color sweeps with barriers.
- **`jacobi_diff`** — Jacobi on displacements with precomputed RHS, in-register iteration.
- **`GS_diff`** — Gauss-Seidel on displacements, per-color with barriers.
- **`precompute_dRHS`** / **`precompute_dRHS_d`** — Precompute RHS vector b' and diagonal Aii' for displacement-form solvers (float and double precision).

Helper functions: `init_hessian`, `accum_vertex_hessian` (3×3 Hessian assembly from spring direction), `solve3x3` (3×3 linear solve via cofactors).

#### GPU Solver (`truss_solver_ocl_new.py`, 252 lines)
- **`TrussSolverOCL`** class extends `OpenCLBase`: clean init pipeline (`_prepare_topology` → `_allocate_buffers` → `_upload_static_data` → `_setup_kernels`). Pre-generated kernel args for zero per-call overhead.
- **Solver callbacks**: `solve_vbd_serial` (VBD with serial kernel), `solve_jacobi_fly` (Jacobi on-the-fly), `solve_jacobi_diff` (Jacobi displacement with precomputed RHS).
- **`step()`**: Python-side MD loop — gravity → GPU solve → velocity update → fixed point constraints.
- **`run()`**: Multi-step runner with optional trajectory tracking.
- **`SOLVERS`** registry for easy solver switching.

#### Older GPU Solver (`truss_ocl.py`, 898 lines)
- Standalone `TrussSolverOCL` (not extending `OpenCLBase`): manual buffer management.
- `build_vbd_chunks()`: packs color groups into VBD workgroups (32 vertices + ≤128 neighbors per chunk).
- `_init_pd_system()`: projective dynamics setup with `make_pd_rhs`, `make_pd_Aii0`.

#### CPU Solver (`truss_solver.py`, 1249 lines)
- **`TrussSolver`** class: implicit Euler with pluggable solver callbacks.
- **Solvers**: `solve_vbd` (VBD), `solve_direct` (dense `np.linalg.solve`), `solve_cholesky` (LDLT), `solve_iterative_momentum` (Heavy Ball), `solve_iterative_chebyshev` (Chebyshev acceleration), `solve_jacobi_diff`, `solve_jacobi_fly`, `solve_gs_diff`, `solve_gs_fly`.
- **`estimate_iterative_spectral_radius()`**: computes ρ(T) for convergence analysis.
- Iteration logging for debugging.

#### Projective Dynamics (`projective_dynamics.py`, 198 lines)
- `make_pd_Aii0()`: diagonal mass/dt² term.
- `make_pd_matrix()`: dense system matrix assembly.
- `make_pd_rhs()`: RHS with spring terms using predicted positions.
- `solve_pd()`: full PD loop with gravity, damping, fixed points, `np.linalg.solve`.

#### Sparse Operations (`sparse.py`, 395 lines)
- `build_neighbor_list()`: bond→point-centered neighbor lists.
- `neigh_stiffness()`: convert bond stiffness to per-neighbor stiffness.
- `neighs_to_dense_arrays()`: pad to dense arrays for OpenCL.
- `make_Aii()`: diagonal of system matrix (with optional PD mass term).
- `dot_sparse()`: sparse matrix-vector product.
- `jacobi_iteration_sparse()`: CPU Jacobi iteration.
- `linsolve_Jacobi_sparse()`, `linsolve_iterative()`: iterative solvers.
- `color_graph()`: Luby's MIS-based graph coloring for parallel GS.
- `gauss_seidel_iteration_colored()`: CPU colored GS.

#### Iterative Linear Solvers (`IterativeLinearSolvers.py`, 195 lines)
- `get_iteration_matrix()`: Jacobi/GS iteration matrix and spectral radius.
- `jacobi_step()`, `gauss_seidel_step()`: single iteration steps.
- `chebyshev_accelerator()`, `momentum_accelerator()`: acceleration schemes.
- `chebyshev_omega_update()`: dynamic ω for Chebyshev.

#### Truss Data Structure (`truss.py`, 387 lines)
- **`Truss`** class: points, bonds, masses, stiffness, fixed points. `build_grid_2d()`, `get_rest_lengths()`, `get_pd_quantities()`, `get_neighbor_list()`.
- **`color_graph()`**: Vivace graph coloring (Luby's MIS algorithm) for parallel Gauss-Seidel. Based on Fratarcangeli et al. SIGGRAPH Asia 2016.

#### Test Harnesses
- **`run_vbd_cloth.py`** (291 lines) — Cloth simulation: grid truss with anchor points, CPU+GPU solvers, trajectory tracking, matplotlib visualization. Graph coloring visualization.
- **`run_solver_debug.py`** (17K) — Extensive solver debugging with multiple configurations.
- **`test_Chebyshev_accel.py`** — Chebyshev acceleration convergence test.
- **`test_Jacobi_Chebyshev_convergence.py`** — Jacobi vs Chebyshev convergence comparison.
- **`test_graph_coloring.py`** — Graph coloring validation.
- **`example_ocl_new.py`** — Usage example for new GPU solver.

#### Documentation
- **`ARCHITECTURE.md`** (327 lines) — Old vs new solver architecture comparison.
- **`IMPLEMENTATION_SUMMARY.md`** (292 lines) — Design summary of new solver.
- **`README_new_solver.md`** — Usage guide.
- **`Vivace_graphColoring_GaussSeidel.md`** — Graph coloring for parallel GS.

### Python — Constraint Solver & CG (`python/`)
- **`constrain_solver.py`** (342 lines) — `stepCG()`, `SolveCG()`: conjugate gradient solver with optional `dotFunc` for matrix-free operation. Used for constraint satisfaction in particle systems.
- **`try_CG.py`** (9K) — CG solver experiments and tests.

### Key Solver Notes
- **Cholesky (Projective Dynamics)**: Pre-computed decomposition of the system matrix. **Breaks when topology changes** (bonds broken by damage). Must re-decompose. Available in C++ and Python.
- **Jacobi**: Iterative, handles topology changes gracefully. No pre-factorization needed. Available in C++, Python CPU, and Python OpenCL.
- **CG**: Iterative with diagonal preconditioner. Also handles topology changes. Available in C++ and Python.
- **VBD (Vertex Block Descent)**: Per-vertex local solve with graph coloring for parallelism. Handles topology changes. Available in C++ and Python OpenCL.
- **Chebyshev/Momentum acceleration**: Speeds up Jacobi/GS convergence. Available in C++ and Python.
- **Graph coloring (Vivace)**: Luby's MIS algorithm partitions vertices for parallel GS. Available in Python.

---

## 4. Slider / Actuator System

### C++
- **`SpaceCraftComponents.h:750-842`** — `Slider` class: inherits `Node`, has `Path path`, motor params (`forceMax`, `powerMax`, `maxSpeed`, `speed`, `Kdv`, `vel`, `mass`), `icontrol` index for user control. `move(dt,l,v,f)` does drive dynamics with force clamping. `updatePath()` extracts path vertices from rail component side. `StructuralComponent::updateSlidersPaths()` initializes slider paths on rings/girders.
- **`SpaceCraft2Mesh2.h:261-308`** — `applySliders2sim()`: maps slider `icontrol` → `control_speed[]` array, updates `EdgeVertBond` in `TrussDynamics`. `sliders2edgeverts()`: creates `EdgeVertBond` constraints with stiffness K=100000.
- **`SpaceCraft.h:382-400`** — `make_Ring2()`: creates Ring with 4 Sliders, each with `icontrol` for user-driven control.

### JavaScript
- **`SpaceCraft.js:88-181`** — `Path` class with `interpolate()` using `bsplineInterpolate()`. `Slider` class with `pAttach` for explicit attachment, `weldDist`, `slidingVertId`.
- **`SpaceCraft2Mesh.js:226-275`** — Slider mesh generation: `pAttach` creates new vertex + welds to sliding component. Standard logic uses `vert_collapse` to find nearest rail vertex. Path built via `buildPathOnEdge` (SDF or stride method).

### Status
- **C++**: Slider motor physics implemented (force/velocity control via `icontrol`). Path extraction works for girders and rings. Attachment via `EdgeVertBond` in simulation.
- **JS**: Slider path interpolation uses cubic B-spline. Two attachment methods: `pAttach` (explicit point + weld) and standard (nearest vertex collapse). Path building via SDF cylinder selection or stride indexing.
- **Known issues**: Slider-to-girder attachment was partially solved in C++. JS version has more robust `pAttach` + `weldDist` mechanism. The `StructuralComponent::updateSlidersPaths` in C++ has debug prints and may have issues with shared paths.

---

## 5. Pipe / Line Attachment System

### C++
- **`SpaceCraftComponents.h:648-659`** — `Pipe` class: has `maxFlow`, `Path path`, pointers to components `a` and `b`. Comment mentions "Also Cables, should be attached to girders". Currently `component_kind()` returns `Rope` (reusing rope type).
- **`SpaceCraft.h:50`** — `std::vector<Pipe*> pipes; // Also Cables, should be attached to girders`

### Status
- **Pipes are defined but not fully integrated**. They have a `Path` (like sliders) but no mesh generation or simulation integration. The comment suggests they should follow girders like cables.
- **No power grid, fuel pipe, or heat pipe simulation exists yet**. The `Pipe` class is a placeholder for future line-type attachments.

---

## 6. Radiosity and Scattering

### C++
- **`cpp/common/dynamics/Radiosity.h`** — `Radiosity` class: extends `TriangleRayTracer` and `LinSolver`. `makeCouplingMatrix()` builds N×N coupling matrix with view factors + occlusion. `processTriangles()` subdivides triangles into surface elements. Iterative solver via `LinSolver`.
- **`cpp/common/dynamics/TriangleRayTracer.h`** — `TriangleRayTracer`: triangle obstacles, `addTriangle()`, `fromMesh()`, `getOcclusion()` ray test. `trinagleToElements2()` subdivides triangles into surface elements for radiosity.

### OpenCL (Shared)
- **`cpp/common_resources/cl/radiosity.cl`** — OpenCL kernels: `occlusion_matrix` (brute-force ray-triangle occlusion between element centers), `coupling_matrix` (view factor + occlusion coupling). 476 lines. Used by both C++ and Python.

### Python — Radiosity Core (`python/Radiosity/`)
- **`TriangleOcclusionRaytracer.py`** — PyOpenCL wrapper: `TriangleOcclusionOCL` class extends `OpenCLBase`. `triangulate_faces()` (fan triangulation), `build_buffers()` (points + tris packed as float4), `run()` (launches `occlusion_matrix` kernel). Test geometries: `three_quads`, `cube2faces`. 3D matplotlib visualization of visible/occluded pairs.
- **`Radiosity3D.py`** — Full 3D radiosity solver (456 lines): `build_face_frame()` (orthonormal basis per face), `project_to_uv()`, `rasterize_face_to_elems()` (hex rasterization per face → surface elements with center, normal, area). `compute_view_factors_3d()` (N×N view factor matrix with Lambertian cosine terms). `solve_radiosity_system_3d()` (solves `(I - W)B = P` with front/back emissivity split). `calculate_temperatures_3d()` (thin-sheet thermal balance → T⁴). 3D visualization with element tiles colored by temperature.
- **`Radiosity.py`** — 2D radiosity solver (420 lines): `compute_view_factors()` (midpoint inverse-square kernel with cosine projections), `solve_radiosity_system()` (front/back signed-cosine split, `(I - F)B = P`). Includes element assembly from face data with emissivity/reflectivity.
- **`PolygonRasterization.py`** — `HexRasterizerNP` class: hexagonal grid rasterization of arbitrary polygons. Clips hex cells to polygon boundary, computes area-weighted centers. Used by `Radiosity3D.py` for face subdivision. 22K.
- **`OctahedralSphereMaping.py`** — Octahedral sphere mapping for direction encoding (radiosity direction sampling).

### Python — Spacecraft Scattering & X-ray (`python/pyScatter/`)
- **`SpacecraftSimulation.cl`** — OpenCL kernels (144 lines): `compute_occlusion` (2D grid of rays between element centers, ray-triangle intersect with self-face skip), `compute_attenuation` (X-ray source → detector grid, optical depth accumulation through thin triangular sheets with material μ).
- **`SimulationPipeline.py`** — Full radiosity pipeline: `compute_occlusion()` (launches `compute_occlusion` kernel), `run_radiosity_pipeline()` (rasterize faces → elements → view factors → occlusion → solve → temperatures). Integrates `Radiosity3D` + `SpacecraftGeometry`.
- **`SpacecraftGeometry.py`** — `generate_christmas_tree()`: procedural spacecraft geometry with tilted radial branches, multi-level, configurable branch count/length/width. Returns vertices + triangular faces.
- **`xray_sim.py`** — Monte-Carlo X-ray simulation (192 lines): random tubes + triangles as primitives, spatial clustering for GPU, source-to-detector ray casting with attenuation. Cluster-based BVH for performance.
- **`AttenuationTest.py`** — Attenuation pipeline test: triangulate spacecraft geometry → `compute_attenuation` kernel → detector image. Visualization with matplotlib.
- **`MonteCarloScattere.md`** — Extensive documentation (85K) of Monte Carlo scattering physics, spacecraft geometry, simulation pipeline.
- **`run_all.py`** — Orchestrates geometry generation → radiosity → attenuation → visualization.

### Status
- **Radiosity coupling matrix computation works** in C++ (CPU), Python (NumPy), and Python (OpenCL GPU).
- **Occlusion calculation** via triangle ray-tracing is implemented and GPU-accelerated in both `radiosity.cl` (shared) and `SpacecraftSimulation.cl` (pyScatter).
- **3D radiosity solver** (`Radiosity3D.py`) is fully functional: face rasterization → view factors → linear solve → temperature calculation. Tested with cube geometry.
- **X-ray attenuation** (`pyScatter/`) computes optical depth through spacecraft geometry from point source to detector grid. GPU-accelerated.
- **Integration with spacecraft**: `EditSpaceCraft.h` includes `Radiosity radiositySolver` but it's not wired into the simulation loop. Python `pyScatter` has a `generate_christmas_tree()` spacecraft geometry but it's a test geometry, not connected to `SpaceCraft` component model.
- **Scattering**: The radiosity framework models thermal radiation scattering. The coupling matrix includes view factors with occlusion. X-ray attenuation models ionizing radiation transport through thin sheets.

---

## 7. Nozzle / Pusher Plate / Damper

### JavaScript (Primary)
- **`constructionBlockTests.js:340-436`** — "TubeSheetBond Hex Ring" test: builds hexagonal outer tube + triangular inner tube (telescopic damper) + **parabolic pusher-plate / plasma nozzle dish** above the damper. Uses `TubeSheet` for tubes, `ParametricParabolaPatch` for dish.
- **`MeshGenTestGUI.js`** — GUI for testing: `ParabolaSheet` (single-layer parabolic truss), `ParabolaSlab` (double-layer with bracing), `ParametricParabola` (annular patch with different inner/outer tessellation), `TubeSheetBond` (two tubes with bond-length analysis).
- **`MeshesUV.js`** — `ParabolaUVfunc`, `ParabolaSheet`, `ParabolaSlab_wrap`, `ParametricParabolaPatch` — UV-surface generators for parabolic dishes.

### Documentation
- **`docs/SpaceCrafting/SpaceCraftConstructionProblems.md`** — Detailed WIP document (729 lines) covering:
  - §1: Telescopic truss recoil damper design (hexagonal outer + triangular inner tube, slider rails, collision avoidance)
  - §1.4: Parabolic dish / pusher-plate generators
  - §1.5: Concrete sub-tasks [T1]-[T5] for implementation
  - §2: Plates on girders and ropes
  - §3: Welds and automatic connections
  - §4: Selection manager, SDF-based selection

### Status
- **JS has working prototypes** for the nozzle/damper system using UV mesh generators.
- **C++ has `TubeSheet`, `SlabTube` etc.** in `DrawUV.h` but the parabolic dish generators may not be ported back.
- The nozzle is built as a **parabolic mesh** using `ParametricParabolaPatch` (thick, for blade shields) and attached to a girder/damper structure.

---

## 8. MHD / Electromagnetic Simulation (Biot-Savart)

### JavaScript
- **`js/mhd_demo/physics.js`** — Full MHD plasma nozzle simulation (874 lines): `coilField(r,z,I,a)` using elliptic integrals K(m) and E(m). `makeParabolicCage()`, `makeSphericalPlasma()`. `stepSimulation()` with control-point method: plasma control points target B=0 (diamagnetic), cage control points target B=B_sc (flux conservation). `sampleFieldAtPoint()` for visualization.
- **`js/mhd_demo/render.js`** — WebGL rendering: GLSL shader with `coilField()` implemented on GPU for real-time field visualization. `ellipticK()` and `ellipticE()` in GLSL.
- **`js/mhd_demo/ui.js`** — UI with mutual inductance calculation, self-inductance, field-from-loop, flux conservation solver.
- **`js/mhd_demo/index.html`** — Interactive demo page.

### Python
- **`doc/python/MHD/inductance_core.py`** — `calc_mutual_inductance()`, `calc_self_inductance()`, `build_inductance_matrix()` using `scipy.special.ellipk/ellipe`. Core library for coil-gun and plasma nozzle.
- **`doc/python/MHD/doc/MHD_Plasma_Nozzle_Simulation.md`** — Comprehensive design document (1101 lines): cylindrical MHD, flux conservation vs control-point method, Biot-Savart with elliptic integrals, plasma dynamics, cage resistance, stability (Rayleigh-Taylor, kink), implementation roadmap.
- **`doc/python/MHD/doc/MHD_GaussGun_demo.md`** — Gauss-gun coil accelerator: driver coils, projectile disk, flux switching, 1D dynamics.
- **`doc/python/MHD/doc/MHD_Plasma_Nozle_python.md`** — Python implementation notes.
- **`doc/python/MHD/MHD_plots.py`** — Plotting utilities.

### Documentation
- **`docs/BiotSavart/MHD_Plasma_Nozzle_webgl.md`** — WebGL MHD plasma nozzle documentation.
- **`docs/BiotSavart/AeroPotentialFlow_simulation_javascript.md`** — Potential flow (related vortex-lattice method).
- **`docs/BiotSavart/PotentialFlow.md`** — Potential flow theory.

### Status
- **JS MHD demo is functional**: Real-time WebGL visualization of magnetic fields from coils, plasma expansion with flux conservation, parabolic cage geometry.
- **Python inductance core is functional**: Mutual/self inductance matrices, Biot-Savart field calculations.
- **Gauss-gun (coilgun) demo exists**: Driver coils with switching, projectile acceleration.
- **Not integrated with spacecraft simulator**: MHD is standalone. No connection to `SpaceCraft` thruster model or `TrussDynamics`.

---

## 9. Space Combat

### C++
- **`cpp/common/Orbital/spaceCombat.h`** — Weapon types: `ProjectileType`, `LaserGunType`, `SpaceGunType`, `SpaceGun`, `SpaceSalvo`, `ProjectedTarget`, `CombatAssembly`, `PulseEngine`, `SpaceShipMobility`. `whippleShieldType` for layered armor. `difractionLimit_spot()` for laser physics.
- **`cpp/common/Orbital/SpaceWorld.h`** — `SpaceWorld`: planets, ships, combatants. ODE integrator for orbital mechanics. Trajectory prediction, combat damage evaluation.
- **`cpp/common/Orbital/SpaceCombatOCL.h`** — Nearly empty wrapper (21 lines), includes only.
- **`cpp/apps/OrbitalWar/spaceTactics.cpp`** — Space combat game with N-body physics, trajectory splines, weapons.

### Encyclopedia
- **`encyclopedia/space_warfare/`** — 7 chapters: tactics, weapons, defense, propulsion, ship construction, physical model.

---

## 10. Damage Model

### Current State
- **`TrussDynamics_d.cpp:1775-1791`** — `evalBondTension()`: computes strain per bond. **TODO comment: "break the bond if strain > maxStrain"** — not yet implemented.
- **`TrussDynamics_d.h:181`** — `Vec2d* maxStrain` — stores pull/push strain limits per bond, populated in `exportSim()` from material properties.
- **`spaceCombat.h`** — `whippleShieldType` with basic impact modeling. `ProjectedTarget` with hitpoint-based damage.
- **`TriangleRayTracer.h`** — Ray-triangle intersection for projectile hit testing.

### What's Missing
- **Bond breaking**: `maxStrain` is stored but never checked/acted upon.
- **Projectile-vs-truss intersection**: No code computes probability of hitting specific truss edges.
- **Dynamic topology change**: When bonds break, Cholesky decomposition is invalidated. Jacobi/CG solvers handle this naturally.
- **Self-collision after disintegration**: `Buckets` collision system exists for point-edge collision but is not used for self-collision of broken truss fragments.

---

## 11. Self-Collision

### C++
- **`TrussDynamics_d.h:208-215`** — Collision infrastructure: `nBBs`, `BBs` (bounding boxes), `pointBBs`, `edgeBBs`, `faceBBs` (Buckets for spatial hashing). `collision_damping` parameter.
- **`TrussDynamics_d.cpp:1740-1773`** — `evalTrussCollisionImpulses_bonds()`: momentum-conserving impulses between bonded points. Only handles bond-length collisions, not self-collision between disconnected fragments.
- **`TrussDynamics_d.h:416`** — `evalTrussCollisionImpulses_bonds()` declared but no general self-collision function.

### Python — Particle Collision with OpenCL (`python/particleCollisions2D/`)
- **`constrains.cl`** (130 lines) — Two OpenCL kernels for constraint-based collision:
  - **`updateJacobi_neighs`**: One Jacobi iteration for projective dynamics with bonds + collisions. Neighbors include both bonded and colliding pairs. Collision bonds use `l0 < 0` sentinel; `Rd` parameter cuts off soft collisions to prevent instability. Weighted average of predicted positions satisfying all constraints.
  - **`updateCollisonNeighbors`**: Broad-phase collision detection using spatial groups. Local memory (32-wide) for neighbor candidate positions. Fills `neighs`/`kngs`/`l0s` arrays with collision pairs. Extended radius `Rp` as safety margin for Jacobi iteration spread.
- **`constrains_cl.py`** (151 lines) — `CLConstrains` class: OpenCL setup, buffer management, `updateJacobi_neighs()` wrapper. Platform/device selection (NVIDIA default).
- **`particle_system.py`** (18K) — CPU particle system with collision detection.
- **`particle_system_new.py`** (25K) — Refactored particle system.
- **`particle_system_o1.py`** (17K) — O(1) spatial hashing particle system.
- **`sparse.py`** (13K) — Sparse neighbor operations (shared with pyTruss).
- **`test_collision.py`** — Collision test harness.
- **`test_water_jacobi.py`** — Water simulation with Jacobi solver.
- **`particle_collisions.md`** — Documentation of particle collision approach.

### Key Collision Concepts
- **Unified bond+collision neighbors**: Both bonds and collisions use the same `neighs`/`kngs`/`l0s` arrays. Collision pairs have `l0 < 0` (negative rest length as sentinel). This allows the same Jacobi/PD solver to handle both bonds and collisions in one pass.
- **Reduced radius `Rd`**: Only hard collisions (deep penetration) are solved by constraint solver. Soft collisions handled by force field. Prevents oscillation at contact boundary.
- **Extended radius `Rp`**: Safety margin for broad-phase — includes pairs within `rc + Rp` because Jacobi iteration may move particles closer.
- **Spatial groups**: Particles organized into groups for broad-phase. Each group has a range of neighbor candidates. Local memory tiling for GPU efficiency.

### What's Missing
- **Self-collision detection for truss fragments**: The particle collision system works for point-point collisions. No code detects collisions between disconnected truss fragments (e.g., after ship disintegration) with edge-edge or edge-point geometry.
- **Point-face collision**: `faceBBs` exists in C++ but no collision response code uses it.
- **Edge-edge collision**: Not implemented in any language.

---

## 12. JS vs C++ vs Python Parity

### Where JS is Ahead
- **Mesh generators**: `ParabolaSheet`, `ParabolaSlab_wrap`, `ParametricParabolaPatch` — parabolic dish generators exist only in JS.
- **Slider attachment**: `pAttach` + `weldDist` mechanism in JS is more robust than C++ `vert_collapse`.
- **Path interpolation**: Cubic B-spline (`bsplineInterpolate`) is active in JS. C++ has `Path::getPos()` but interpolation method unclear.
- **Selection system**: `SelectionBanks` with SDF-based selection in JS. C++ has selection but less developed.
- **MeshGenTestGUI**: Interactive parameter playground for UV mesh generators — no C++ equivalent.
- **Duplicate detection**: `scanDuplicateVerts`, `scanDuplicateEdges` in JS MeshBuilder.

### Where C++ is Ahead
- **TrussDynamics**: Full physics simulation with 14 solver methods, Cholesky decomposition, FIRE optimization. JS has no dynamics. Python has dynamics but fewer solver methods.
- **SpaceCraft2Mesh export**: `exportSim()` converts mesh to `TrussDynamics_d` with material properties, mass, stiffness. JS has no simulation export. Python has no spacecraft-to-mesh conversion.
- **Lua scripting**: Full Lua API for spacecraft design. JS uses command-based API. Python has none.
- **Space combat**: `spaceCombat.h`, `SpaceWorld.h` with orbital mechanics, weapons, damage. JS has none. Python has Jupyter notebooks only.
- **Radiosity**: C++ has `Radiosity` class integrated with `EditSpaceCraft.h`. JS has no radiosity. Python has standalone radiosity but not integrated with spacecraft.
- **MHD**: Neither C++ spacecraft nor JS spacecraft integrate MHD, but JS has standalone demo. Python has inductance core but not integrated.

### Where Python is Ahead
- **OpenCL truss solvers**: `truss.cl` (849 lines, 12 kernels) — most complete OpenCL truss solver suite. VBD, Jacobi (3 variants), GS (3 variants), momentum mixing, precomputed RHS. C++ has OpenCL dynamics but less documented. JS has none.
- **Graph coloring (Vivace)**: Luby's MIS algorithm for parallel Gauss-Seidel. Not in C++ or JS.
- **Chebyshev acceleration**: Convergence analysis and spectral radius estimation. Not in C++ or JS.
- **3D radiosity solver**: `Radiosity3D.py` with face rasterization, view factors, temperature calculation. More complete than C++ for standalone use.
- **X-ray attenuation**: `pyScatter/` with Monte Carlo X-ray simulation through spacecraft geometry. Not in C++ or JS.
- **Particle collision OpenCL**: `constrains.cl` with unified bond+collision neighbor system. C++ has `Buckets` but no GPU collision kernel. JS has none.
- **Conjugate gradient**: `constrain_solver.py` with matrix-free CG (`dotFunc`). Not in C++ or JS as standalone.

### Porting Considerations
- JS MeshBuilder is a near-complete port of C++ `Mesh::Builder2` (2017 lines).
- JS `BuildCraft_blocks_js` mirrors C++ `BuildCraft_blocks` but with enhancements (plates, pAttach sliders).
- **Porting JS mesh generators back to C++** would benefit from C++ visual debugging (OpenGL renderer, `SpaceCraftDraw.h`).
- A **simple mesh viewer in C++** that can examine parameters and groups would help bridge the gap.
- **Python `pyTruss` OpenCL kernels** could be ported to C++ `TrussDynamics_d` OpenCL path (`spaceCraftDynamicsOCL.cpp`).
- **Python `constrains.cl`** unified bond+collision approach could improve C++ collision handling.
- **Python `Radiosity3D.py`** could be integrated with C++ `Radiosity.h` for validation.

---

## 13. File Cross-Reference Index

| Subsystem | C++ | JS | Python | OpenCL |
|-----------|-----|-----|--------|--------|
| Component model | `SpaceCraft.h`, `SpaceCraftComponents.h` | `SpaceCraft.js` | — | — |
| Lua scripting | `EditSpaceCraft.h` | — | — | — |
| Mesh builder | `MeshBuilder2.h` | `MeshBuilder.js` | — | — |
| Mesh generation | `SpaceCraft2Mesh2.h`, `SpaceCraft2Mesh_blocks.h` | `SpaceCraft2Mesh.js`, `MeshGenerators.js`, `MeshesUV.js` | — | — |
| UV generators | `DrawUV.h` (C++ apps) | `MeshGenerators.js`, `MeshesUV.js` | — | — |
| Truss structure | `Truss.h` | — | `pyTruss/truss.py` | — |
| Truss dynamics | `TrussDynamics_d.h/.cpp`, `TrussDynamics_f.h/.cpp` | — | `pyTruss/truss_solver.py`, `pyTruss/projective_dynamics.py` | `pyTruss/truss.cl` |
| GPU truss solver | `spaceCraftDynamicsOCL.cpp` (experimental) | — | `pyTruss/truss_solver_ocl_new.py`, `pyTruss/truss_ocl.py` | `pyTruss/truss.cl` (12 kernels) |
| Cholesky solver | `Cholesky.h`, `Lingebra.h` | — | `pyTruss/truss_solver.py` (`solve_cholesky`) | — |
| CG solver | — | — | `constrain_solver.py`, `try_CG.py` | — |
| Sparse operations | — | — | `pyTruss/sparse.py`, `particleCollisions2D/sparse.py` | — |
| Iterative solvers | `TrussDynamics_d.cpp` (Jacobi, GS, Momentum, FIRE) | — | `pyTruss/IterativeLinearSolvers.py`, `pyTruss/truss_solver.py` | `pyTruss/truss.cl` (Jacobi, GS, VBD, Momentum) |
| Graph coloring | — | — | `pyTruss/truss.py` (`color_graph`), `pyTruss/sparse.py` | `pyTruss/truss.cl` (`gauss_seidel_iteration_colored`) |
| Chebyshev accel | `TrussDynamics_d.cpp` | — | `pyTruss/IterativeLinearSolvers.py`, `pyTruss/test_Chebyshev_accel.py` | — |
| Slider/actuator | `SpaceCraftComponents.h`, `SpaceCraft2Mesh2.h` | `SpaceCraft.js`, `SpaceCraft2Mesh.js` | — | — |
| Spline interpolation | — | `SplineCubic.js` | — | — |
| Pipes/lines | `SpaceCraftComponents.h` (stub) | — | — | — |
| Radiosity | `Radiosity.h`, `TriangleRayTracer.h` | — | `Radiosity/Radiosity3D.py`, `Radiosity/Radiosity.py` | `radiosity.cl` (shared) |
| Radiosity occlusion | `TriangleRayTracer.h` | — | `Radiosity/TriangleOcclusionRaytracer.py` | `radiosity.cl` (`occlusion_matrix`) |
| Polygon rasterization | — | — | `Radiosity/PolygonRasterization.py` (`HexRasterizerNP`) | — |
| X-ray attenuation | — | — | `pyScatter/xray_sim.py`, `pyScatter/AttenuationTest.py` | `pyScatter/cl/SpacecraftSimulation.cl` (`compute_attenuation`) |
| Spacecraft geometry (test) | — | — | `pyScatter/SpacecraftGeometry.py` | — |
| Nozzle/damper | — | `constructionBlockTests.js`, `MeshGenTestGUI.js` | — | — |
| MHD/Biot-Savart | — | `mhd_demo/physics.js`, `render.js`, `ui.js` | `doc/python/MHD/inductance_core.py` | — |
| Space combat | `spaceCombat.h`, `SpaceWorld.h` | — | `projects/SpaceCombat/` | `SpaceCombatOCL.h` (stub) |
| Damage | `TrussDynamics_d.cpp` (TODO), `spaceCombat.h` | — | — | — |
| Collision (truss) | `TrussDynamics_d.h` (partial, `Buckets`) | — | — | — |
| Collision (particle) | — | — | `particleCollisions2D/particle_system*.py` | `particleCollisions2D/constrains.cl` (`updateJacobi_neighs`, `updateCollisonNeighbors`) |
| Rendering | `SpaceCraftDraw.h` | `mhd_demo/render.js` | `Radiosity/Radiosity3D.py` (matplotlib) | — |
| GUI | `SpaceCraftGUI.h` | `MeshGenTestGUI.js` | — | — |
| Construction docs | — | — | — | — |
| Construction problems | `SpaceCraftConstructionProblems.md` | — | — | — |
