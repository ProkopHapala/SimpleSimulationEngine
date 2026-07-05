# Python Scripts Audit — Duplication & Cleanup

**Created**: 2025-07-04
**Status**: Audit complete, awaiting decisions

---

## Summary

The `python/` and `doc/python/` directories contain **~200 Python files** across ~30 directories. Many are iterative versions of the same experiment (v1, v2, v3...) with no clear indication of which is current or working. Cross-directory duplication exists in shared utilities (sparse.py, GaussAtom files).

**Key problems**:
- Numbered version chains with no "current" marker (e.g. `PipeValvelesEngine.py` through `PipeValvelesEngine8.py`)
- Cross-directory copy-paste of shared utilities (`sparse.py`, `BoysFunc.py`, etc.)
- Dead files (`python/p` — empty file, `python/test.json` — test artifact)
- Old experiment superseded by newer refactor but not removed (GLCL → GLCL2, terrain_ocl → terrain_ocl_2)
- Bug in `particle_system_new.py:11` — `ir*irE_Coulomb` (missing `*`)

---

## Directory-by-Directory Audit

### 1. `python/GLCL/` vs `python/GLCL2/` — GLCL Integration

**GLCL/** (3 files): `GLGUI.py` (845 lines), `NBody_glcl.py`, `NBody_glsl.py`
- Monolithic PyQt5+OpenGL+OpenCL N-body demo
- `GLGUI.py` is a huge monolithic OpenGL helper, imports everything inline

**GLCL2/** (6 files + subdirs): `GLCLGUI.py`, `GLCLBrowser.py`, `ModularGL.py`, `OCLsystem.py`, `OGLsystem.py` + `cl/`, `doc/`, `scripts/`, `shaders/`
- Properly refactored modular version: separate `OGLsystem.py` (GL), `OCLsystem.py` (CL), `GLCLGUI.py` (widget)
- `ModularGL.py` has stub classes that `raise RuntimeError` — incomplete but structured

**Verdict**: GLCL2 is the successor. **GLCL/ is obsolete** — `NBody_glcl.py` and `NBody_glsl.py` are standalone demos superseded by GLCL2's modular architecture. `GLGUI.py` functionality is split into `OGLsystem.py` + `GLCLGUI.py`.

**Action**: Delete `python/GLCL/`. Keep `python/GLCL2/`.

---

### 2. `python/terrain_ocl/` vs `python/terrain_ocl_2/` — Terrain Erosion

**terrain_ocl/** (12 items): `terrain.py` (680 lines), `terrain.cl` (40K), 4 test scripts, docs, `OUT` artifact
- `test_errosion.py` — **entirely commented out** (100 lines, all `#`)
- `test_errosion_2.py` — imports `GPUShadertoyErosion` from `terrain` (may not exist in current terrain.py)
- `test_terrain.py`, `test_watershed.py`, `verify_real.py` — additional test scripts

**terrain_ocl_2/** (5 items): `terrain.py` (416 lines), `terrain.cl` (3K), `test_terrain.py`, `verify_logic.cl`, `tmp.md`
- Much smaller, refactored. `terrain.py` has shared utilities (hash, smoothstep) without GPU solver class
- `terrain.cl` is 10x smaller (3K vs 40K) — likely a stripped-down or different approach

**Verdict**: terrain_ocl_2 is a refactor/iteration. terrain_ocl has more complete docs (`terrain_ocl.md` 37K, `Watershed_Transform.md`, `Accumulative_Signed_Distance_Field_Errosion.md`). The docs should be preserved regardless.

**Action**: Keep terrain_ocl_2 as current. Move docs from terrain_ocl to terrain_ocl_2 or to `doc/`. Delete `test_errosion.py` (all commented out). Verify `test_errosion_2.py` still works with terrain_ocl_2's terrain.py.

---

### 3. `python/OrbitalWar/` — Trajectory Optimization

**5 files**: `trajectoryOptmization.py` (6.9K), `.cl` (9.3K), `.md` (77K), `test_trajectoryOptmization.py` (6K), `spacetactics.py` (16K)

- `trajectoryOptmization.py` + `.cl` — OpenCL-based orbital trajectory optimizer
- `spacetactics.py` — space combat tactics model (separate concern, 16K)
- `test_trajectoryOptmization.py` — test script
- `.md` is very detailed (77K)

**Verdict**: Self-contained, no duplication issues. `spacetactics.py` may belong in `pyCombatModels/` instead.

**Action**: Keep. Consider moving `spacetactics.py` to `pyCombatModels/`.

---

### 4. `python/pyTruss/` — Truss/FEM Solver (WORST duplication)

**23 files** — massive version proliferation:

**Core data model**:
- `truss.py` (387 lines) — `Truss` class with `build_rope`, `build_cloth` etc. **Used by all other files**

**CPU solvers**:
- `truss_solver.py` (1249 lines) — main CPU solver with VBD, Projective Dynamics, Jacobi, Gauss-Seidel. **Current/best**
- `projective_dynamics.py` (198 lines) — standalone PD functions (imported by truss_ocl.py)
- `projective_dynamics_iterative.py` (170 lines) — iterative PD with Chebyshev acceleration

**GPU/OCL solvers** — 3 versions!:
- `truss_ocl.py` (898 lines) — **old** OCL solver, imports from `projective_dynamics` and `sparse`
- `truss_solver_ocl.py` (355 lines) — **current** OCL solver, used by `run_vbd_cloth.py`
- `truss_solver_ocl_new.py` (252 lines) — **newest** OCL solver, uses `OpenCLBase` from pyMolecular, used by `run_vbd_cloth_new.py`

**Runner scripts** — 3 versions:
- `run_vbd_cloth_old.py` (213 lines) — uses `truss_ocl` (old solver). **Obsolete**
- `run_vbd_cloth.py` (291 lines) — uses `truss_solver_ocl` (current). **Active**
- `run_vbd_cloth_new.py` (301 lines) — uses `truss_solver_ocl_new` (newest, OpenCLBase). **Active**

**Shared utilities**:
- `sparse.py` (394 lines) — neighbor lists, graph coloring, Jacobi/GS sparse solvers. **Duplicated** in `particleCollisions2D/sparse.py` (348 lines, diverged)
- `IterativeLinearSolvers.py` (7.4K) — standalone iterative solvers
- `plot_utils.py` (2.3K) — plotting helper

**Tests**:
- `test_Chebyshev_accel.py`, `test_Jacobi_Chebyshev_convergence.py`, `test_graph_coloring.py`
- `run_solver_debug.py` (17K) — debug runner

**Docs**:
- `ARCHITECTURE.md`, `IMPLEMENTATION_SUMMARY.md`, `README_new_solver.md`, `Vivace_graphColoring_GaussSeidel.md`

**Verdict**: The evolution is: `truss_ocl.py` → `truss_solver_ocl.py` → `truss_solver_ocl_new.py`. Similarly: `run_vbd_cloth_old.py` → `run_vbd_cloth.py` → `run_vbd_cloth_new.py`. The old versions should be removed. `sparse.py` is duplicated with particleCollisions2D.

**Action**:
- Delete: `truss_ocl.py`, `run_vbd_cloth_old.py` (clearly superseded)
- Keep: `truss_solver.py`, `truss_solver_ocl.py`, `truss_solver_ocl_new.py`, `run_vbd_cloth.py`, `run_vbd_cloth_new.py`
- Consolidate: `sparse.py` should be shared — extract to `python/common/` or similar
- Clarify: Are both `truss_solver_ocl.py` and `truss_solver_ocl_new.py` still needed, or is `_new` the final version?

---

### 5. `python/pyMolecular/` — Molecular Simulation

**39 items** — well-structured but with some duplication:

**Core modules**: `GaussAtom.py`, `BoysFunc.py`, `GaussKinetic.py`, `GaussProduct.py`, `GaussProductDerivs.py`, `quadratureCoefs.py`, `splines.py`, `integral_NumCyl.py`, `atomicUtils.py`, `elements.py`
- **All 8 of these are identical copies** of files in `python/pyGaussAtom/` (verified by diff)

**CLCFGO variants**:
- `CLCFGO.py` (524 lines) — main module with ctypes C bindings
- `CLCFGO_coulomb_derivs.py` (467 lines) — coulomb derivative calculations
- `CLCFGO_normalization_derivs.py` (247 lines) — normalization derivatives
- `CLCFGO_normalization_derivs_2.py` (216 lines) — **v2 of above**, adds matplotlib plotting. Nearly identical header.

**eFF variants**:
- `eFF.py` (188 lines) — ctypes wrapper for C eFF library
- `eFF_terms.py` (374 lines) — energy term calculations
- `eFF_terms_derivs.py` (388 lines) — **derivative version** of eFF_terms.py (same header/docstring)
- `eFF_KineticAndOverlap.py` (176 lines) — kinetic energy & overlap calculations

**OCL subsystem**: `OCL/` — `MMFF.py`, `MMparams.py`, `MolecularDynamics.py`, `OpenCLBase.py`, `clUtils.py`, `run_scanNonBond.py`, `cl/`

**Other**: `MMFF-ref.py` (70K — reference implementation), `ReactiveFF.py`, `RigidMol.py`, `c_interface.py`, `plotUtils.py`, `testing.py`, `examples/`

**Verdict**: pyGaussAtom is a **pure subset** of pyMolecular — all 8 files are identical. The CLCFGO_normalization_derivs_2.py is a minor iteration on _derivs.py. eFF_terms.py and eFF_terms_derivs.py are related but distinct (values vs derivatives).

**Action**:
- Delete: `python/pyGaussAtom/` entirely (subset of pyMolecular)
- Merge: `CLCFGO_normalization_derivs_2.py` into `CLCFGO_normalization_derivs.py` (keep the newer one)
- Keep: everything else

---

### 6. `python/pyGaussAtom/` — DUPLICATE of pyMolecular subset

**8 files**, all identical copies of files in `pyMolecular/`:
- `BoysFunc.py`, `GaussAtom.py`, `GaussKinetic.py`, `GaussProduct.py`, `GaussProductDerivs.py`, `integral_NumCyl.py`, `quadratureCoefs.py`, `splines.py`

Only difference: `GaussAtom.py` in pyMolecular has 4 extra lines (`r = None`, `S = None`, `k_h2m = 0.1`).

**Verdict**: **Pure duplicate. Delete entirely.**

---

### 7. `python/pyShaderToy/` vs `python/pySymGLSL/` — Shader Experiments

**pyShaderToy/** (55 items): `GUI.py`, `ModelerGUI.py`, `ModelerGUI_new.py`, `FormulaNode.py`, `ExpressionParser.py`, `CSGLang.py`, `GLUtils.py` + 43 shader files
- `GUI.py` and `ModelerGUI.py` — nearly identical headers (same TODO, same sources)
- `ModelerGUI_new.py` — uses `moderngl` instead of raw OpenGL, imports `FormulaNode`
- Large shader collection in `shaders/` (43 items)

**pySymGLSL/** (18 items): `GLSL_GUI.py`, `GLSL_Simulation.py`, `__init__.py`, `pipelines/`, `shaders/`, `doc/`
- More structured, has `__init__.py` (importable package)
- `GLSL_Simulation.py` — simulation pipeline
- `GLSL_GUI.py` — GUI widget

**Verdict**: Both are active but serve different purposes. pyShaderToy is a ShaderToy-like explorer, pySymGLSL is a simulation framework. Within pyShaderToy, `GUI.py` and `ModelerGUI.py` overlap significantly, and `ModelerGUI_new.py` supersedes `ModelerGUI.py`.

**Action**:
- Delete: `pyShaderToy/ModelerGUI.py` (superseded by `ModelerGUI_new.py`)
- Keep: both directories otherwise
- Investigate: whether `GUI.py` is still used or superseded

---

### 8. `python/pySimE/` — Simulation Engine (109 items)

**Structure**:
- `__init__.py`, `constants.py`, `common_physics.py`, `utils.py`
- `chemistry/` — `common.py`, `data/`, `tests/`
- `math/` — 1 item
- `physics/` — 1 item
- `space/` — **92 items** (the bulk)

**`space/` breakdown**:
- `KosmoSuit_main.py`, `Rocket.py`, `railGun.py`, `__init__.py`, `debug.py`
- `data/` — `materials.py`, `table_ChemicalFuels.py`, `table_Materials.py`, `table_SolarSystem.py`
- `tests/` (7 files) — `OrionProjectBalistic.py`, `SkyHookNumeric.py`, `SkyHookNumeric_2.py`, `sapce_elevator.py`, `sapce_elevator_sympy.py`, `tetherSling.py`, `integrationTest.py`
- `exp/` (75 items) — **two experiment subdirectories**:

**`exp/OrbitalTransferOpt/`** (48 items):
- Parent dir has ~28 files: various optimization approaches (Simplex, polynomial, spline, Cos basis, cylindrical MC, etc.)
- `Simplex2.py` vs `Simplex2_prokop.py` — nearly identical (1 function name diff: `MCBias2b_Run` vs `MCBias2_Run`)
- `Poly4th_numeric.py` duplicated in `OrbOpt_Cos/` (identical)
- `Random_optimization.py` duplicated in `OrbOpt_Cos/` (diverged)
- `Simplex_optimization.py` duplicated in `OrbOpt_CubicSplineSimplex/` and `OrbOpt_map/` (identical copies)
- `CubicSpline.py` vs `OrbOpt_map/ClampedCubicSpline.py` — diverged
- Subdirs `OrbOpt_Cos/`, `OrbOpt_CubicSplineSimplex/`, `OrbOpt_map/` are experiment snapshots with duplicated dependencies

**`exp/pykep/`** (27 files):
- Uses external `PyKEP` library (orbital mechanics)
- `lambert_Fit.py` through `lambert_Fit_5.py` — 5 iterations of Lambert problem fitting
- `analyze_Asteroids.py` vs `analyze_Asteroids2.py` — v1 reads from `/home/prokop/Desktop/`, v2 reads from `D:\know\...` — **both have hardcoded paths, likely broken**
- `analyze_AsteroidsSize.py` vs `analyze_Asteroids_size.py` — near-duplicate names
- `sort_asteroids_aLxLy.py` vs `sort_asteroids_aLxLy_2.py` — versioned

**Verdict**: The `exp/` directory is a research graveyard. Many files are experiment snapshots with hardcoded paths. The OrbitalTransferOpt subdirectories are frozen experiment states with duplicated dependencies.

**Action**:
- Delete: `Simplex2_prokop.py` (trivial rename of `Simplex2.py`)
- Delete: duplicated `Simplex_optimization.py` copies in subdirs (identical)
- Delete: duplicated `Poly4th_numeric.py` in `OrbOpt_Cos/` (identical)
- Fix or delete: `analyze_Asteroids.py` and `analyze_Asteroids2.py` (hardcoded paths)
- Consolidate: `lambert_Fit_5.py` is likely the final version — archive or delete earlier versions
- Consider: moving `exp/` to an `archive/` subdirectory

---

### 9. `doc/python/Burn1D/` — 1D Combustion/Flow (34 items, WORST version chain)

**5 independent experiment chains**, each with many versions:

**EulerianTube chain** (4 files):
- `EulerianTube.py` (418 lines), `EulerianTube2.py` (391 lines), `EulerianTube3.py` (452 lines), `EulerianTube3_ustream.py` (466 lines)
- v1→v2: minor refactor (same class signature). v3: changed constructor signature (added `u_amb`, `bc_mode`). `_ustream`: variant of v3.
- **Latest**: `EulerianTube3.py` / `EulerianTube3_ustream.py`

**PipeValvelesEngine chain** (11 files!):
- `PipeValvelesEngine.py` (166 lines) → `PipeValvelesEngine2.py` → ... → `PipeValvelesEngine8.py` (450 lines)
- Plus `PipeValvelesEngine_intertia.py`, `_intertia2.py`, `_intertia3.py` (3 more)
- v1: simple script, no CLI. v8: full CLI with argparse, LineCollection, diffusion. Clear evolution.
- **Latest**: `PipeValvelesEngine8.py` (most features, CLI, diffusion)

**SphericalImplosion chain** (4 files):
- `SphericalImplosion.py` (377 lines) → v2 → v3 → `SphericalImplosion4.py` (412 lines)
- v4 adds argparse, numba `@njit`, PolyCollection visualization
- **Latest**: `SphericalImplosion4.py`

**TimeDependentRadiationTransfer chain** (6 files):
- `TimeDependnetRadiationTransferSpherical.py` (6.9K) → v2 → v3 → v4 → v5 → v6 (18K)
- **Latest**: `TimeDependnetRadiationTransferSpherical6.py`

**wave_in_tube chain** (4 files):
- `wave_in_tube.py` (173 lines) → v2 → v3 → `wave_in_tube4.py` (254 lines)
- v4 adds argparse, LineCollection, asymmetric geometry
- **Latest**: `wave_in_tube4.py`

**Verdict**: 34 files, only ~5 are "latest". The rest are intermediate versions. The `.md` docs are valuable and should be kept.

**Action**: For each chain, keep only the latest version + docs:
- Keep: `EulerianTube3.py`, `EulerianTube3_ustream.py` (if still needed separately)
- Keep: `PipeValvelesEngine8.py` only
- Keep: `SphericalImplosion4.py` only
- Keep: `TimeDependnetRadiationTransferSpherical6.py` only
- Keep: `wave_in_tube4.py` only
- Keep: all `.md` files
- Delete: all intermediate versions (v1-v7 etc.)
- Delete: `OUT` artifact

---

### 10. `doc/python/PotentialFlow/` — Aerodynamics (14 items)

**3 experiment chains**:

**SlenderBodyTheory chain** (7 files):
- `SlenderBodyTheory.py` (116 lines) → `SlenderBodyTheory2.py` → v3 → v4 → `SlenderBodyTheory5-a.py` (466 lines), `SlenderBodyTheory5-b.py` (17090 bytes)
- v5-a imports from `SlenderBodyTheory2` (builds on v2). v5-b is a companion.
- **Latest**: `SlenderBodyTheory5-a.py` + `SlenderBodyTheory5-b.py`

**airfoil chain** (4 files):
- `airfoil.py` (546 lines) → `airfoil2.py` → `airfoil3.py` → `airfoil4.py` (525 lines)
- v4 adds argparse, Legendre polynomials, UnivariateSpline
- **Latest**: `airfoil4.py`

**simple_vlm chain** (2 files):
- `simple_vlm.py` (367 lines) → `simple_vlm_2.py` (416 lines)
- v2 adds argparse, finite-core vortex model (Vatistas/Scully)
- **Latest**: `simple_vlm_2.py`

**Also**: `PotentialFlow.py` (13K) — standalone, `PotentialFlow.md` (150K — massive doc)

**Verdict**: 14 files, ~5 are latest. Similar pattern to Burn1D.

**Action**: Keep only latest + docs:
- Keep: `SlenderBodyTheory2.py` (needed by v5-a), `SlenderBodyTheory5-a.py`, `SlenderBodyTheory5-b.py`
- Keep: `airfoil4.py`
- Keep: `simple_vlm_2.py`
- Keep: `PotentialFlow.py`, `PotentialFlow.md`
- Delete: `SlenderBodyTheory.py`, `SlenderBodyTheory3.py`, `SlenderBodyTheory4.py`, `airfoil.py`, `airfoil2.py`, `airfoil3.py`, `simple_vlm.py`

---

### 11. `doc/python/MHD/` — Magnetohydrodynamics (20 items)

**8 Python files** + `doc/` (6 items) + `tmp.md`:

- `demo_dipole_gemini.py` (115 lines) — simple demo using `ellipk`, `ellipe`
- `demo_dipole_gemini2.py` (303 lines) — **completely different** (superconducting nozzle + diamagnetic bubble). Not a version chain — different demos despite the name.
- `demo_coil_motion_flux.py`, `demo_diamagnetic_dipole_bubble.py`, `demo_gauss_gun_flux.py`, `demo_plasma_dynamics_flux.py`, `demo_plasma_sphere_flux.py` — distinct demos
- `check_B_consistency.py`, `check_B_kernels.py` — validation scripts
- `MHD_plots.py` — plotting utilities
- `fast_eliptke_Abramowitz.py` — elliptic integral implementation
- `inductance_core.py` — inductance calculations
- `vector_potentil_coil.py` — vector potential (note: typo in name "potentil")

**Verdict**: Less duplication than other dirs. Each demo is a distinct physics scenario. The `gemini` suffix suggests AI-generated. `tmp.md` (32K) likely should be formalized or removed.

**Action**:
- Keep: all demos (they're distinct, not version chains)
- Fix: rename `vector_potentil_coil.py` → `vector_potential_coil.py`
- Review: `tmp.md` — formalize into proper doc or delete

---

### 12. `python/particleCollisions2D/` — 2D Particle System (12 items)

**3 particle system versions**:
- `particle_system.py` (441 lines) — `Particle` class, `Bond` class, basic system
- `particle_system_o1.py` (431 lines) — O(1) variant, nearly identical class structure
- `particle_system_new.py` (671 lines) — **different approach** (functional, not class-based). Has **BUG on line 11**: `ir*irE_Coulomb` should be `ir*ir*E_Coulomb`

**Other**: `constrains.cl`, `constrains_cl.py`, `sparse.py` (348 lines — **diverged copy** of pyTruss/sparse.py), `plot_utils.py`, `visualizer.py`, `test_collision.py`, `test_water_jacobi.py`, `sparse.py`

**Verdict**: 3 versions of particle system. `_new` has a bug. `sparse.py` is a diverged duplicate.

**Action**:
- Fix: `particle_system_new.py:11` — `ir*irE_Coulomb` → `ir*ir*E_Coulomb`
- Consolidate: `sparse.py` with pyTruss version
- Clarify: which particle_system is current?

---

### 13. `python/Radiosity/` — Radiosity (11 items)

- `Radiosity.py` (17K), `Radiosity3D.py` (17K), `PolygonRasterization.py` (22K), `OctahedralSphereMaping.py` (7.5K), `TriangleOcclusionRaytracer.py` (11K)
- Docs: `Radiosity.md`, `Radiosity3D.md`, `PolygonRasterization.md`, `OctahedralSphereMaping.md`
- Artifacts: `OUT-polyraster`, `OUT-radiosity3D`, `OUT-radiosity3D.png`

**Verdict**: Well-organized, each file is a distinct module. No duplication.

**Action**: Keep. Delete `OUT-*` artifacts (gitignored anyway).

---

### 14. `python/pyCombatModels/` — Combat Models (5 items)

- `CombatModel.py` (4.3K), `AirCombat.py` (1.7K), `LandCombat.py` (4.7K), `SpaceCombat.py` (3.9K), `__init__.py`
- Clean, no duplication.

**Action**: Keep as-is.

---

### 15. `python/pyRay/` — Ray Tracer (11 items)

- `GUI.py`, `common.py`, `image.py`, `ocl.py`, `scene.py`, `__init__.py`, `cl/`, `tests/`
- Clean modular structure.

**Action**: Keep as-is.

---

### 16. `python/pyScatter/` — X-ray Scattering (9 items)

- `xray_sim.py`, `SimulationPipeline.py`, `SpacecraftGeometry.py`, `AttenuationTest.py`, `run_all.py`
- `cl/` — OpenCL kernels
- Images: `attenuation_preview.png`, `geom_preview.png`, `radiosity_heatmap.png` (gitignored)

**Verdict**: Clean, no duplication.

**Action**: Keep as-is.

---

### 17. Other small directories

| Directory | Files | Status | Action |
|-----------|-------|--------|--------|
| `python/pyApprox/` | 4 | Clean — integration, radial, rational, taylor | Keep |
| `python/pyFlight/` | 2 | Just `c_interface.py` + `__init__.py` | Keep |
| `python/pyMeta/` | 3 | `cpp_utils.py`, `metaprograming.py` | Keep |
| `python/pyMetaC/` | 4 | `base.py`, `metaC.py`, `tests.py` | Keep |
| `python/pyPapers/` | 2 | `bibtex.py`, `mendelay.py` | Keep |
| `python/pyVis3D/` | 2 | Just `c_interface.py` + `__init__.py` | Keep |
| `python/pyplot3D/` | 2 | Just `main.py` + `__init__.py` | Keep |
| `python/pyShock/` | 1 | Just `c_interface.py` (no `__init__.py`!) | Keep, add `__init__.py` |
| `python/subPixelContour/` | 2 | `subPixelContour.py` (42K) + `.md` | Keep |
| `python/benchPython/` | 2 | Benchmark scripts | Keep |
| `python/patterns/` | 1 | `patterns.py` (12.6K) | Keep |
| `python/LandCraft/` | 1 | `landcraft.py` (14K) | Keep |
| `python/EulerianImpacFluid/` | 8 | Fluid sim + `.cl` + docs | Keep |
| `python/GUI/` | 1 | `GUITemplate.py` | Keep |
| `python/utils/` | 3 | C++ header extraction utilities | Keep |

---

### 18. Root-level `python/` files

**test_*.py scripts** (17 files):
- `test_Fly.py`, `test_FlyView.py`, `test_FlyViewKey.py` — flight simulation tests
- `test_KosmoSuite_*.py` (5 files) — KosmoSuite tests (fission pulse, space launch, elmag, nbody, ship accel)
- `test_Molecular.py`, `test_MolecularDist.py`, `test_RigidMol.py`, `test_RigidMol_Campher.py` — molecular tests
- `test_ReactiveFF.py`, `test_SurfConf.py`, `test_Vis3D.py` — other tests
- `test_chem_entalpy.py`, `test_chem_fuels2table.py`, `test_chem_reaction.py` — chemistry tests
- `test_Multipole.py`, `test_Metaprograming.py`, `test_common_physics.py`, `test_physics_relativisticKineticEnergy.py`

**Other root files**:
- `BaseGUI.py` (5.8K) — base GUI class
- `PGS.py` (2K) — Projected Gauss-Seidel
- `constrain_solver.py` (10.3K) — constraint solver
- `try_CG.py` (9.4K) — conjugate gradient experiment
- `try_Linearized_Move.py` (5.6K) — linearized movement experiment

**Dead files**:
- `python/p` — **empty file, 0 bytes. Delete.**
- `python/test.json` — test artifact for pySymGLSL pipeline. Delete or move to test data.
- `python/OUT-scanNonBond` — output artifact (gitignored). Delete.

**Verdict**: The 17 test scripts are scattered at root level — they test modules from various subdirectories (pyMolecular, pySimE, etc.) but live in `python/` root. This is messy but functional.

**Action**:
- Delete: `python/p` (empty)
- Delete: `python/OUT-scanNonBond` (artifact)
- Move: `test.json` to `python/pySymGLSL/test_data/` or delete
- Consider: moving test scripts closer to their modules

---

## Cross-Directory Duplication Summary

| Duplicated file | Locations | Status |
|----------------|-----------|--------|
| `sparse.py` | `pyTruss/` (394 lines), `particleCollisions2D/` (348 lines) | Diverged copies — should be shared |
| `BoysFunc.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `GaussAtom.py` | `pyGaussAtom/`, `pyMolecular/` | Near-identical (pyMolecular has 4 extra lines) |
| `GaussKinetic.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `GaussProduct.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `GaussProductDerivs.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `integral_NumCyl.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `quadratureCoefs.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `splines.py` | `pyGaussAtom/`, `pyMolecular/` | Identical |
| `Simplex_optimization.py` | `OrbitalTransferOpt/`, `OrbOpt_CubicSplineSimplex/`, `OrbOpt_map/` | Identical copies |
| `Poly4th_numeric.py` | `OrbitalTransferOpt/`, `OrbOpt_Cos/` | Identical |
| `Random_optimization.py` | `OrbitalTransferOpt/`, `OrbOpt_Cos/` | Diverged |

---

## Priority Action List

### Immediate deletes (dead/superseded files)
1. `python/p` — empty file
2. `python/OUT-scanNonBond` — artifact
3. `python/pyGaussAtom/` — entire directory (subset of pyMolecular)
4. `python/GLCL/` — entire directory (superseded by GLCL2)
5. `python/terrain_ocl/test_errosion.py` — all commented out
6. `python/pyTruss/truss_ocl.py` — superseded by truss_solver_ocl.py
7. `python/pyTruss/run_vbd_cloth_old.py` — superseded
8. `python/pyShaderToy/ModelerGUI.py` — superseded by ModelerGUI_new.py

### Version chain cleanup (keep only latest)
9. `doc/python/Burn1D/` — delete intermediate versions (keep latest of each chain + docs)
10. `doc/python/PotentialFlow/` — delete intermediate versions (keep latest + docs)
11. `python/pySimE/space/exp/OrbitalTransferOpt/` — delete duplicated files in subdirs

### Bug fixes
12. `python/particleCollisions2D/particle_system_new.py:11` — `ir*irE_Coulomb` → `ir*ir*E_Coulomb`

### Consolidation
13. Extract `sparse.py` to shared location (used by pyTruss + particleCollisions2D)
14. Move `python/OrbitalWar/spacetactics.py` to `python/pyCombatModels/`
15. Fix hardcoded paths in `pySimE/space/exp/pykep/analyze_Asteroids*.py`
16. Rename `doc/python/MHD/vector_potentil_coil.py` → `vector_potential_coil.py`

### Needs decision
17. `python/pyTruss/` — are both `truss_solver_ocl.py` and `truss_solver_ocl_new.py` needed?
18. `python/particleCollisions2D/` — which particle_system is current?
19. `python/pySimE/space/exp/` — archive the entire exp/ directory?
20. `python/test.json` — keep as test data or delete?
