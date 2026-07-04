# Topical Audit — Master Index

Cross-language maps connecting implementations of the same scientific topic across the repo.
One file per topic, following the format defined in `doc/AGENTs/skills/doc-audit/SKILL.md`.

## How to Use This Index

- Each topic below links to its audit file (if created) or lists **entry points** to investigate.
- **Status**: `created` = audit file exists, `pending` = not yet audited, `n/a` = too minor for full audit.
- **Entry points** are the most relevant files/docs to read first when diving into a topic.
- Filter out `node_modules/`, `*_bak/`, and ShaderToy `.glsl.md` snippets — these are not primary sources.

## Table of Contents

- **[Methods & Tools](#methods--tools)** — algorithms, math, data structures (the "means")
  - [Numerical Math](#numerical-math)
    - [Linear Algebra & Iterative Solvers](linear-algebra-solvers.md)
    - [Optimization](optimization.md)
    - [ODE Integration](ode-integration.md)
    - [Numerical Stability & Precision](numerical-stability.md)
    - [Fast Math Approximations](fast-math.md)
    - [Fourier Transform](fourier-transform.md)
  - [Spatial Data Structures](#spatial-data-structures)
    - [Grids & Rulers](grids-rulers.md)
    - [Spatial Hashing & NN Search](spatial-hashing.md)
- **[Domains](#domains)** — fields of physics & computing where methods are applied
  - [Physics Simulations](#physics-simulations)
    - [Rigid Body Dynamics](rigid-body-dynamics.md)
    - [Soft Body & Truss Dynamics](soft-body-truss-dynamics.md)
    - [Molecular Mechanics & MD](molecular-mechanics.md)
    - [Collision Detection](collision-detection.md)
    - [Aerodynamics & Hydrodynamics](aerodynamics-hydrodynamics.md)
    - [Continuum Mechanics & Impact Fluid](continuum-mechanics-impact.md)
    - [MHD & Plasma Physics](mhd-plasma.md)
    - [Fluid Dynamics (other)](fluid-dynamics.md)
    - [Projective Dynamics](projective-dynamics.md)
    - [Orbital Mechanics & Trajectory Optimization](orbital-mechanics-trajectory.md)
  - [GPU & Parallel Computing](#gpu--parallel-computing)
    - [OpenCL & GPU Computing](opencl-gpu-computing.md)
    - [Parallel Particle-to-Cell](parallel-particle-cell.md)
    - [pyOpenCL Workflows](pyopencl-workflows.md)
  - [Computer Graphics](#computer-graphics)
    - [OpenGL Rendering](opengl-rendering.md)
    - [Radiosity, Ray Tracing & Scattering](radiosity-raytracing-scattering.md)
    - [Solid Modeling & CSG](solid-modeling-csg.md)
    - [Noise & Procedural Generation](noise-procedural.md)
    - [WebGL / WebGPU](webgl-webgpu.md)
  - [Terrain & World Generation](#terrain--world-generation)
    - [Terrain & Hydraulic Erosion](terrain-hydraulics.md)
    - [LandCraft (2D World Sim)](landcraft.md)
- **[Applications](#applications)** — end-user apps, games, and interactive tools
  - [Spacecraft Design, Construction & Combat](spacecraft-design-combat.md)
  - [Combat Models (Melee, Damage, Aero, Land)](combat-models.md)
  - [LandTactics (Tactical Sim)](land-tactics.md)
  - [FormationTactics (Battle Sim)](formation-tactics.md)
  - [Tanks](tanks.md)
  - [Molecular Editor](molecular-editor.md)
- **[Infrastructure](#infrastructure)** — build system, bindings, GUI framework
  - [GUI Components](gui.md)
  - [Python-C/C++ Bindings](python-cpp-bindings.md)
  - [pySymGLSL](pysymglsl.md)
  - [Build System](build-system.md)

## Topic Hierarchy

## Methods & Tools

### Numerical Math

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| Linear Algebra & Iterative Solvers | [linear-algebra-solvers.md](linear-algebra-solvers.md) | created | `cpp/common/math/Lingebra.h`, `cpp/common/math/SparseMatrix.h`, `doc/NumricalStabilityLinearSolvers.md`, `doc/ChebyshevAndInterialAccelerationOfLinearSolvers.md`, `doc/Parallel_GaussSeidel.md`, `python/constrain_solver.py`, `python/PGS.py` |
| Optimization | [optimization.md](optimization.md) | created | `cpp/common/optimization/lineSearch.h`, `cpp/common/optimization/optimizer_random.h`, `cpp/common/math/DynamicOpt.cpp` (FIRE), `cpp/apps/MolecularEditor/doc/GlobalOptimization.md` |
| ODE Integration | [ode-integration.md](ode-integration.md) | created | `cpp/common/dynamics/ODEintegrator.h` (RKF45), `python/pySimE` (various ODE tests) |
| Orbital Mechanics & Trajectory Optimization | [orbital-mechanics-trajectory.md](orbital-mechanics-trajectory.md) | created | `cpp/common/Orbital/SpaceBodies.h`, `cpp/common/Orbital/SpaceWorld.h`, `cpp/common/Orbital/TrajectoryVariation.h`, `cpp/common/dynamics/ODEintegrator.h`, `cpp/common/math/SplineManager.h`, `cpp/common/math/CubicBSpline.h`, `cpp/common/dynamics/DynamicOpt.h`, `cpp/apps/OrbitalWar/test_OptContinuousThrust.cpp`, `projects/SpaceCombat/` |
| Numerical Stability & Precision | [numerical-stability.md](numerical-stability.md) | created | `doc/NumricalStabilityLinearSolvers.md`, `doc/High_Precision_Physics_on_singlepoint_GPU_hardware.md` |
| Fast Math Approximations | [fast-math.md](fast-math.md) | created | `cpp/common/math/fastmath.h`, `cpp/common/math/gonioApprox.h`, `cpp/common/math/functions.h` |
| Fourier Transform | [fourier-transform.md](fourier-transform.md) | created | `cpp/common/math/` (FFT implementation) |

### Spatial Data Structures

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| Grids & Rulers | [grids-rulers.md](grids-rulers.md) | created | `cpp/common/dataStructures/SimplexGrid.h`, `cpp/common/dataStructures/Ruler2DFast.h`, `cpp/common/dataStructures/Grid.h`, `cpp/common/maps/CubicRuler.h`, `docs/LandCraft/Ruler2DFast.md`, `docs/Buckets.md` |
| Spatial Hashing & NN Search | [spatial-hashing.md](spatial-hashing.md) | created | `cpp/common/dataStructures/HashMat.h`, `cpp/common/dataStructures/Buckets.h` |

## Domains

### Physics Simulations

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| Rigid Body Dynamics | [rigid-body-dynamics.md](rigid-body-dynamics.md) | created | `cpp/common/dynamics/Body2D.h`, `cpp/common/dynamics/Body.h`, `cpp/common/dynamics/AeroCraft.h`, `docs/RigidBodyDynamics/RigidBodyProperly.md`, `docs/RigidBodyDynamics/rigid_body_proof.py` |
| Soft Body & Truss Dynamics | [soft-body-truss-dynamics.md](soft-body-truss-dynamics.md) | created | `cpp/common/dynamics/SoftBody.h`, `cpp/common/dynamics/TrussDynamics_d.h`, `python/pyTruss/`, `docs/TrussGeneration/`, `doc/VertexBlockDescent.md`, `doc/Markdown/cpp/common/dynamics/ProjectiveDynamics.md` |
| Molecular Mechanics & MD | [molecular-mechanics.md](molecular-mechanics.md) | created | `cpp/common/dynamics/MMFF.h`, `cpp/apps/MolecularEditor/`, `cpp/apps/MolecularEditor2/`, `python/pyMolecular/`, `doc/AGENTs/protocols/domain/` (forcefields, noncovalent, intramolecular, quantum) |
| Collision Detection | [collision-detection.md](collision-detection.md) | created | `cpp/common/engine/Notes/BroadPhaseCollision.md`, `cpp/common/engine/Notes/AABBTree.md`, `cpp/common/engine/Notes/BoxSweepAndPrune.md`, `cpp/common/dataStructures/HashMat.h`, `java/Common/CellSort.java` |
| Aerodynamics & Hydrodynamics | [aerodynamics-hydrodynamics.md](aerodynamics-hydrodynamics.md) | created | `cpp/common/dynamics/AeroSurf.h`, `cpp/common/dynamics/AeroCraft.h`, `cpp/common/math/PotentialFlow.h`, `cpp/common/CombatModels/AirCombatModel.h`, `cpp/sketches_SDL/3D/test_VortexLattice.cpp`, `docs/BiotSavart/` (4 files) |
| Continuum Mechanics & Impact Fluid | [continuum-mechanics-impact.md](continuum-mechanics-impact.md) | created | `python/EulerianImpacFluid/EulerianImpacFluid.py`, `python/EulerianImpacFluid/EulerianImpacFluid.cl`, `python/EulerianImpacFluid/EulerianImpacFluid.md`, `doc/Continum_Mechanics_Solver.md` |
| MHD & Plasma Physics | [mhd-plasma.md](mhd-plasma.md) | created | `doc/python/MHD/`, `docs/BiotSavart/MHD_Plasma_Nozzle_webgl.md`, `docs/BiotSavart/PotentialFlow.md`, `docs/BiotSavart/VortexParticleMethod.md` |
| Fluid Dynamics (other) | [fluid-dynamics.md](fluid-dynamics.md) | created | `js/FlowField/`, `docs/BiotSavart/AeroPotentialFlow_simulation_javascript.md` |
| Projective Dynamics | [projective-dynamics.md](projective-dynamics.md) | created | `doc/Markdown/cpp/common/OCL/ProjectiveDynamicsOCL.md`, `cpp/common/dynamics/Notes/ProjectiveDynamics.md` |

### GPU & Parallel Computing

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| OpenCL & GPU Computing | [opencl-gpu-computing.md](opencl-gpu-computing.md) | created | `cpp/common/OCL/`, `docs/OpenCLBase.md`, `cpp/common/OCL/OpenCL_API_opt.md`, `python/GLCL2/`, `python/pyRay/`, `python/terrain_ocl/`, `doc/High_Precision_Physics_on_singlepoint_GPU_hardware.md` |
| Parallel Particle-to-Cell | [parallel-particle-cell.md](parallel-particle-cell.md) | created | `doc/Parallel_Particle_To_Cell_accumulation.md` |
| pyOpenCL Workflows | [pyopencl-workflows.md](pyopencl-workflows.md) | created | `python/GLCL2/doc/GLCL_manifest.md`, `python/GLCL2/doc/GLCL_performance_refactor.md`, `doc/CodingRules/workflows/python/python_pyopencl_workflow.md` |

### Computer Graphics

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| OpenGL Rendering | [opengl-rendering.md](opengl-rendering.md) | created | `cpp/common_SDL/SDL2OGL/Draw2D.h`, `cpp/common_SDL/SDL2OGL/Draw3D.h`, `docs/OpenGL_Rendering/` (5 files), `doc/Markdown/GUI.md` |
| Radiosity, Ray Tracing & Scattering | [radiosity-raytracing-scattering.md](radiosity-raytracing-scattering.md) | created | `cpp/common/dynamics/TriangleRayTracer.h`, `cpp/common/dynamics/Radiosity.h`, `cpp/common/dynamics/Scatterer2.h`, `RayTracing.md`, `docs/RadiosityRaytracing/` (4 files), `python/Radiosity/`, `python/pyScatter/`, `cpp/common_resources/cl/radiosity.cl` |
| Solid Modeling & CSG | [solid-modeling-csg.md](solid-modeling-csg.md) | created | `js/GLSL_solid_modeling/`, `python/pyRay/`, `docs/SolidModeling/SolidModeling_js.md` |
| Noise & Procedural Generation | [noise-procedural.md](noise-procedural.md) | created | `cpp/common/math/Noise.h`, `js/LandCraft_web/ShaderToy_inspiration/Noise/`, `js/LandCraft_web/doc/BlueNoise.md` |
| WebGL / WebGPU | [webgl-webgpu.md](webgl-webgpu.md) | created | `doc/js/WebGL_Sim.md`, `doc/js/WebGL_Grok.md`, `doc/js/WebGPU_Chrome_Linux.md`, `js/NBody2D_WebGPU/` |

### Terrain & World Generation

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| Terrain & Hydraulic Erosion | [terrain-hydraulics.md](terrain-hydraulics.md) | created | `cpp/common/maps/TerrainHydraulics.h`, `cpp/common/maps/TerrainSimplex.md`, `python/terrain_ocl/`, `docs/LandCraft/TerrainHydraulics.md`, `docs/LandCraft/hydraulics1D.md`, `docs/LandCraft/BasinFilling.md`, `js/LandCraft_web/` |
| LandCraft (2D World Sim) | [landcraft.md](landcraft.md) | created | `docs/LandCraft/` (15+ files), `js/LandCraft_web/`, `docs/LandCraft/UserGuide.md`, `docs/LandCraft/LandCraft_main.md` |

## Applications

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| Spacecraft Design, Construction & Combat | [spacecraft-design-combat.md](spacecraft-design-combat.md) | created | `cpp/common/Orbital/SpaceCraft.h`, `cpp/common/Orbital/SpaceCraftComponents.h`, `cpp/apps/OrbitalWar/`, `docs/SpaceCrafting/SpaceCrafting_new.md`, `docs/SpaceCrafting/SpaceCraftConstructionProblems.md`, `docs/SpaceTactics/SpaceTactics.md`, `encyclopedia/space_warfare/`, `js/spacecraft_editor/` |
| Combat Models (Melee, Damage, Aero, Land) | [combat-models.md](combat-models.md) | created | `cpp/common/CombatModels/AirCombatModel.h`, `cpp/apps/AeroCombat/`, `cpp/apps/FormationTactics/`, `cpp/apps/LandTactics/`, `cpp/apps/Tanks/`, `docs/FormationTactics/DamagePhysicsModel.md`, `docs/FormationTactics/PolearmsPhysicsModel.md`, `docs/LandTactics/AirCombatModel.md` |
| LandTactics (Tactical Sim) | [land-tactics.md](land-tactics.md) | created | `cpp/apps/LandTactics/`, `docs/LandTactics/` (8 files) |
| FormationTactics (Battle Sim) | [formation-tactics.md](formation-tactics.md) | created | `cpp/apps/FormationTactics/`, `docs/FormationTactics/` (DamagePhysicsModel, PolearmsPhysicsModel) |
| Tanks | [tanks.md](tanks.md) | created | `cpp/apps/Tanks/`, `cpp/apps/Tanks/notes/ArmorPenetration.md` |
| Molecular Editor | [molecular-editor.md](molecular-editor.md) | created | `cpp/apps/MolecularEditor/`, `cpp/apps/MolecularEditor2/`, `cpp/apps_OCL/MolecularEditorOCL/`, `docs/MolGUI_web.md` |

## Infrastructure

| Topic | Audit File | Status | Key Entry Points |
|-------|-----------|--------|-----------------|
| GUI Components | [gui.md](gui.md) | created | `cpp/common_SDL/SDL2OGL/GUI.h`, `python/BaseGUI.py`, `python/GLCL/GLGUI.py`, `doc/Markdown/GUI.md` |
| Python-C/C++ Bindings | [python-cpp-bindings.md](python-cpp-bindings.md) | created | `doc/CodingRules/workflows/python/python_ctypes_workflow.md`, `doc/AGENTs/skills/ctypes-bindings/SKILL.md`, `python/pySimE/` (wraps C++ core) |
| pySymGLSL | [pysymglsl.md](pysymglsl.md) | created | `doc/Markdown/pySymGLSL.md`, `python/pySymGLSL/` |
| Build System | [build-system.md](build-system.md) | created | `cpp/CMakeLists.txt`, `C/build.sh`, `doc/AGENTs/skills/cpp-build/SKILL.md` |

## Source Quality Notes

- **`docs/`** — auto-generated by LLM, may be stale/redundant. Use as secondary reference, verify against source code.
- **`doc/Markdown/cpp/`** — auto-generated API docs from C++ headers. Useful for navigation but not authoritative.
- **`doc/AGENTs/protocols/domain/`** — high-quality domain protocols (forcefields, QM, topology). Authoritative for molecular simulation topics.
- **`doc/AGENTs/skills/`** — high-quality skill definitions. Authoritative for development workflows.
- **`doc/python/`** — LLM chat transcripts and notes for Python simulation modules. Mixed quality.
- **`js/LandCraft_web/ShaderToy_inspiration/*.glsl.md`** — downloaded ShaderToy snippets, not documentation. Exclude from audits.
- **`*_bak/` directories** — backups, exclude from audits.
- **`node_modules/`** — third-party, exclude.

## Audit File Format

See `doc/AGENTs/skills/doc-audit/SKILL.md`. Each audit file uses YAML frontmatter:
```yaml
---
type: TopicalAudit
title: <Topic Name>
tags: [topic, cross-language]
---
```
Body sections: Summary, Implementations (table), Parity Status, Open Issues.
