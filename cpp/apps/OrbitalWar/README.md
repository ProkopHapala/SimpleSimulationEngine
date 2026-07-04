# OrbitalWar

Spacecraft design, construction, and combat simulation. The flagship application of SSE. Design spacecraft from components (girders, plates, thrusters, guns, shields, radiators) using Lua scripts, convert to truss mesh, simulate dynamics (N-body gravity + truss elasticity + damage), and fight in space combat with railguns, lasers, and nuclear weapons.

## Features

- **Spacecraft design**: Lua-scripted component model → `Mesh::Builder2` truss/block mesh
- **N-body orbital mechanics**: planets, moons, ships with RKF45 integration
- **Truss dynamics**: spacecraft mesh simulated as elastic truss (SoftBody/TrussDynamics)
- **Weapons**: railguns, lasers, nuclear warheads, Whipple shields
- **Damage**: projectile-mesh intersection via ray tracing, radiative damage via scattering
- **Trajectory optimization**: continuous thrust trajectory planning
- **Time scrubbing**: replay simulation forward/backward
- **GPU acceleration**: OpenCL version of spacecraft dynamics
- **Mesh export**: export spacecraft mesh to `.obj`/`.truss` files

## Files

### Main applications

- **spaceTactics.cpp** — space combat game: N-body physics, trajectory splines, weapons, time scrubbing
- **spaceCraftEditor.cpp** — interactive spacecraft editor with Lua ship loading
- **SpaceCraftEditorNew.cpp** — newer version of the spacecraft editor
- **spaceCraftDynamics.cpp** — spacecraft dynamics simulation (truss physics, weapons, damage)
- **spaceCraftDynamicsOCL.cpp** — GPU-accelerated (OpenCL) spacecraft dynamics
- **spaceCraftMeshExport.cpp** — headless mesh export to `.obj`/`.truss`
- **constructionBlockApp.cpp** — block-based spacecraft construction demo
- **SolarSystemMap.cpp** — solar system map viewer
- **orbitEditor.cpp** — orbital trajectory editor
- **trussSimBatch.cpp** — batch truss simulation (headless)

### Test applications

- **test_OptContinuousThrust.cpp** — trajectory optimization with continuous thrust
- **test_SpaceFlightODE.cpp** — spacecraft flight ODE integration test

### Headers

- **SpaceCraftDynamicsApp.h** — spacecraft dynamics application framework
- **spaceCraftSimulator.h** — CPU spacecraft simulator
- **spaceCraftSimulatorOCL.h** — OpenCL spacecraft simulator

### Other

- **CMakeLists.txt** — build targets for all applications
- **data/** — Lua ship scripts (22 files: ship_ICF_*, ship_NTR_*, ship_NFPP_*, colony_*)
- **Notes/** — design notes
- **py/** — Python analysis scripts
- **consrain_solver_testing.txt** — constraint solver testing notes

## Related Documentation

### Spacecraft Design & Construction

- [docs/SpaceCrafting/SpaceCrafting.md](../../docs/SpaceCrafting/SpaceCrafting.md) — spacecraft crafting system overview
- [docs/SpaceCrafting/SpaceCrafting_new.md](../../docs/SpaceCrafting/SpaceCrafting_new.md) — newer version of the crafting design
- [docs/SpaceCrafting/SpaceCraftingWithBlocks.md](../../docs/SpaceCrafting/SpaceCraftingWithBlocks.md) — block-based spacecraft construction
- [docs/SpaceCrafting/SpaceCraftingWithBlocks_new.md](../../docs/SpaceCrafting/SpaceCraftingWithBlocks_new.md) — newer block-based construction
- [docs/SpaceCrafting/SpaceCraftConstructionProblems.md](../../docs/SpaceCrafting/SpaceCraftConstructionProblems.md) — construction problems and solutions
- [docs/SpaceCrafting/SpaceCrafting_cpp_devnotes.md](../../docs/SpaceCrafting/SpaceCrafting_cpp_devnotes.md) — C++ development notes
- [docs/SpaceCrafting/spaceCraft_mesh_export_cli.md](../../docs/SpaceCrafting/spaceCraft_mesh_export_cli.md) — mesh export CLI documentation
- [docs/SpaceCrafting/SpaceCraft_web.md](../../docs/SpaceCrafting/SpaceCraft_web.md) — web version of spacecraft editor
- [docs/SpaceCrafting/SpaceCraft_web_new.md](../../docs/SpaceCrafting/SpaceCraft_web_new.md) — newer web version
- [docs/SpaceCrafting/SpaceCraft_simulation_javascript.md](../../docs/SpaceCrafting/SpaceCraft_simulation_javascript.md) — JavaScript simulation version

### Space Tactics & Combat

- [docs/SpaceTactics/SpaceTactics.md](../../docs/SpaceTactics/SpaceTactics.md) — space tactics game design overview
- [docs/SpaceTactics/spaceCombat.md](../../docs/SpaceTactics/spaceCombat.md) — space combat mechanics: weapons, shields, damage
- [docs/SpaceTactics/SpaceWorld.md](../../docs/SpaceTactics/SpaceWorld.md) — space world model: bodies, physics
- [docs/SpaceTactics/SpaceBodies.md](../../docs/SpaceTactics/SpaceBodies.md) — celestial body definitions
- [docs/SpaceTactics/appliedPhysics.md](../../docs/SpaceTactics/appliedPhysics.md) — applied physics models for space combat
- [docs/SpaceTactics/SplineManager.md](../../docs/SpaceTactics/SplineManager.md) — trajectory spline management
- [docs/SpaceTactics/asteroidEngineering.md](../../docs/SpaceTactics/asteroidEngineering.md) — asteroid engineering concepts

### Truss Generation & Mesh

- [docs/TrussGeneration/MeshBuilder2.md](../../docs/TrussGeneration/MeshBuilder2.md) — Mesh::Builder2 dynamic mesh builder documentation
- [docs/TrussGeneration/ConstructionBlock.md](../../docs/TrussGeneration/ConstructionBlock.md) — construction block system for spacecraft
- [docs/TrussGeneration/truss_low_to_hi.md](../../docs/TrussGeneration/truss_low_to_hi.md) — low-res to high-res truss conversion
- [docs/TrussGeneration/MeshFaceGenerationDiscrepancy.md](../../docs/TrussGeneration/MeshFaceGenerationDiscrepancy.md) — mesh face generation issues
- [docs/TrussGeneration/DrawUV.md](../../docs/TrussGeneration/DrawUV.md) — UV mapping for mesh rendering
