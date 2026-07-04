# cpp/apps — Applications & Games

Interactive SDL2/OpenGL applications built on the SSE engine. Each subdirectory is a separate app with its own `CMakeLists.txt`. Build via CMake with `WITH_SDL=ON`. Run via `tests_bash/apps/*.sh` scripts.

## Application Index

### Combat & Tactics

| App | Description | Key Physics |
|-----|-------------|-------------|
| [AeroCombat](AeroCombat/README.md) | Flight simulator with realistic aerodynamics and air combat | AeroSurf, Body.h, terrain collision |
| [FormationTactics](FormationTactics/README.md) | TotalWar-like historical battle simulation with formations, morale, melee | Soldier AI, formation cohesion, faction AI |
| [LandTactics](LandTactics/README.md) | WWII tactical land battle with terrain cover, line-of-sight | LOS, terrain mobility, armor/penetration |
| [MinimalTactics](MinimalTactics/README.md) | Minimal tactical combat for testing simplified combat models | Basic unit combat |
| [StrategicBattleModels](StrategicBattleModels/README.md) | Strategic-level battle model definitions (headers only, WIP) | Armor penetration, damage model |
| [Tanks](Tanks/README.md) | Tank sim with turret control, spline terrain, armor penetration | Armor per-polygon, gun ports |
| [MultiFight3D](MultiFight3D/README.md) | 3D multiplayer melee and projectile combat arena | Rigid body, collision detection |
| [SwordPlay](SwordPlay/README.md) | Soft-body sword fighting (Mortal Kombat / Blade Symphony style) | SoftPolyLine, muscle-skeleton |
| [NonInertial](NonInertial/README.md) | Liero/Worms-like combat in rotating non-inertial frame | Coriolis, centrifugal forces |

### Naval

| App | Description | Key Physics |
|-----|-------------|-------------|
| [SailWar](SailWar/README.md) | Age-of-sail naval combat with sailing aerodynamics | Wind on sails, buoyancy, cannon ballistics |
| [SailWar_Multi](SailWar_Multi/README.md) | Multiplayer SailWar with UDP networking (requires `WITH_NET=ON`) | Same + SDL2_net state sync |
| [NavalBattle](NavalBattle/README.md) | Battleship with multiple cannon turrets | Buoyancy, turret rotation, projectile ballistics |

### Space

| App | Description | Key Physics |
|-----|-------------|-------------|
| [OrbitalWar](OrbitalWar/README.md) | Spacecraft design, construction, and combat. The flagship app. | N-body gravity, truss dynamics, weapons, damage, OpenCL |

### Building & Creation

| App | Description | Key Physics |
|-----|-------------|-------------|
| [LandCraft](LandCraft/README.md) | Economy/logistics sim (Transport Tycoon / Factorio style) | Terrain hydraulics, road pathfinding, economy |
| [BlockHouseTactics](BlockHouseTactics/README.md) | Tactical blockhouse construction from 3D-hashed boxes | 3D spatial hashmap, box faces |
| [CastleBuilder](CastleBuilder/README.md) | Castle/building construction from modular components | Modular building blocks |
| [CAD](CAD/README.md) | 2D CAD tool and analytical geometry demo | 2D geometric constraints |
| [ShapePainter](ShapePainter/README.md) | Drawing program with parametric shape brushes | Hybrid pixel/vector, parametric shapes |

### Science & Visualization

| App | Description | Key Physics |
|-----|-------------|-------------|
| [MolecularEditor](MolecularEditor/README.md) | Interactive molecular editor and MD simulator | MMFF force field, FIRE optimization, conformational search |
| [MusicVizualizer](MusicVizualizer/README.md) | Real-time music visualizer with OpenGL 3 shaders (requires `WITH_MUSIC=ON`) | FFT, shader effects |
| [DemoCrat](DemoCrat/README.md) | Plugin-based demo loader using `dlopen`/`dlsym` | Dynamic shared library loading |

### Archived

| App | Description |
|-----|-------------|
| [temp](temp/README.md) | Legacy/experimental apps: Fortress, NBodyWorld, Ship2D, SpaceCraft, SailWar_h_bak |

## Build

All apps require `WITH_SDL=ON` in CMake. Additional options:
- `WITH_NET=ON` — enables SailWar_Multi (requires SDL2_net)
- `WITH_MUSIC=ON` — enables MusicVizualizer (requires SDL2_mixer)
- `WITH_OPENCL=ON` — enables GPU-accelerated OrbitalWar dynamics

## Run

Use the shell scripts in `tests_bash/apps/` to compile and run each app:
```bash
cd tests_bash/apps
./Tanks.sh        # compiles and runs Tanks
./AeroCombat.sh   # compiles and runs AeroCombat
```

## Documentation

Detailed design documents for many apps live in the [docs/](../../docs/) directory:

- [docs/AeroCombat/](../../docs/AeroCombat/) — aerodynamic combat model, AeroSurf force calculation
- [docs/FormationTactics/](../../docs/FormationTactics/) — damage physics model, polearms physics model
- [docs/LandCraft/](../../docs/LandCraft/) — LandCraft architecture, economy, roads, hydraulics, terrain, pathfinding, user guide
- [docs/LandTactics/](../../docs/LandTactics/) — LandTactics world model, units, factions, shelters, air combat
- [docs/SpaceCrafting/](../../docs/SpaceCrafting/) — spacecraft crafting system, block construction, mesh export, web versions
- [docs/SpaceTactics/](../../docs/SpaceTactics/) — space tactics game design, combat mechanics, physics, trajectory splines
- [docs/TrussGeneration/](../../docs/TrussGeneration/) — Mesh::Builder2, construction blocks, truss conversion, UV mapping
- [docs/RigidBodyDynamics/](../../docs/RigidBodyDynamics/) — rigid body dynamics formulation and proof
- [docs/RadiosityRaytracing/](../../docs/RadiosityRaytracing/) — radiosity, scattering, sphere sampling, ray tracing
- [docs/BiotSavart/](../../docs/BiotSavart/) — potential flow, vortex particle method, MHD plasma nozzle
- [docs/SolidModeling/](../../docs/SolidModeling/) — solid modeling with GLSL (CSG)
- [docs/OpenGL_Rendering/](../../docs/OpenGL_Rendering/) — OpenGL 3 rendering, mesh rendering, text rendering
- [docs/MolGUI_web.md](../../docs/MolGUI_web.md) — molecular GUI web version
- [docs/Buckets.md](../../docs/Buckets.md) — Buckets data structure (used in mesh building)
- [docs/OpenCLBase.md](../../docs/OpenCLBase.md) — OpenCL base class for GPU acceleration
