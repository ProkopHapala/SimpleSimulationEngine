---
type: TopicalAudit
title: Spacecraft Design, Construction & Combat
tags: [topic, cross-language, spacecraft, orbital, combat, games]
---

## Summary

Spacecraft design pipeline: Lua scripts → `SpaceCraft` component model → `Mesh::Builder2` truss/block mesh → truss dynamics, radiosity/scattering, damage/hit-testing. Multiple apps: `OrbitalWar` (space combat sim), `spaceCraftEditor` (interactive design), `spaceCraftDynamics` (physics sim), `spaceCraftDynamicsOCL` (GPU). Web versions in JS. Encyclopedia of space warfare physics. Python Jupyter notebooks for combat calculations.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/Orbital/SpaceCraft.h` | active | Central ship model: nodes, girders, ropes, rings, sliders, welds, radiators, shields, collectors, tanks, thrusters, guns |
| C++ | `cpp/common/Orbital/SpaceCraftComponents.h` | active | Component hierarchy: `CatalogItem` → `ShipComponent` → `StructuralComponent` → specific parts |
| C++ | `cpp/common/Orbital/EditSpaceCraft.h` | active | C API for Lua scripting of spacecraft design |
| C++ | `cpp/common/Orbital/SpaceCraft2Truss.h` | active | Convert spacecraft to truss mesh |
| C++ | `cpp/common/Orbital/SpaceCraft2Mesh2.h` | active | Older truss-style mesh conversion |
| C++ | `cpp/common/Orbital/SpaceCraft2Mesh_blocks.h` | active | Newer block-based mesh conversion via `Mesh::Builder2` |
| C++ | `cpp/common/Orbital/spaceCraftEditorUtils.h` | active | Editor utilities |
| C++ | `cpp/apps/OrbitalWar/spaceTactics.cpp` | active | Space combat game: N-body physics, trajectory splines, weapons (railguns, lasers, Whipple shields), time scrubbing |
| C++ | `cpp/apps/OrbitalWar/spaceCraftEditor.cpp` | active | Interactive spacecraft editor |
| C++ | `cpp/apps/OrbitalWar/spaceCraftDynamics.cpp` | active | Spacecraft dynamics simulation |
| C++ | `cpp/apps/OrbitalWar/spaceCraftDynamicsOCL.cpp` | experimental | GPU-accelerated dynamics |
| C++ | `cpp/common/Orbital/SpaceWorld.h` | active | N-body gravitational simulation, RKF45 integrator, planets + ships |
| JS | `js/spacecraft_editor/` | active | Web-based spacecraft editor |
| JS | `js/spacecraft_editor/multiagent_plan.md` | active | Multi-agent development plan for web editor |
| Python | `projects/SpaceCombat/` | active | Jupyter notebooks: combat calculations, maneuverability, industry models |
| Doc | `docs/SpaceCrafting/SpaceCrafting_new.md` | active | High-level overview: Lua → mesh → simulation pipeline (35KB, comprehensive) |
| Doc | `docs/SpaceCrafting/SpaceCraftConstructionProblems.md` | active | Detailed construction problems: telescopic truss recoil damper, plates, welds, SDF selection (49KB, WIP) |
| Doc | `docs/SpaceCrafting/SpaceCraft_web.md` | active | Web spacecraft editor documentation |
| Doc | `docs/SpaceCrafting/SpaceCraft_web_new.md` | active | Newer web editor doc |
| Doc | `docs/SpaceCrafting/SpaceCraftingWithBlocks.md` | active | Block-based construction system |
| Doc | `docs/SpaceCrafting/SpaceCraftingWithBlocks_new.md` | active | Newer block construction |
| Doc | `docs/SpaceCrafting/SpaceCraft_simulation_javascript.md` | active | JS simulation documentation |
| Doc | `docs/SpaceCrafting/SpaceCraft_gen_web_2025_12_29.md` | active | Web generation doc |
| Doc | `docs/SpaceCrafting/spaceCraft_mesh_export_cli.md` | active | CLI mesh export documentation |
| Doc | `docs/SpaceCrafting/SpaceCrafting_cpp_devnotes.md` | active | C++ development notes |
| Doc | `docs/SpaceTactics/SpaceTactics.md` | active | Comprehensive code map: N-body, weapons, propulsion, shields |
| Doc | `docs/SpaceTactics/SpaceBodies.md` | active | Celestial body definitions |
| Doc | `docs/SpaceTactics/SpaceWorld.md` | active | World model documentation |
| Doc | `docs/SpaceTactics/appliedPhysics.md` | active | Applied physics models |
| Doc | `docs/SpaceTactics/asteroidEngineering.md` | active | Asteroid mining/engineering |
| Doc | `docs/SpaceTactics/spaceCombat.md` | active | Combat model documentation |
| Doc | `docs/SpaceTactics/SplineManager.md` | active | Trajectory spline management |
| Doc | `doc/OrbitalWar_wiki.md` | active | Wiki-style overview |
| Doc | `doc/Markdown/cpp/common/Orbital/SpaceCraft.h.md` | active | Auto-generated API doc |
| Doc | `doc/Markdown/cpp/common/Orbital/SpaceCraftComponents.h.md` | active | Auto-generated API doc |
| Encyclopedia | `encyclopedia/space_warfare/00_Overview.md` | active | Space warfare overview |
| Encyclopedia | `encyclopedia/space_warfare/01_Tactics.md` | active | Tactical doctrine |
| Encyclopedia | `encyclopedia/space_warfare/02_Weapons.md` | active | Weapon systems: kinetics, lasers, nukes |
| Encyclopedia | `encyclopedia/space_warfare/03_Defense.md` | active | Defense: Whipple shields, armor |
| Encyclopedia | `encyclopedia/space_warfare/04_Propulsion.md` | active | Propulsion: chemical, nuclear pulse, electric |
| Encyclopedia | `encyclopedia/space_warfare/05_Ship_construction.md` | active | Ship construction principles |
| Encyclopedia | `encyclopedia/space_warfare/07_Physical_model_and_equations.md` | active | Physical equations: damage, radiation, heat transport |

## Parity Status

- **C++ spacecraft model ↔ JS web editor**: Both implement same component hierarchy. JS is newer/active development per construction problems doc. C++ is the reference.
- **C++ `spaceCraftDynamics` ↔ `spaceCraftDynamicsOCL`**: CPU vs GPU dynamics. OCL version experimental.
- **Truss dynamics**: spacecraft mesh feeds into `TrussDynamics_d` (see [soft-body-truss-dynamics.md](soft-body-truss-dynamics.md))
- **Radiosity/scattering**: spacecraft surfaces feed into radiosity and scattering solvers (see [radiosity-raytracing-scattering.md](radiosity-raytracing-scattering.md))
- **Damage model**: projectile clouds as line segments, intersected with mesh via `TriangleRayTracer`. Nuclear/radiative damage reuses radiosity/scattering tools.

## Open Issues

- Telescopic truss recoil damper (P1) — main focus per construction problems doc, hexagonal outer + triangular inner tube
- Plates on girders/ropes (P2) — triangular/quad plates between polylines, tessellation matching
- Welds and automatic connections (P3) — weld manager, rigid vs soft connections
- `spaceCraftDynamicsOCL` status unclear — experimental GPU version
- Multiple `docs/SpaceCrafting/` files overlap significantly (Bak-GPT5, Bak-K2 versions are backups)
- JS editor has many completed tasks (J0-J8) but J2, J7 still pending
- Spacecraft mesh export CLI documented but integration with simulation pipeline incomplete
- See `docs/SpaceCrafting/SpaceCraftConstructionProblems.md` for full TODO checklist
