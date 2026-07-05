---
name: spacecraft-debug
description: Navigate and debug OrbitalWar / SpaceCraft — design, truss sim, combat, radiosity, orbits; hub to docs and tests
trigger:
  glob:
    - "**/Orbital/**"
    - "**/OrbitalWar/**"
    - "**/spacecraft_editor/**"
    - "**/SpaceCraft*"
    - "**/EditSpaceCraft*"
    - "**/ship_*.lua"
    - "**/tests_bash/Orbital/**"
    - "**/docs/SpaceCrafting/**"
---

## What This System Is

**OrbitalWar** is SSE’s integrated stack for **designing, building, and fighting with spacecraft** in a physically motivated model — not a single mesh format or editor, but a pipeline from ship definition to combat outcomes.

**Problems it addresses:**

- **Structural design** — ships as typed components (girders, rings, ropes, sliders, shields, radiators, tanks, thrusters, guns), authored in **Lua** (and optionally sketch geometry), held in **`SpaceCraft` as SSOT**
- **Geometry** — convert logical structure to **`Mesh::Builder2`** truss/block meshes for rendering, export, and simulation
- **Elastic dynamics** — **`TrussDynamics`** (CPU/OpenCL/Python): bonds, projective dynamics, sliders as actuators, collision impulses; topology changes when damage breaks bonds (partially implemented)
- **Thermal / radiation** — **`Radiosity`**, triangle ray-tracing, Python/OpenCL scattering and X-ray attenuation on ship surfaces
- **Orbital mechanics & trajectories** — N-body worlds, continuous-thrust optimization, trajectory splines and time scrubbing
- **Space combat** — railguns, lasers, nuclear weapons, Whipple shields; hit testing and damage coupling to the truss model
- **Tooling** — C++ apps (editor, dynamics, mesh export), JS web editor, headless `tests_bash` regression scripts

Debugging usually means tracing **which layer** failed: logical model → mesh build → truss export → solver → weapons/thermal/orbit coupling. Details live in linked docs and file-header **caveats** — this skill routes you there.

## Pipeline (mental model)

```
Ship definition (Lua / editor / sketch import)
        → SpaceCraft (components + materials)
        → Mesh::Builder2 (blocks / truss sticks)
        → TrussDynamics (+ sliders, damage, collision)
        → coupled: Radiosity / ray hits / N-body orbit / weapons
```

Parallel **JavaScript** editor (`js/spacecraft_editor/`) and **Python** truss/radiosity solvers mirror parts of the C++ path — check parity when symptoms differ by front-end.

## Module Router (C++ — start here)

| Concern | Where |
|--------|--------|
| SSOT container, deferred import, `fulfillTag` | `cpp/common/Orbital/SpaceCraft.h` |
| Components, `sideToPath`, sliders, plates | `cpp/common/Orbital/SpaceCraftComponents.h` |
| Lua API (`Girder`, `fromObj`, …) | `cpp/common/Orbital/EditSpaceCraft.h` |
| OBJ sketch import, `BuildCraft_sketch` | `cpp/common/Orbital/SpaceCraftFromOBJ.h` |
| Full mesh + sliders (editor/dynamics) | `cpp/common/Orbital/SpaceCraft2Mesh2.h` |
| CLI / headless truss export | `cpp/common/Orbital/SpaceCraft2Mesh_blocks.h` |
| Flat mesh arrays | `cpp/common/geometry/MeshBuilder2.h` |
| Truss physics, bond strain, solvers | `cpp/common/dynamics/TrussDynamics_d.h` |
| Radiosity, surface coupling | `cpp/common/dynamics/Radiosity.h`, `TriangleRayTracer.h` |
| Space combat, weapons, shields | `cpp/common/Orbital/spaceCombat.h` |
| N-body world, orbits, combat eval | `cpp/common/Orbital/SpaceWorld.h` |
| Trajectory / continuous thrust | `cpp/apps/OrbitalWar/test_OptContinuousThrust.cpp`, `optContinuousThrust.sh` |
| Space combat game | `cpp/apps/OrbitalWar/spaceTactics.cpp` |
| Solar system / orbit editor | `cpp/apps/OrbitalWar/SolarSystemMap.cpp`, `orbitEditor.cpp` |
| Headless CLI entry | `cpp/apps/OrbitalWar/spaceCraftMeshExport.cpp` |
| Interactive editor / sim apps | `cpp/apps/OrbitalWar/spaceCraftEditor*.cpp`, `spaceCraftDynamics*.cpp` |
| Example ship scripts | `cpp/apps/OrbitalWar/data/ship_*.lua`, `sketch_*.obj` |
| App index | `cpp/apps/OrbitalWar/README.md` |

**Before changing code:** read target file header + `doc-read-navigate` skill.

## Docs Router

| Topic | Document |
|-------|----------|
| Full codebase map (C++ + JS + Python) | `docs/SpaceCrafting/Codebase_Reference_Map.md` |
| Pipeline overview | `docs/SpaceCrafting/SpaceCrafting_new.md` |
| Physical / combat model (encyclopedia) | `encyclopedia/space_warfare/` |
| Mesh export CLI | `docs/SpaceCrafting/spaceCraft_mesh_export_cli.md` |
| LOD / sketch import (subset) | `docs/SpaceCrafting/SpaceCRaft_LODs.chat.md` |
| Construction pitfalls | `docs/SpaceCrafting/SpaceCraftConstructionProblems.md` |
| JS / web editor | `docs/SpaceCrafting/SpaceCraft_web_new.md` |
| Python truss / radiosity | `python/pyTruss/`, `python/pyScatter/`, `python/Radiosity/` |

## Running Tests & Demos

Scripts live in `tests_bash/Orbital/`. Pattern (see `tests_bash/README.md`):

1. `cd tests_bash/Orbital`
2. Script rebuilds via `cpp/Build/apps/OrbitalWar`, symlinks `*.x`, sets `LD_PRELOAD` to ASan
3. Symlinks `data` → `cpp/apps/OrbitalWar/data`

| Script | Use when debugging |
|--------|-------------------|
| `spaceCraftMeshExport.sh` | Headless Lua → OBJ + truss |
| `spaceCraftDynamics.sh` | Truss sim + weapons + damage (GUI) |
| `spaceTactics.sh` | Full space combat + N-body + time scrubbing |
| `optContinuousThrust.sh` | Trajectory optimization with thrust |
| `spaceCraftEditor.sh` / `SpaceCraftEditorNew.sh` | Interactive Lua ship load / edit |
| `constructionBlock.sh` | Mesh generators (parabola, quad slab, …) |
| `trussSimBatch.sh` | Batch truss from exported `.truss` |
| `spaceCraftDynamicsOCL.sh` | GPU dynamics path |
| `SolarSystemMap.sh` / `orbitEditor` apps | Orbit / system map tooling |

**Fast headless check** (after build):

```bash
cd tests_bash/Orbital
./spaceCraftMeshExport.x -s data/ship_ICF_marksman_2.lua -o /tmp/ship.obj -v 1
```

## Common Debug Themes (see file-header caveats)

- **Mesh build** — `sideToPath` / slider rails need tessellated host; `build_order` not fully wired; two `BuildCraft_blocks` variants (CLI vs editor)
- **Truss / damage** — bond breaking TODO in `TrussDynamics`; Cholesky invalid after topology change; use iterative solvers when debugging damage
- **Combat / hits** — ray-triangle and projected targets in `spaceCombat.h`; coupling to truss damage still incomplete
- **Radiosity / thermal** — coupling matrix + occlusion work; integration with live sim varies by app
- **Cross-language** — C++ vs JS editor vs Python solvers may diverge; use `numerical-parity` skill
- **Sketch / OBJ import** (optional path) — deferred fields, `fulfillTag`, sketch vs blocks LOD — see `SpaceCraftFromOBJ.h`

## Related Skills

- `cpp-build` — cmake, ASAN vs opt, `cpp/Build` symlink
- `cpp-memory` — ASAN crashes in mesh/truss code
- `visual-debug` — SVG export, review artifacts (`.obj`/`.svg` paths)
- `numerical-parity` — C++ vs JS or ref mesh counts
- `code-reuse` / `doc-read-navigate` — before adding parallel implementations
- `doc-task-summary` — after fixing a bug, update header caveats
