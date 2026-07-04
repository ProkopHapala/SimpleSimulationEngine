# LandCraft

Economy and logistics simulation game inspired by Transport Tycoon, Cities: Skylines, Rise of Industry, and Factorio. The goal is to efficiently use terrain and natural resources to build civilization. Features procedural terrain generation with hydraulic erosion, road/building construction, vehicle logistics, and economic resource chains.

## Simulation Features

- Procedural terrain with hydraulic erosion simulation
- Road network construction and pathfinding
- Vehicle-based resource transport logistics
- Economic model: resource extraction, processing, trade
- Lua scripting for scenario definitions
- 2D top-down view with tiled rendering

## Files

- **LandCraft_main.cpp** — main application: world view, interaction, building, economy
- **LandCraftWorld.cpp / .h** — core world model: terrain, hydraulics, roads, vehicles, economy
- **Economy.h** — economic model: resources, production, consumption, trade
- **Roads.h** — road network: construction, pathfinding, vehicle routing
- **LandCraftLua.h** — Lua scripting interface for scenario/terrain generation
- **hydraulics1D.h** — 1D hydraulic erosion simulation
- **CMakeLists.txt** — build target: `LandCraft_main`
- **data/** — terrain data, Lua scripts, configuration
- **notes/** — design notes

## Related Documentation

- [docs/LandCraft/LandCraft.md](../../docs/LandCraft/LandCraft.md) — LandCraft overview and architecture
- [docs/LandCraft/UserGuide.md](../../docs/LandCraft/UserGuide.md) — user guide for the LandCraft game
- [docs/LandCraft/Economy.md](../../docs/LandCraft/Economy.md) — economic model: resources, production, trade
- [docs/LandCraft/Roads.md](../../docs/LandCraft/Roads.md) — road network construction and routing
- [docs/LandCraft/PathFinder.md](../../docs/LandCraft/PathFinder.md) — pathfinding algorithm
- [docs/LandCraft/TerrainHydraulics.md](../../docs/LandCraft/TerrainHydraulics.md) — hydraulic erosion simulation
- [docs/LandCraft/TerrainCubic.md](../../docs/LandCraft/TerrainCubic.md) — cubic-spline terrain representation
- [docs/LandCraft/hydraulics1D.md](../../docs/LandCraft/hydraulics1D.md) — 1D hydraulics model
- [docs/LandCraft/TiledView.md](../../docs/LandCraft/TiledView.md) — tiled view rendering system
- [docs/LandCraft/Ruler2DFast.md](../../docs/LandCraft/Ruler2DFast.md) — fast 2D ruler/grid data structure
- [docs/LandCraft/BasinFilling.md](../../docs/LandCraft/BasinFilling.md) — basin filling / water accumulation simulation
- [docs/LandCraft/Land2D_simulation_javascript.md](../../docs/LandCraft/Land2D_simulation_javascript.md) — JavaScript web version of the simulation
- [docs/LandCraft/LandCraft_main.md](../../docs/LandCraft/LandCraft_main.md) — main application design notes
