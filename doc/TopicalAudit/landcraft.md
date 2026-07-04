---
type: TopicalAudit
title: LandCraft
tags: [topic, cpp, javascript, webgpu, terrain, hydraulics, economy, rivers, roads, vehicles, pathfinder, simulation, game]
---

## Summary

LandCraft is a terrain simulation and world-building application with both C++ (SDL2+OpenGL) and JS (WebGPU) implementations. C++ version: `LandCraftWorld` holds terrain/hydraulics, rivers, roads, vehicles, economy, and pathfinder. Terrain generation via bisection noise + droplet erosion. Water relaxation via raster scanning and contour-based basin filling (Bellman-Ford, Dijkstra, priority-flood). River finding via rain accumulation. Road building with profile extraction. Economy with technology chains and factory production. Config-driven initialization via `world.ini`. JS version: WebGPU-based terrain generation with noise/FBM, colorization, and viewport rendering.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/apps/LandCraft/LandCraftWorld.h` | active | Central world state: `HydraulicGrid2D hydro`, roads, vehicles, economy (techs, factories), `PathFinder pf`. Config-driven init via `loadConfig()`. |
| C++ | `cpp/apps/LandCraft/LandCraftWorld.cpp` | active | Implementation: terrain generation (bisecNoise + erosion), river finding, road building, vehicle spawning, economy, pathfinder binding, config parsing. |
| C++ | `cpp/apps/LandCraft/LandCraft_main.cpp` | active | `LandCraftApp(AppSDL2OGL)`: GUI, rendering, event handling, command dispatch, draw layers, hydro relaxation, vehicle update. |
| C++ | `cpp/apps/LandCraft/Economy.h` | active | `Commodity`, `Technology` (consumes/produces maps), `Factory` (produce with resource checking). `str2map()` parser. |
| C++ | `cpp/common/maps/TerrainHydraulics.h` | active | `HydraulicGrid2D(Grid2DAlg)`: ground/water arrays, lakes, rivers, droplet erosion, water relaxation (raster X/Y), basin filling (Bellman-Ford, Dijkstra, contour two-wave), `findAllRivers()`, `traceDroplet()`. See `fluid-dynamics.md`. |
| C++ | `cpp/common/maps/HydraulicGrid2D.h` | active | Hydraulic grid extensions. See `fluid-dynamics.md`. |
| C++ | `cpp/libs/CombatModels/LandCraftLib.cpp` | active | C API wrapper for LandCraftWorld: terrain, hydraulics, basin filling, pathfinder, rivers, roads, vehicles, economy. Used by Python ctypes. See `python-cpp-bindings.md`. |
| C++ | `cpp/apps/LandCraft/Roads.h` | active | `Road`, `RoadBuilder`, `RoadVehicleType`, `RoadVehicle` — road path building and vehicle movement. |
| C++ | `cpp/common/maps/PathFinder.h` | active | Path finding on terrain grid. `set()`, `bind()`, `allocate()`. |
| JS+WGSL | `js/LandCraft_web/` | active | WebGPU terrain generation. See `webgl-webgpu.md`, `noise-procedural.md`. |
| JS+WGSL | `js/LandCraft_web/main.js` | active | Main loop: generation pass, colorize pass, render pass, readback. Progressive accumulation. |
| JS+WGSL | `js/LandCraft_web/generators.js` | active | WGSL terrain generators (FBM, ridged multifractal). See `noise-procedural.md`. |
| JS+WGSL | `js/LandCraft_web/noise-lib.js` | active | WGSL noise functions. See `noise-procedural.md`. |
| JS+WGSL | `js/LandCraft_web/gpu-core.js` | active | `GPUContext`: WebGPU device, canvas, uniforms, camera. See `webgl-webgpu.md`. |
| JS+WGSL | `js/LandCraft_web/compute.js` | active | `MapAlgorithm`: compute pipeline. See `webgl-webgpu.md`. |
| JS+WGSL | `js/LandCraft_web/renderers.js` | active | `ViewportRenderer`, `LineRenderer`, `TextRenderer`. See `webgl-webgpu.md`. |
| JS+WGSL | `js/LandCraft_web/LandCraft.html` | active | Entry point HTML. |
| JS+WGSL | `js/LandCraft_web/chemistry.html` | active | Chemistry simulation page. |
| JS+WGSL | `js/LandCraft_web/data/` | active | Data files (CSV). |
| JS+WGSL | `js/LandCraft_web/doc/` | active | Documentation (FractalTerrain.md). |
| Python | `tests_bash/LandCraft/test_landcraft.py` | active | Python test harness via ctypes C API: terrain generation, basin filling (Bellman/Dijkstra/contour), save/load, plotting. |
| Bash | `tests_bash/LandCraft/LandCraft.sh` | active | Shell test runner. |
| Bash | `tests_bash/LandCraft/LandCraftLib.sh` | active | Library test runner. |
| Doc | `docs/LandCraft/Land2D_simulation_javascript.md` | doc | JS simulation design. |
| Doc | `docs/LandCraft/BasinFilling.md` | doc | Basin filling algorithms. |
| Doc | `docs/LandCraft/Economy.md` | doc | Economy model design. |
| Doc | `docs/LandCraft/LandCraft.md` | doc | General LandCraft design. |
| Doc | `js/LandCraft_web/doc/FractalTerrain.md` | doc | Terrain generation algorithms. See `noise-procedural.md`. |

## Sub-topics

### Terrain Generation

C++ (`LandCraftWorld::generateTerrain()`):
1. `makeTerrainBisec(seed)` — bisection noise on power-of-2 grid
2. `droplerErosion(niter, nDrops, nStepMax, ...)` — hydraulic erosion via droplet simulation
3. `normalize_ground(maxHeight)` — scale to [0, maxHeight], water = ground

JS (WebGPU `generators.js`):
1. FBM with configurable octaves, noise flavor (sin/value/simplex/gradient/orbit)
2. Domain warping, ridged multifractal
3. Derivative-based attenuation for elevated terrain
4. Colorization pass (height-based color mapping)

### Water Relaxation & Basin Filling

C++ (`HydraulicGrid2D`):
- `relaxWaterRasterX/Y()` — scan along rows/columns, find water blocks, flatten to average level
- `basinFill_BellmanFord_Boundary()` — iterative relaxation from boundary seeds
- `basinFill_Dijkstra_Boundary()` — priority-queue based filling
- `basinFill_Contour_beginFlood/floodStep/beginDrain/drainStep/finish()` — two-wave contour method (flood + drain)
- `sumWater()` — conservation check (total water, potential energy)

### River Finding

- `gatherRain(minSinkFlow)` — accumulate flow from sinks
- `findAllRivers(minFlow)` — trace rivers from high-flow cells to outlets
- `River` class: `path` (cell indices), `flow` (per-cell flow), `mouth` (parent river)
- Rivers sorted by length in `makeRivers()`

### Road Building & Vehicles

- `makeRoadStraight(p1, p2)` — `RoadBuilder::pushStright()` generates path
- `roadProfileHeights()` — sample ground/water along road
- `RoadVehicle` — moves along road path, `ipath`/`idir`/`onWay` state
- `vehicleStepAll(dt)` — advance all vehicles

### Economy

- `Technology` — name, cycle_time, unit_space, consumes/produces maps
- `Factory` — `stored` map, `setTechnology()`, `produce(N)` — checks resources, computes max production, consumes/produces
- `loadTechnologies(fname)` — parse tech file

### Pathfinder

- `PathFinder::set(hydro)` — bind to terrain grid
- `pf.bind(ground, nullptr)` — cost from ground height
- `pf.allocate()` — allocate search buffers

### Config-Driven Init

`loadConfig(fname)` parses `world.ini`:
- `techs <file>` — load technology table
- `terrain <size> <step>` — allocate + generate/load terrain
- `rivers <min_flow>` — find rivers
- `straight_road <x0> <y0> <x1> <y1>` — build road
- `vehicles <n> [type]` — spawn vehicles
- `load_terrain` / `save_terrain` — cache I/O

## Parity Status

- **C++ terrain ↔ JS terrain**: Different algorithms (bisecNoise+erosion vs FBM+WGSL). No parity test. C++ has erosion, JS does not.
- **C++ water relaxation ↔ JS**: JS has no water simulation. C++ only.
- **C++ economy ↔ JS**: JS has no economy. C++ only.
- **C++ roads/vehicles ↔ JS**: JS has no roads/vehicles. C++ only.
- **C++ basin filling algorithms**: Bellman-Ford, Dijkstra, contour — tested via `test_landcraft.py` with comparison of `moveCost` output.

## Open Issues

- `LandCraftWorld::generateTerrain()` requires power-of-2 grid for bisecNoise — exits on non-power-of-2
- `HydraulicGrid2D::relaxWaterRasterY()` has `wsum*n.x/(i-istart)` — may be incorrect for non-square grids
- `Economy::Factory::produce()` doesn't check storage capacity
- `RoadBuilder::pushStright()` — typo in name (should be `pushStraight`)
- `LandCraftApp` has many commented-out members being migrated to `LandCraftWorld` — incomplete migration
- JS LandCraft has no hydraulics, economy, roads, or rivers — feature gap
- `PathFinder` not fully audited — algorithm unknown
- `HydraulicGrid2D::overlap_Line()` and `overlap_Triangle()` not implemented
- `loadConfig()` exits on file not found — should return false
- `droplerErosion()` uses `rand()` — not deterministic without `srand()`
- `LandCraftWorld::makeRivers()` sorts rivers by length but doesn't update `mouth` pointers
