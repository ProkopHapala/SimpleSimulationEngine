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

---

## Algorithm Review: Path-finding, Flooding, Rivers, Routes

This section reviews algorithms relevant to terrain processing (flooding, river finding, path-finding between cities, road optimization). For each algorithm we describe: merit, use cases, strengths, weaknesses, GPU parallelization potential, and implementation status across platforms (C++, Python, OpenCL, JavaScript).

### A. Flooding / Basin Filling (Lakes)

#### A1. Bellman-Ford Wave Relaxation

- **Merit**: Simple iterative relaxation of spill levels. Each cell updates `spill = min(current, max(ground, neighbor_spill))` over all neighbors, repeated until convergence.
- **Use cases**: Basin filling / depression filling on terrain grids. Finding lake water levels.
- **Strengths**: Trivial to implement. No priority queue needed. Works with any neighbor topology (4/8/6-hex). Cache-friendly raster scan.
- **Weaknesses**: O(N · iters) — convergence can be slow for large or convoluted basins (information propagates one cell per sweep per direction). May need many iterations on large maps.
- **GPU parallelization**: Good — maps to Jacobi-style kernel where each thread updates one cell from neighbors, with global barrier between iterations. Convergence is the bottleneck; can be accelerated with jump-flooding or hierarchical approaches. O(log N) iterations with JFA, O(N) without.
- **Implemented**:
  - C++: `HydraulicGrid2D::basinFill_BellmanFord()`, `basinFill_BellmanFord_step()` in `cpp/common/maps/TerrainHydraulics.cpp:360-403`. Tested via `test_landcraft.py`.
  - OpenCL: `solve_tiles` kernel in `python/terrain_ocl/terrain.cl:12-146` uses Bellman-Ford relaxation in 16×16 shared-memory tiles (64 iterations).
  - Python: `solve_relaxation_numpy()` in `python/terrain_ocl/terrain.py` — vectorized CPU reference.
  - JS: Not implemented.

#### A2. Dijkstra Priority-Flood

- **Merit**: Optimal single-pass basin filling using a priority queue. Pop lowest spill cell, propagate `max(ground, spill)` to neighbors. Each cell visited exactly once.
- **Use cases**: Basin filling where correctness and speed matter. Reference algorithm for CPU terrain processing.
- **Strengths**: O(N log N) — optimal for CPU. Visits each cell once. Mathematically correct spill levels.
- **Weaknesses**: Priority queue (heap) operations are inherently sequential and random-access. Poor cache behavior on large grids. Does not parallelize well.
- **GPU parallelization**: Poor — priority queues are globally synchronized, causing massive thread contention on GPU. This is the classic "GPU-unfriendly" algorithm. Alternatives: iterative relaxation (A1), hierarchical tile-graph (A4), or Borůvka MST (A5).
- **Implemented**:
  - C++: `HydraulicGrid2D::basinFill_Dijkstra()` in `cpp/common/maps/TerrainHydraulics.cpp:406-440`. Uses `std::priority_queue`. Tested.
  - Python: `global_dijkstra()` in `python/terrain_ocl/terrain.py:276` — CPU reference for verification.
  - JS: Not implemented.

#### A3. Contour Two-Wave (Flood + Drain)

- **Merit**: Novel hybrid: Phase 1 floods from seeds via contour expansion (like BFS wavefront) up to a level cap. Phase 2 drains back from boundary, monotonically lowering spill levels to correct values. No priority queue needed.
- **Use cases**: Basin filling with progressive/steppable visualization. Allows incremental rendering of flood propagation.
- **Strengths**: No PQ needed. Steppable (can pause/resume). Contour expansion is cache-friendly. Good for interactive visualization.
- **Weaknesses**: Two-phase complexity. May require many steps for large basins. The flood phase may over-estimate, requiring drain phase to correct.
- **GPU parallelization**: Moderate — contour expansion can be parallelized (each frontier cell processed independently), but the frontier size varies per step. Drain phase is similar. Better than PQ-Dijkstra but worse than JFA or hierarchical.
- **Implemented**:
  - C++: `HydraulicGrid2D::basinFill_Contour_beginFlood/floodStep/beginDrain/drainStep/finish()` in `cpp/common/maps/TerrainHydraulics.h:406-416`, implementation in `TerrainHydraulics.cpp`. Exposed via C API `lc_basin_contour_*` in `cpp/libs/CombatModels/LandCraftLib.cpp:175-179`.
  - Python/OpenCL/JS: Not implemented.

#### A4. Hierarchical Tile-Graph Flooding

- **Merit**: Divide map into 16×16 (or 32×32) tiles. Solve each tile locally (intra-tile relaxation/Dijkstra in shared memory). Build a graph where nodes = tile sinks, edges = lowest-cost boundary crossings (gates). Solve flooding on the coarse graph. Broadcast results back to pixels.
- **Use cases**: Large/megapixel terrain maps where full-grid algorithms are too slow. GPU-friendly.
- **Strengths**: O(N) with O(log N) hierarchical levels. GPU-friendly: intra-tile solve in shared memory, inter-tile on coarse graph. Scales to 4k×4k+ maps.
- **Weaknesses**: Tile boundary pathology — narrow passes on tile boundaries may be missed. Requires careful gate computation. Implementation complexity.
- **GPU parallelization**: Excellent — this is the primary GPU strategy. Intra-tile: shared memory relaxation (like `solve_tiles` kernel). Inter-tile: small graph on CPU or parallel Borůvka. Broadcast: simple parallel kernel.
- **Implemented**:
  - OpenCL: `solve_tiles` kernel in `python/terrain_ocl/terrain.cl:12-146`. Gate finding in `terrain.py:compute_gates_simple()`, `compute_gates_diag()`. Tile graph Dijkstra in `terrain.py:tile_graph_shortest()`. Tested via `test_terrain.py`, `verify_real.py`.
  - C++: Not implemented (could port from OpenCL).
  - JS: Not implemented (design discussed in `js/LandCraft_web/doc/Flooading.md`).

#### A5. Watershed via Steepest Descent + Pointer Jumping

- **Merit**: Each pixel points to its steepest-descent neighbor (flow map). Then pointer jumping propagates basin IDs: `basin[i] = basin[flow[i]]`, repeated until convergence. All pixels end up pointing to their sink.
- **Use cases**: Watershed segmentation, basin identification, finding drainage basins. Does not compute water levels — only basin membership.
- **Strengths**: O(N · log N) iterations on GPU. Embarrassingly parallel. Simple kernels. No priority queue.
- **Weaknesses**: Only finds basin membership, not spill levels or water depths. Requires separate pass for depression filling. Steepest descent can create artifacts (single-pixel channels).
- **GPU parallelization**: Excellent — two simple parallel kernels (`calc_flow`, `propagate_basins`), each thread processes one pixel. Pointer jumping converges in O(log N) iterations.
- **Implemented**:
  - OpenCL: `calc_flow` kernel in `python/terrain_ocl/terrain.cl:148-225`. `propagate_basins` kernel in `terrain.cl:230-270`. Python wrapper in `terrain.py:run_calc_flow()`, `run_propagate_basins()`. Tested via `test_watershed.py`.
  - C++: Not implemented (C++ uses `gatherRain()` instead — see B1).
  - JS: Not implemented.

#### A6. Borůvka Minimum Spanning Tree (Design Only)

- **Merit**: Parallel MST algorithm. Every node simultaneously picks cheapest edge, contracts connected components, repeats. O(log V) iterations.
- **Use cases**: Finding spill connectivity between basins on large maps. Replaces iterative relaxation with fast parallel contraction.
- **Strengths**: Most parallel MST algorithm. O(log V) iterations. Perfect for GPU.
- **Weaknesses**: Implementation complexity. Requires graph construction (tile-level or basin-level).
- **GPU parallelization**: Excellent — designed for parallel execution.
- **Implemented**: Not implemented. Discussed in `js/LandCraft_web/doc/Flooading.md` as a design option for GPU flooding at megapixel scale.

#### A7. Jump Flooding Algorithm (JFA) (Design Only)

- **Merit**: Propagate information across grid in O(log N) steps using power-of-two step sizes. Step 1: check neighbors at distance 2^(n-1), step 2: 2^(n-2), ... until distance 1.
- **Use cases**: Fast approximate flooding on GPU. Can bridge large distances quickly.
- **Weaknesses**: May skip thin barriers (tunneling through walls). Needs conservative rasterization or multi-sampling to avoid.
- **GPU parallelization**: Excellent — designed for GPU. Each step is a simple parallel kernel.
- **Implemented**: Not implemented. Discussed in `js/LandCraft_web/doc/Flooading.md`.

---

### B. River Finding

#### B1. Rain Accumulation + Upstream Tracking

- **Merit**: Sort all cells by height (descending). Each cell passes its accumulated water to the lowest neighbor. Cells with no lower neighbor become sinks. Then trace rivers upstream from sinks, following highest-flow neighbors.
- **Use cases**: River network extraction from terrain. Produces river trees with flow values.
- **Strengths**: Simple. Produces flow accumulation map + river paths + river tree. Handles branching via recursive tracking.
- **Weaknesses**: O(N log N) due to sort. Sequential processing (each cell depends on higher cells being processed first). Does not parallelize easily. Steepest-descent only — no multi-direction flow (D∞ or D8 variants).
- **GPU parallelization**: Poor for the accumulation step (sequential dependency chain). The steepest-descent flow map can be computed in parallel (like `calc_flow` kernel), but accumulation requires following chains. Possible approaches: (1) pointer jumping to find sinks, then atomic-add flow along paths; (2) hierarchical reduction. River tracing is inherently sequential but could be done per-river in parallel.
- **Implemented**:
  - C++: `HydraulicGrid2D::gatherRain()` in `cpp/common/maps/TerrainHydraulics.cpp:152-181`. `traceDroplet()` at `:183-207`. `trackRiver()` at `:209-251`. `trackRiverRecursive()` at `:253-277`. `findAllRivers()` at `:279-289`. Exposed via C API `lc_gather_rain()`, `lc_find_all_rivers()`, `lc_trace_droplet()`, `lc_river_*` in `LandCraftLib.cpp:94-101`.
  - Python/OpenCL/JS: Not implemented. The OpenCL watershed (`calc_flow`) computes the flow map but does not accumulate flow or trace rivers.

#### B2. GPU Flow Accumulation (Not Implemented)

- **Merit**: Use the GPU-computed steepest-descent flow map (`calc_flow` kernel) and accumulate flow in parallel. Each pixel starts with 1 unit of water. Use pointer jumping or parallel reduction along flow paths to accumulate.
- **Use cases**: Fast river finding on large maps.
- **GPU parallelization**: Moderate — flow accumulation along paths is inherently sequential, but can be parallelized per-path or using hierarchical reduction. Atomic adds on sink cells. This is an open research area.
- **Implemented**: Not implemented. The flow map exists (`calc_flow` kernel) but accumulation does not.

---

### C. Path-finding / Optimal Routes Between Cities

#### C1. Multi-Source Contour Expansion (Bellman-Ford-like)

- **Merit**: All city centers seeded simultaneously with cost=0. Each step expands the frontier: for each frontier cell, check neighbors, compute `cost = moveCost + neigh_dist + heightCost(dh) + terrain_cost`. If better, add to new frontier. Produces a Voronoi-like partition (each cell assigned to nearest city) with optimal costs.
- **Use cases**: Finding optimal routes between all pairs of cities simultaneously. Territory partitioning. Network connectivity analysis.
- **Strengths**: Multi-source — solves all-pairs in one pass. Simple to implement. No priority queue. Cost function is flexible (height, terrain, distance). After expansion, `findConnections()` finds optimal border crossings between city regions, and `track()` reconstructs full paths.
- **Weaknesses**: O(N · diameter) — contour expansion is BFS-like, cells may be revisited many times. Slower than PQ-Dijkstra for single-pair queries. No heuristic (unlike A*). No early termination.
- **GPU parallelization**: Good — each frontier cell can be processed independently. Frontier size varies but is embarrassingly parallel within each step. Similar to basin-filling contour expansion. Could use JFA for acceleration.
- **Implemented**:
  - C++: `PathFinder` class in `cpp/common/maps/PathFinder.h:14-213`. `pepare()` at `:72-88`, `path_step()` at `:93-113`, `extend_path()` at `:115-126`, `findConnections()` at `:128-163`, `track()` at `:165-187`, `makePaths()` at `:189-195`. Used in `cpp/apps/LandTactics/LTWorld.cpp:300-309`. Bound to terrain in `LandCraftWorld::BindPathFinderToMap()` at `LandCraftWorld.cpp:433-439`. Exposed via C API `lc_pf_*` in `LandCraftLib.cpp:143-161`.
  - Python/OpenCL/JS: Not implemented.

#### C2. Priority-Queue Dijkstra (Single-Source)

- **Merit**: Standard Dijkstra with min-heap. Optimal O(N log N) for single-source shortest path.
- **Use cases**: Single-pair route queries (e.g., "road from city A to city B"). Faster than multi-source for one-off queries.
- **Strengths**: Optimal CPU algorithm. Each cell visited once. Well-understood.
- **Weaknesses**: Priority queue is sequential, poor cache behavior. Does not parallelize well. Only finds one path per run.
- **GPU parallelization**: Poor — same issues as A2. PQ is GPU-unfriendly.
- **Implemented**:
  - Python: `global_dijkstra()` in `python/terrain_ocl/terrain.py:276` — CPU reference for verification.
  - C++: Not implemented in `PathFinder` (uses contour expansion instead). The `basinFill_Dijkstra` code in `TerrainHydraulics.cpp` demonstrates the PQ pattern but is for basin filling, not path-finding.
  - JS: Not implemented.

#### C3. A* Path-finding (Not Implemented)

- **Merit**: Dijkstra + heuristic (e.g., Euclidean distance to goal). Dramatically faster for single-pair queries — explores only cells likely on the optimal path.
- **Use cases**: Real-time road building between two cities. Interactive route planning.
- **Strengths**: Much faster than Dijkstra for single-pair. Optimal if heuristic is admissible. Well-known, widely used.
- **Weaknesses**: Single-pair only (not multi-source). Heuristic must be admissible for optimality. Memory for open/closed sets.
- **GPU parallelization**: Poor — same PQ issue as Dijkstra. However, parallel A* variants exist (e.g., PA*, HPA*). For GPU, hierarchical approaches (C4) are better.
- **Implemented**: Not implemented anywhere. Natural next step for single-pair road routing in C++.

#### C4. Hierarchical Path-finding (HPA* / Tile-Graph)

- **Merit**: Divide map into tiles. Solve path within each tile (local Dijkstra/relaxation in shared memory). Build coarse graph between tile boundaries (gates). Find path on coarse graph, then refine within tiles.
- **Use cases**: Large maps where full-grid path-finding is too slow. GPU path-finding.
- **Strengths**: O(N) with hierarchical levels. GPU-friendly (intra-tile in shared memory). Scales to megapixel maps.
- **Weaknesses**: Tile boundary pathology. May not find true optimal path (approximation). Implementation complexity.
- **GPU parallelization**: Excellent — intra-tile solve in shared memory (like `solve_tiles` kernel), inter-tile on small graph. This is the standard GPU path-finding approach.
- **Implemented**:
  - OpenCL: `solve_tiles` kernel in `python/terrain_ocl/terrain.cl:12-146` solves intra-tile paths. `tile_graph_shortest()` in `terrain.py:163` finds inter-tile path. `get_pixel_path()` in `terrain.py:260` reconstructs full path. Tested via `test_terrain.py`, `verify_real.py`.
  - C++: Not implemented (could port from OpenCL).
  - JS: Not implemented.

---

### D. Road Building

#### D1. Straight-Line Bresenham (Grid Line Drawing)

- **Merit**: Draws a line of grid cells from A to B using a Bresenham-like algorithm adapted for hex/square grids.
- **Use cases**: Quick straight roads. Simple path generation without terrain awareness.
- **Strengths**: O(length) — very fast. Simple.
- **Weaknesses**: Not terrain-aware — goes straight regardless of slopes, water, or obstacles. No cost optimization.
- **GPU parallelization**: N/A — too simple to warrant GPU.
- **Implemented**:
  - C++: `RoadBuilder::pushStright()` in `cpp/apps/LandCraft/Roads.h:32-137`. Used by `LandCraftWorld::makeRoadStraight()`. Exposed via C API `lc_road_*` in `LandCraftLib.cpp:108-113`.
  - Python/OpenCL/JS: Not implemented.

#### D2. Terrain-Aware Road Optimization (Not Implemented)

- **Merit**: Use path-finding (A* or multi-source Dijkstra) with a cost function (slope, terrain difficulty, water crossing cost) to find the optimal road route between two points.
- **Use cases**: Realistic road building that avoids steep terrain, follows valleys, minimizes construction cost.
- **Strengths**: Produces realistic, terrain-aware roads. Can use existing `PathFinder` cost function (`heightCost` with `ch2`, `chplus`, `chminus`).
- **Weaknesses**: Requires path-finding infrastructure (A* or Dijkstra). More expensive than straight line.
- **GPU parallelization**: Use hierarchical approach (C4) for large maps. Single road on CPU is fine.
- **Implemented**: Not implemented. The `PathFinder` infrastructure exists (cost function, neighbor expansion) but is not connected to `RoadBuilder`. This is a natural integration point: replace `pushStright()` with a call to `PathFinder` for terrain-aware routing.

---

### E. Water Relaxation (Local)

#### E1. Scanline Water Flattening

- **Merit**: Scan along rows (X) or columns (Y), find contiguous water blocks, flatten to average water level.
- **Use cases**: Fast approximate water relaxation. Pre-processing for basin filling.
- **Strengths**: Very fast, cache-friendly. O(N) per scan.
- **Weaknesses**: Only 1D relaxation — does not handle 2D flow correctly. Known bug in `relaxWaterRasterY()` for non-square grids (line 136: `wsum*n.x/(i-istart)` — should be `wsum/(count)` where count accounts for stride).
- **GPU parallelization**: Good — each row/column is independent, can be processed in parallel.
- **Implemented**:
  - C++: `HydraulicGrid2D::relaxWaterRasterX()` in `TerrainHydraulics.h:75-111`, `relaxWaterRasterY()` at `:113-149`.
  - Python/OpenCL/JS: Not implemented.

#### E2. Local Neighbor Water Leveling

- **Merit**: For each cell, load neighbors, sort by height, find common water level by distributing excess water from lowest to highest. Sets water level for cell + all neighbors.
- **Use cases**: Physically-based local water relaxation. More accurate than scanline.
- **Strengths**: Handles 2D flow. Conserves water volume. Works with any neighbor topology.
- **Weaknesses**: O(N · nneigh · log(nneigh)) per sweep. Needs multiple sweeps for convergence. Sequential (each cell modifies neighbors — race conditions if parallelized naively).
- **GPU parallelization**: Poor — each cell writes to neighbors, causing race conditions. Would need double-buffering or red-black ordering. Jacobi-style variant possible but convergence slower.
- **Implemented**:
  - C++: `HydraulicGrid2D::relaxWater()` in `TerrainHydraulics.h:201-327`. `relaxWater2cells()` at `:151-175`. `relaxWaterHexCells()` at `:177-199`. Full grid sweep `relaxWater()` at `:329-335`.
  - Python/OpenCL/JS: Not implemented.

---

### F. Implementation Gap Summary

| Algorithm | C++ | Python | OpenCL | JS | Gap |
|-----------|-----|--------|--------|-----|-----|
| Bellman-Ford basin fill | ✅ | ✅ (ref) | ✅ (tile) | ❌ | JS needs implementation |
| Dijkstra basin fill | ✅ | ✅ (ref) | ❌ | ❌ | OpenCL/JS gap (GPU needs alternatives) |
| Contour two-wave | ✅ | ❌ | ❌ | ❌ | Only C++ |
| Hierarchical tile flooding | ❌ | ✅ | ✅ | ❌ | C++ could port from OpenCL |
| Watershed (pointer jumping) | ❌ | ✅ | ✅ | ❌ | C++ could use GPU approach |
| Borůvka MST | ❌ | ❌ | ❌ | ❌ | Design only (`Flooading.md`) |
| Jump Flooding (JFA) | ❌ | ❌ | ❌ | ❌ | Design only (`Flooading.md`) |
| Rain accumulation + rivers | ✅ | ❌ | ❌ | ❌ | GPU river finding missing |
| GPU flow accumulation | ❌ | ❌ | ❌ | ❌ | Open research area |
| Multi-source contour pathfind | ✅ | ❌ | ❌ | ❌ | Only C++ |
| PQ Dijkstra pathfind | ❌ | ✅ (ref) | ❌ | ❌ | C++ could add single-source mode |
| A* pathfind | ❌ | ❌ | ❌ | ❌ | Natural next step |
| Hierarchical pathfind (HPA*) | ❌ | ✅ | ✅ | ❌ | C++ could port |
| Straight road (Bresenham) | ✅ | ❌ | ❌ | ❌ | Only C++ |
| Terrain-aware road | ❌ | ❌ | ❌ | ❌ | PathFinder + RoadBuilder integration |
| Scanline water relaxation | ✅ | ❌ | ❌ | ❌ | Only C++ (has bug) |
| Local water leveling | ✅ | ❌ | ❌ | ❌ | Only C++ |
