# 2D Terrain & Tactics – JavaScript Plan (WIP)

This document summarizes the **2D top‑view, map‑centric projects** that we may want to reimplement in JavaScript/WebGL, and lists the relevant C++ files/functions.

Focus:

- Modern / medieval **armies and sieges** (LandTactics, FormationTactics, CastleBuilder).
- **Terrain, hydraulics, roads, and economy** (LandCraft).
- Shared sub‑problems around 2D maps: terrain fields, hydraulics/erosion, linear structures, visibility, movement, collisions, and danger/cover heat‑maps.

---

## 1. Projects and Grid Types

### 1.1 LandCraft – geographic economy on hex grid (+ square helpers)

- Docs:
  - `docs/LandCraft/LandCraft_main.md`
  - `docs/LandCraft/LandCraft.md`
  - `docs/LandCraft/UserGuide.md`
- Apps / code:
  - `cpp/apps/LandCraft/LandCraft_main.cpp` – main SDL/OpenGL app (`LandCraftApp`).
  - `cpp/apps/LandCraft/LandCraftWorld.cpp/.h`
  - `cpp/apps/LandCraft/Roads.h`
  - `cpp/common/maps/TerrainHydraulics.cpp/.h` (`HydraulicGrid2D`).
  - `cpp/common/maps/PathFinder.h`.
  - `cpp/common_SDL/SDL2OGL/TiledView.cpp/.h` – tiled rendering helper.
- **Grid:**
  - Primary: **hexagonal grid** via `SimplexRuler` (`ruler`).
  - Secondary: **square grid helper** via `Ruler2DFast` (`square_ruler`) for some calculations.

### 1.2 LandTactics – modern/WWII tactics on shared terrain

- Docs:
  - `docs/LandTactics/LandTactics.md`
- Apps / code:
  - `cpp/apps/LandTactics/LandTactics_main.cpp`
  - `cpp/apps/LandTactics/LTWorld.cpp/.h`
  - `cpp/apps/LandTactics/LTSquad.cpp/.h`
  - `cpp/apps/LandTactics/LTUnit.cpp/.h`
  - `cpp/apps/LandTactics/LTSurroundings.cpp/.h`
  - `cpp/apps/LandTactics/LTShelter.h`
  - `cpp/apps/LandTactics/LTdraw.h`, `LTrender.h`
- **Grid:**
  - Shares terrain infra with LandCraft:
    - `SimplexRuler` (hex grid)
    - `Ruler2DFast` (square grid)
  - World terrain is backed by `HydraulicGrid2D` (same as LandCraft).

### 1.3 FormationTactics – formation‑level melee / line battles

- Docs:
  - `docs/FormationTactics` (folder – details not read yet, but conceptually related).
- Apps / code:
  - `cpp/apps/FormationTactics/FormationTactics_main.cpp`
  - `cpp/apps/FormationTactics/FormationWorld.cpp/.h`
  - `cpp/apps/FormationTactics/Formation.cpp/.h`
  - `cpp/apps/FormationTactics/Faction.cpp/.h`
  - `cpp/apps/FormationTactics/Soldier.cpp/.h`, `SoldierType.cpp/.h`
  - `cpp/apps/FormationTactics/BattleLine.cpp/.h`, `PolyLineFormation.h`
- **Grid:**
  - Uses continuous 2D positions (`Vec2d`) for soldiers and formations. Terrain may still be sampled from a 2D field (height/cost), but the core interaction is **particle‑like** in continuous space.

### 1.4 CastleBuilder – static structures on terrain

- Apps / code:
  - `cpp/apps/CastleBuilder/CastleBuilder_main.cpp`
  - `cpp/apps/CastleBuilder/CastleWorld.cpp/.h`
- **Grid:**
  - Also uses a 2D map; likely reuses similar height‑field ideas for castle placement and walls (needs a closer read later).

---

## 2. Shared Sub‑Problems for 2D JS/WebGL Reimplementation

Across these projects we see the same core building blocks. For JS we should treat them as **reusable modules**.

### 2.1 Terrain representation and rendering

- **Height fields / scalar maps**
  - LandCraft: `HydraulicGrid2D` with `ground[]` (terrain height) and `water[]` (water level).  
  - LandTactics: `LTWorld` holds `HydraulicGrid2D hydraulics; double* ground;` and uses `SimplexRuler`/`Ruler2DFast` for indexing.

- **Tiled rendering**
  - `cpp/common_SDL/SDL2OGL/TiledView.cpp/.h`:
    - `TiledView::renderAll`, `TiledView::shiftRender`, `TiledView::draw_raw`
    - Maintains a grid of tiles and only redraws those entering the viewport.

- **JS/WebGL direction:**
  - Represent terrain as **2D textures / height maps**.
  - Use fragment shaders to color by height, water, cost, etc.
  - Implement a **tile/LOD system** similar to `TiledView` for large worlds.

### 2.2 Hydraulics, erosion, flooding

- **Core class:** `HydraulicGrid2D` in `TerrainHydraulics.cpp/.h`.
  - Droplet erosion: `initDroplet`, droplet stepping (erosion / sediment).
  - Basin fill / flooding algorithms:
    - `basinFill_Contour_*` (contour-based flood).
    - `basinFill_BellmanFord` (+ `_step`).
    - `basinFill_Dijkstra` – priority‑queue flood.
  - Water relaxation and rainfall / river extraction (see LandCraft docs & app).

- **LandCraft integration:**
  - `LandCraftApp` functions:
    - `generateTerrain()` – noise + erosion.
    - `hydroRelaxUpdate()`, `hydro1D_update()`.
    - `makeRivers()`, `traceDroplet`, inflow/outflow commands.

- **JS/WebGL direction:**
  - Start with **CPU height‑field kernels** that mimic a subset of `HydraulicGrid2D` (e.g. simple relaxation, one basin‑fill variant).
  - Later port to compute shaders / fragment passes for:
    - Water relaxation on a texture.
    - Droplet erosion (random walks on texture grid).

### 2.3 Linear structures: rivers, roads, walls

- **Roads / transport:** `cpp/apps/LandCraft/Roads.h`
  - Data structures:
    - `RoadTile { uint16_t ia, ib; double height; }` – grid indices + elevation along road.
    - `Road` – array of `RoadTile`s.
    - `RoadBuilder` – interactive construction (`pushStright` for hex grid paths).
    - `RoadVehicle`, `RoadVehicleType` – vehicles moving along `road->path`, speed vs slope.

- **Pathfinding:** `cpp/common/maps/PathFinder.h`
  - `PathFinder : Grid2DAlg` with:
    - `height[]`, `terrain_cost[]`.
    - `moveCosts[]`, `toBasin[]`, `toTile[]`.
  - Algorithms for multi‑center path networks (ways between basins/centers), building `Way` objects with `path` indices.

- **Walls / linear obstacles:**
  - In LandTactics / CastleBuilder, linear objects (walls, ditches, hedges) are represented as `LTLinearObject` or similar (line segments with width and cover properties).

- **JS/WebGL direction:**
  - Hex or square grid **edge graph** for rivers, roads, walls.
  - Simple **road editor** in JS (click‑and‑drag to create polylines on grid; export as `RoadTile`‑like arrays).
  - Pathfinding based on `height` + terrain penalties → used both for AI movement and procedural road routing.

### 2.4 Visibility and line-of-sight on terrain

- **LandTactics visibility:**
  - Docs: `docs/LandTactics/LandTactics.md` – visibility/LOS as core tactical mechanic.
  - Code: `LTWorld`, `LTrender`, `LTdraw` probably include LOS checks and debug draws (exact functions not yet enumerated).

- **JS/WebGL direction:**
  - LOS over height field via ray marching in height map:
    - For each candidate direction from a unit, sample terrain heights along ray, check if any point rises above eye line.
  - Visualizations:
    - "Visible" overlay for a unit by painting a coverage texture.
    - Combined visibility of a whole faction.

### 2.5 Movement and collisions of units on terrain

- **LandCraft:**
  - Vehicles moving along roads (`RoadVehicle::move`), slope‑dependent speeds (TODO in `getSpeed`).

- **LandTactics:**
  - `LTSquad`, `LTUnit` – squads and units moving in continuous 2D space.
  - `LTWorld::simulationStep(double dt)` calls `f->update(dt)` for each `LTSquad`.
  - Terrain and hydraulics (`HydraulicGrid2D`) are available to modulate movement costs (e.g. uphill/downhill, soft ground) – details to be confirmed inside squad/unit update logic.

- **FormationTactics:**
  - `Formation.cpp` – continuous 2D positions for `Soldier`s, with:
    - Bounding boxes: `Formation::update_bbox()`.
    - Interactions:
      - `Formation::interact(Formation* fb)` – pairwise soldier interactions between formations, separated into friend/enemy cases.
      - `Formation::interactInside()` – cohesion/spacing within one formation.
    - Targeting and reorientation: `Formation::setTarget(...)`, center/rotation updates.
  - `Salvo` logic for volley fire (hit distribution along a line).

- **JS/WebGL direction:**
  - Simple **2D particle dynamics** for units with:
    - Terrain‑dependent move cost (sample slope / softness / road presence).
    - Collision avoidance / cohesion (boids‑like forces or FormationTactics‑style all‑pairs interaction within formations).
  - Optionally offload some per‑unit updates to GPU (compute or vertex shader) for large crowds.

### 2.6 Danger / cover heat-maps (LandTactics)

This is the most interesting for GPU visualization.

- **Key types & functions:**
  - `cpp/apps/LandTactics/LTSurroundings.h`:
    - `class LTsurrounding` holds:
      - `std::vector<LTStaticObject*> objs;`
      - `std::vector<LTLinearObject*> lobjs;` – linear features like walls, hedges, obstacles, each with `width` and `cover`.
      - `std::vector<LTUnit*> enemies;`
      - `std::vector<LTUnit*> coleagues;`
    - Parameters:
      - Search radii: `RSearchLocal`, `RSearchGlobal`.
      - Deployment distance: `Rdeploy`.
      - Weights: `cover2E`, `overlap2E`.
      - Optional constraint region: `bConstr`, `ConstrPos`, `ConstrRad`, `ConstrE`.
    - `double unitPosFittness(const LTUnit* u, const Vec2d& p) const`:
      - Computes **fitness / utility** of putting unit `u` at position `p`:
        - For each `LTLinearObject* lo` in `lobjs`:
          - Compute squared distance `r2 = lo->dp(p).norm2()` to segment.
          - If `r2 < width^2`, increase `cover` up to `lo->cover * (1 - r2/r2max)`.
        - For each colleague `uu` in `coleagues` (except `u`):
          - Add `overlap += 1 / (1 + (dist^2 / Rdeploy^2))` using `goal_pos`.
        - Combine: `E = overlap2E * overlap + cover2E * cover`.
        - Optional constraint penalty around `ConstrPos`.
        - TODO: add explicit **danger from enemies**.
    - `void tryFindBetterPlace(LTUnit* u, int n, bool bJumpToGoal) const`:
      - Stochastic search around `u->goal_pos`:
        - With probability 0.8: local search in radius `RSearchLocal`.
        - Otherwise: global search in radius `RSearchGlobal`.
      - Accepts only moves that **increase E** (greedy hill climb) and updates `u->goal_pos`.

  - `cpp/apps/LandTactics/LTWorld.cpp`:
    - `void LTWorld::getSurroundings(LTsurrounding& sur, const LTFaction* fac, const Vec2d& pos, double R)` – populates `sur` for a given position.
    - `void LTWorld::prepareSurroundings(...)` – sets up search radii & constraint region for deployment.
    - `void LTWorld::optimizeDeployment(LTSquad* s, double R, int n, int m, bool bJumpToGoal)` – high‑level deployment optimization over all units in a squad using `tryFindBetterPlace` repeatedly.

- **JS/WebGL direction:**
  - CPU prototype:
    - Implement a JS `Surroundings` struct with wall/cover segments and units.
    - Evaluate `unitPosFitness` on a **grid of sample points** → produce a scalar field (danger/cover heat‑map).
  - GPU version (GLSL):
    - Represent walls/cover segments and squad positions as uniforms/buffers.
    - In a fragment shader over the terrain texture, compute:
      - Distance to nearest cover segment; convert to local cover value.
      - Penalty for being too close to friendly units.
      - Potential future term for enemy LOS / weapon ranges.
    - Output **heat‑map color** directly.
  - This can then be used both for **visualization** and for driving a JS solver that picks better deployment positions.

---

## 3. JS/WebGL Reimplementation: High-Level Plan

For a future JS/TS module (most likely under a `land2d` sub‑namespace or a dedicated editor mode), we can structure features roughly as:

- **Terrain module**
  - Height field + scalar layers (water, cost, danger).
  - Texture‑based rendering, TiledView‑like LOD.

- **Hydraulics module**
  - Simple relaxation + basin fill in JS.
  - Optional GPU implementation for water flooding & erosion.

- **Linear structures module**
  - Roads, rivers, walls as polylines over hex/square grid.
  - Road/path editing + simple pathfinding (PathFinder‑inspired).

- **Units & formations module**
  - Particle‑like units with terrain‑aware movement.
  - Formation behaviors (from FormationTactics) for melee/line battles.

- **Visibility & danger module**
  - LOS over height maps.
  - Danger/cover heat‑maps from LandTactics `LTsurrounding` logic.

All of these share the same **2D map representation**, so careful design of the map/grid/core data structures in JS will pay off across LandCraft‑style economy sims, LandTactics‑style combat, and FormationTactics‑style battles.
