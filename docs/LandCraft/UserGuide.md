# LandCraft User Guide

## Overview
LandCraft is a geographic economy simulation where terrain, water, roads, and vehicles interact. You build infrastructure that respects geographic constraints (slope, rivers, basins) and observe how water flow and transport efficiency emerge from the underlying systems.

This guide focuses on features and how to test them during play, not the code internals.

### Where in code
- Core app class: `LandCraftApp` in `cpp/apps/LandCraft/LandCraft_main.cpp`
- 2D hydraulics: `HydraulicGrid2D` in `cpp/common/maps/TerrainHydraulics.h`
- Roads/vehicles: `RoadBuilder`, `Road`, `RoadVehicle` in `cpp/apps/LandCraft/Roads.h`
- Pathfinding (planned): `PathFinder` in `cpp/common/maps/PathFinder.h`


## World and Water

- Terrain generation
  - Press `N` to generate a new procedural terrain. Elevations are synthesized and eroded to produce realistic features.
  - Code: `LandCraftApp::generateTerrain()`; invoked via `LandCraftApp::commandDispatch(Actions::generateTerrain)` from `LandCraftApp::eventHandling()`.

- Terrain editing
  - Left-click picks the current paint height from the tile under the cursor.
  - Left-drag paints ground elevation to the picked height on hovered tiles.
  - Code: pick with `LandCraftApp::eventHandling()` on `SDL_MOUSEBUTTONDOWN` (left); paint in `LandCraftApp::mouseHandling()` modifying `hydraulics.ground[ihex]`.

- Water relaxation (lakes/levels)
  - Press `R` to relax water everywhere; water locally equalizes subject to surrounding heights.
  - Press `E` to relax water only at the hex under the cursor (fine control).
  - Code: global `HydraulicGrid2D::relaxWater()` and local `HydraulicGrid2D::relaxWater(const Vec2i&)` called via `LandCraftApp::commandDispatch(Actions::relaxWater / Actions::relaxWaterHex)`.

- Inflow/Outflow (dams/canals experiments)
  - Press `O` to set an outflow boundary at the cursor tile (drains to ground level).
  - Press `I` to set an inflow source at the cursor tile (raises local water).
  - Combine with relaxation to observe drainage and flooding.
  - Code: handled in `LandCraftApp::commandDispatch(Actions::outflow / Actions::inflow)` using `HydraulicGrid2D` outflow buffers (`allocate_outflow()`, `contour2`, `isOutflow`).

- Rivers
  - Press `G` to simulate uniform rainfall and extract river networks.
  - A river dropdown appears; select a river to view its elevation and flow profile.
  - Code: `LandCraftApp::makeRivers()` calls `HydraulicGrid2D::gatherRain()` and `HydraulicGrid2D::findAllRivers()`; rivers are `River` structs with `path` and `flow`.
  - River plots: `LandCraftApp::addRiverPlot()`; dropdown built in `makeRivers()`.

- Droplet tracing
  - Press `T` at a tile to trace a droplet downhill. Toggle drawing of the trace via GUI.
  - Code: `LandCraftApp::commandDispatch(Actions::traceDroplet)` calls `HydraulicGrid2D::traceDroplet()`; drawn by `LandCraftApp::drawDropletTrace()` when enabled.


## Roads and Vehicles

- Build roads
  - Right-button press sets the start tile; right-button release sets the end tile.
  - A straight road on the hex grid is built between the two points.
  - Code: mouse right press/release in `LandCraftApp::eventHandling()`; road creation via `LandCraftApp::addRoadStright()` which uses `RoadBuilder::pushStright()` and `RoadBuilder::writeIt()`.

- Road profiles
  - The HUD shows road profiles: terrain height vs. water level along the road.
  - Toggle the road profile panel via GUI.
  - Code: computed in `LandCraftApp::addRoadPlot()` using `Plot2D`; rendered in `LandCraftApp::drawHUD()`.

- Vehicles
  - Vehicles can be spawned on roads (one demo vehicle is created by default).
  - Toggle running and drawing of vehicles in the GUI.
  - Code: created in `LandCraftApp::makeVehicles()`; movement via `RoadVehicle::move()` and `RoadVehicle::moveStep()` using `RoadVehicleType::getSpeed()`; toggles in `registerDrawLayers()`.


## GUI Toggles (Panels)

Enable or disable:
- Run 2D water relaxation
- Run 1D hydraulics demo
- Draw rivers
- Draw droplet trace
- Draw roads
- Run vehicles
- Draw vehicles
- Show road profile
- Show river profile

- Code: toggles are members on `LandCraftApp` (`bRunHydroRelax`, `bRunHydro1D`, `bDrawRivers`, `bDrawTraceDroplet`, `bDrawRoads`, `bRunVehicles`, `bDrawVehicles`, `bRoadProfile`, `bRiverProfile`) and are wired in `LandCraftApp::registerDrawLayers()` using `CheckBoxList` (file `cpp/apps/LandCraft/LandCraft_main.cpp`). The HUD draw calls occur in `LandCraftApp::drawHUD()`; the main draw uses flags in `LandCraftApp::draw()`.


## Controls Summary

- Keyboard
  - `N`: Generate terrain
  - `R`: Relax water everywhere
  - `E`: Relax water at cursor hex
  - `S`: Save terrain/water state to data/*.bin
  - `L`: Load terrain/water state (or generate if missing/corrupt)
  - `O`: Set outflow at cursor
  - `I`: Set inflow at cursor
  - `G`: Gather rain and find rivers
  - `M`: Toggle terrain visualization mode
  - `T`: Trace droplet from cursor

  - Code: bindings set in `LandCraftApp::registerCommands()`; dispatch handled in `LandCraftApp::eventHandling()` and routed to `LandCraftApp::commandDispatch()` switching over `LandCraftApp::Actions` (file `cpp/apps/LandCraft/LandCraft_main.cpp`).

- Mouse
  - Left-click: pick paint height from tile
  - Left-drag: paint terrain to picked height
  - Right-press: begin road
  - Right-release: finish road and build

  - Code: implemented in `LandCraftApp::eventHandling()` (press/release) and `LandCraftApp::mouseHandling()` (drag painting). Road building calls `LandCraftApp::addRoadStright()`.


## Visualization Modes

- Mode 1: Terrain elevation colored by height and water depth shading.
- Mode 2: River/flow visualization with grayscale water and red markers for discovered rivers.

- Code: selected via `LandCraftApp::commandDispatch(Actions::terrainViewMode)` toggling `terrainViewMode`; actual color selection in `LandCraftApp::terrainColor()` used by `LandCraftApp::drawTerrain()`.


## Suggested Test Scenarios

- Terrain editing and lakes
  - Pick a hilltop height and paint a basin; press `R` and observe lake formation.

- Rain to rivers
  - Press `G` and inspect the river dropdown; choose a river to plot its ground and flow.

- Roads vs. terrain
  - Build a road across varied terrain; view the road profile to see climbs and descents.

- Inflow/outflow control
  - Mark inflow and outflow points and relax water to simulate canal/dam behavior.

- Droplet sanity
  - Trace a droplet from a slope shoulder; verify the trace moves consistently downhill.


## Saving, Loading, and Command Scripts

- Save/Load
  - Use `S`/`L` to save/load `data/ground.bin` and `data/water.bin`.
  - On load failure or missing files, the game regenerates a new terrain.

  - Code: handled in `LandCraftApp::commandDispatch(Actions::save / Actions::load)` using `saveBin`/`loadBin` (file `cpp/apps/LandCraft/LandCraft_main.cpp`). Initial load/cached generation logic also in `LandCraftApp::makeMap()`.

- Command scripts
  - `data/comands.ini` can contain scripted commands executed at startup for repeatable setups.

  - Code: executed via `cmdPars.execFile("data/comands.ini")` in `LandCraftApp` constructor. Commands are defined by `Commander keybinds` and GUI `CommandList` (`registerCommands()`).


## Advanced: Pathfinding and Economy (Current State)

- Pathfinding
  - A pathfinding module exists to compute minimum-cost routes across basins, accounting for slope and terrain penalties. Integration for interactive path planning is planned.

- Economy
  - Technologies (recipes of inputs→outputs with cycle time and unit space) are loaded at startup. Factory simulation APIs exist but are not yet wired to the game loop.


## Tips

- Use GUI toggles to isolate subsystems for testing (e.g., relax water without vehicle updates).
- Keep an eye on water conservation qualitatively: large relax steps should not create or destroy water except at inflow/outflow boundaries.
- When building roads, consider flatter routes for better transport efficiency.


## Roadmap Highlights

- Expose APIs for Python and Lua scripting to automate tests:
  - Terrain/water access, relaxers, and river extraction
  - Road building and path profiles
  - Vehicles stepping and state queries
  - Pathfinding solve steps and path retrieval
  - Economy: load technologies, create factories, set stocks, produce

- Integrate pathfinding for slope/terrain-aware road planning.
- Couple transport with economy to simulate supply chains.
