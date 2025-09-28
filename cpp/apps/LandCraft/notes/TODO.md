
# Simple

- Floading terrain
- critical angle of repose (sypny uhel)
  - some materials does not have angle of repose
- Edit terrain height - save, load


- contruction using water - build power plant, watermil
- Slopes with respect to sun control plants which can grow there
- Wind controled by terrain (coarse grained)

#### Which facilities depend on terrain?
- Harbor - cover ships from waves
- dam/powerplant - max height differece and gradient with minimal barriers
- irrigation
- channel/watterway with channel lock
- road/train - minimal slope, minimal lengh, minimal cost

- Basins - map of all pixels associated to particular sink
    - each sink can set water level independently
    - each sink knows height of saddlepoints to other sinks (basin connectivity network)
    - How to Do it:
        - random droplet tracking which writes to know[] resp. contour
        - when we have pixel map of basis, we can just trace boundaries of these basins

- River Building
	- downward rainfall gathering
    		1. make sort pixels by ground height
    		2. trace rainfall on each pixel downwards according to height sorted order. After this operation each pixel knows total amout of rainwater which flow throught it per unit of time.
    		3. 
	- upward river tracking 
		-  identify sink-less pixelss (those which does not have lower laying neighbor
		- from each sinkless pixel build 1D river upwards
		- for each pixel find neighbor with highest flow, this is the main branch, if there are other neighbors with flow higher than treshhold, they are other side branches (feeders)

# Modular world initialization and separation of concerns (SoC)

- **Goal**: The World class (`LandCraftWorld`) owns all simulation state and steps. The App (`LandCraft_main.cpp`) owns only GUI, input, and rendering.
- **Init sequence (config-driven, similar to `LandTactics/LTWorld`)**:
  - `W.init(dataPath)`
  - `W.loadTechnologies(data/Technologies.txt)`
  - `W.makeMap(size, step, useCache)`
  - If no cache or `useCache=false`: `W.generateTerrain()` and save caches
  - `W.makeRivers()`
  - `W.makeRoadStraight(seedA, seedB)` (or load from file)
  - `W.makeVehicles()`
- **Config file** (future): `data/world.ini` or `data/world.lua` listing steps with parameters and file paths.
  - map.size, map.step, terrain.cache.on, terrain.cache.paths
  - river.min_flow, river.mode
  - roads.from=auto|file, vehicles.default_type
  - techs.file
- **Runtime update separation**:
  - Simulation steps in world:
    - `W.hydroRelaxStep()`
    - `W.vehiclesStep(dt)`
  - Rendering in app:
    - `drawTerrain()`, `drawRivers()`, `drawRoad()`, `drawVehicles()`
  - App should not modify simulation arrays directly except via well-defined world methods.
- **Cleanup TODOs**:
  - Replace `drawRoad(roads[0])` with `if(!W.roads.empty()) drawRoad(W.roads[0])` in `draw()`.
  - Remove/comment legacy app init functions once world init compiles and runs (`makeMap`, `generateTerrain`, `makeRivers`, `makeRoads`, `makeVehicles`).
  - Strip debug prints/rendering from simulation code paths; keep diagnostics as return values or behind verbosity flags.
