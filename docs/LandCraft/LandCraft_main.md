# LandCraft_main.cpp - Main Game Controller

## Purpose
The central orchestrator that integrates all LandCraft subsystems into a cohesive simulation game. Manages the game loop, user interaction, and system coordination.

## Core Architecture

### Main Application Class: LandCraftApp
- **Inherits**: AppSDL2OGL - Provides SDL2/OpenGL rendering foundation
- **Purpose**: Coordinates terrain, hydraulic, economic, and transportation systems

### Key Systems Integration

#### 1. Geographic Data Management
- **map_center**: Reference point for the game world
- **ruler**: Hexagonal grid system for terrain representation using SimplexRuler
- **square_ruler**: Alternative rectangular grid via Ruler2DFast for specific calculations
- **maxHeight**: Maximum terrain elevation (500.0 units)

#### 2. Hydraulic Simulation
- **hydraulics**: 2D water flow simulation (HydraulicGrid2D)
- **ground**: Terrain elevation array shared with hydraulics system
- **water**: Water level array for hydraulic calculations
- **hydro1d**: 1D hydraulic simulation for river profiles

#### 3. Economic Framework
- **commodities**: Hash map of all tradeable goods (Commodity objects)
- **technologies**: Hash map of production technologies (Technology objects)
- Loaded from "data/Technologies.txt" via loadTechnologies()

#### 4. Transportation Network
- **roadBuilder**: Tool for constructing roads between points
- **roads**: Collection of all constructed roads
- **vehicleTypes**: Different vehicle specifications for transportation
- **vehicles**: Active vehicles moving along road networks

### Game Loop Functions

#### Initialization
- **makeMap()**: Creates terrain grid with specified size and resolution
- **generateTerrain()**: Procedural generation using noise and erosion
- **makeRivers()**: Identifies and creates river systems
- **makeRoads()**: Constructs initial road networks
- **makeVehicles()**: Populates roads with transport vehicles

#### Update Cycle
- **hydroRelaxUpdate()**: Advances 2D hydraulic simulation
- **hydro1D_update()**: Advances 1D river flow simulation
- **updateVehicles()**: Moves vehicles along road networks

#### Rendering Pipeline
- **draw()**: Main render pass showing terrain, water, roads, and vehicles
- **drawHUD()**: Interface elements including river/road profiles
- **drawTerrain()**: Hexagonal grid rendering with elevation-based coloring
- **drawRivers()**: Visualizes river networks with flow-based line widths
- **drawRoad()**: Renders road networks as continuous paths
- **drawVehicles()**: Shows active transport vehicles

### User Interaction

#### Input Handling
- **eventHandling()**: Processes keyboard and mouse events
- **mouseHandling()**: Manages terrain editing via mouse input
- **commandDispatch()**: Routes user commands to appropriate systems

#### Interactive Features
- **Terrain Editing**: Click-drag to modify elevation
- **Water Control**: Add/remove water sources at specific locations
- **Road Construction**: Right-click drag to build roads between points
- **River Analysis**: Select and visualize river profiles
- **Vehicle Monitoring**: Track transport efficiency

### Data Persistence
- **save/load**: Binary storage of terrain and water states
- **Command Scripts**: Automated execution via "data/comands.ini"

### Visualization Tools
- **riverProfile**: 2D plot showing river elevation and flow profiles
- **roadProfile**: 2D plot showing road elevation along routes
- **GUI Panels**: Interactive controls for simulation parameters

## Integration Patterns

### System Coupling
- **Terrain → Hydraulics**: Elevation data drives water flow
- **Hydraulics → Transportation**: Rivers affect optimal road placement
- **Economy → Infrastructure**: Resource needs drive construction priorities
- **Transportation → Economy**: Road networks enable resource distribution

### Real-time Updates
- **Frame-based**: Visual updates synchronized with rendering
- **Physics-based**: Hydraulic simulation runs at fixed time steps
- **Event-driven**: User interactions trigger immediate system responses

This file serves as the central nervous system of LandCraft, ensuring all subsystems work together to create a cohesive geographic economy simulation experience.

## Technical (from code)
The following items are grounded directly in `cpp/apps/LandCraft/LandCraft_main.cpp`:

### Main Application Class: `LandCraftApp`
- Inherits `AppSDL2OGL` (SDL2/OpenGL base).
- Holds systems: `SimplexRuler` hex grid (`ruler`), `Ruler2DFast` square grid (`square_ruler`), `HydraulicGrid2D hydraulics`, `Hydraulics1D hydro1d`, roads/vehicles, basic economy containers, GUI.

### Geographic Data
- `map_center` set from hex grid size; `maxHeight=500.0` used to scale terrain.
- Terrain/water arrays are aliases to `hydraulics.ground`/`hydraulics.water`.

### Initialization sequence (constructor)
- Loads font textures; registers commands/UI layers.
- `loadTechnologies("data/Technologies.txt")` populates `technologies` via `Technology::fromFile`.
- `makeMap(128,50,false)`: sets hex grid, allocates hydraulics, loads ground/water from `data/*.bin` unless `newMap`.
- `makeRivers()`: runs `hydraulics.gatherRain(100)`, `findAllRivers(50)`, builds a GUI dropdown and default plots.
- `makeRoads()`: calls `addRoadStright({10,15},{55,38})`, initializes and renders `roadProfile`.
- `makeVehicles()`: creates one `RoadVehicleType` and one `RoadVehicle` bound to the road.
- Initializes `hydro1d` arrays with noise and constant water; fills random water on 2D grid; calls one `hydraulics.relaxWater({36,39})`.

### Update & draw
- `draw()`: draws terrain via `drawTerrain()`; optionally runs `hydroRelaxUpdate()` [need-check: function body not present in this file]; `hydro1D_update()`; draws rivers/trace/roads/vehicles when toggles are set.
- `drawHUD()`: draws `roadProfile`/`riverProfile` and GUI.

### Interaction
- `eventHandling()`: key-bindings via `Commander` and `CommandList` dispatch to `commandDispatch()`; mouse left edits terrain height under cursor; right button builds a new straight road between press/release.
- `commandDispatch()`: implements actions `generateTerrain`, `relaxWater`, `relaxWaterHex`, `save`, `load`, `outflow`, `inflow`, `findRivers`, `terrainViewMode`, `traceDroplet` by delegating to hydraulics helpers and file IO.

### Rendering helpers
- `drawTerrain()` renders a hex grid using triangle strips and `terrainColor()` modes (height/water and known-river overlay).
- `drawRiver()` renders river segments with width proportional to `flow[ii]`.
- `addRiverPlot()` and `addRoadPlot()` push series into `Plot2D` widgets for inspection.

### Roads & vehicles
- `RoadBuilder.pushStright()` used to generate a hex-grid path; `writeIt()` copies to a contiguous `RoadTile*` array.
- `RoadVehicle.move(t)` advances along the `road->path` using placeholder `getSpeed()` (returns 1.0) and slope from adjacent tile `height` values.

### Persistence & scripts
- `saveBin`/`loadBin` used for `data/ground.bin` and `data/water.bin`.
- `cmdPars.execFile("data/comands.ini")` executes command script lines.

### Integration intent
- Terrain → Hydraulics: elevation drives water flow and river detection. (Supported by `HydraulicGrid2D` usage.)
- Hydraulics → Transportation: rivers can inform road placement [need-check: no direct coupling in code beyond visualization].
- Economy → Infrastructure: production/consumption should drive construction [need-check: economy not yet integrated in `LandCraftApp` beyond loading technologies].
- Transportation → Economy: roads enable distribution [need-check: no commodity flows implemented here].

### Update cadence
- Frame-based rendering updates happen each `draw()`.
- Physics time stepping for hydraulics is toggled by UI flags [need-check: fixed timestep not enforced in this file].

This file currently serves as an integration prototype; several systems are visualized and interactively manipulated, while deeper economic coupling remains to be implemented.