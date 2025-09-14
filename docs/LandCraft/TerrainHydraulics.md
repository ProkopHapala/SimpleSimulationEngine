# TerrainHydraulics.h - 2D Water Flow and Terrain Modification System

## Purpose
Provides comprehensive hydraulic simulation for realistic water flow, erosion, and terrain modification in LandCraft. Enables natural river formation, lake creation, and large-scale water management projects.

## Core Components

### 1. Hydraulic Grid System

#### HydraulicGrid2D Class
- **Purpose**: 2D cellular automaton for water flow simulation
- **Inherits**: Grid2DAlg - Provides hexagonal grid foundation
- **Key Data**:
  - `ground`: Terrain elevation array
  - `water`: Water surface elevation array
  - `lakes`: Continuous water bodies for realistic lake simulation

### 2. Water Flow Algorithms

#### Relaxation Methods
- **Purpose**: Simulates water finding equilibrium state
- **Multiple Approaches**:
  - `relaxWaterRasterX/Y`: Row/column-based relaxation for efficiency
  - `relaxWaterHexCells`: Hexagonal neighbor processing for accuracy
  - `relaxWater`: Comprehensive local equilibrium calculation

#### Local Equilibrium Algorithm
1. **Neighbor Collection**: Gathers 6 hexagonal neighbors plus center cell
2. **Water Accounting**: Calculates total water volume in local area
3. **Height Sorting**: Orders cells by ground elevation
4. **Level Calculation**: Determines common water level conserving volume
5. **Distribution**: Assigns new water levels based on sorted heights

### 3. River Network Analysis

#### River Formation
- **Purpose**: Identifies natural drainage patterns from hydraulic simulation
- **Key Classes**:
  - `River`: Represents complete river with path and flow data
  - `Lake`: Defines continuous water bodies at specific elevations

#### Flow Accumulation
- **Algorithm**: `gatherRain()` - Simulates rainfall and tracks water accumulation
- **Process**:
  1. Assign uniform rainfall across terrain
  2. Trace water flow downhill to lowest points (sinks)
  3. Accumulate flow along drainage paths
  4. Identify major rivers based on flow thresholds

#### River Tracing
- **Methods**:
  - `traceDroplet()`: Follows individual water particles downstream
  - `trackRiver()`: Builds complete river from sink to source
  - `findAllRivers()`: Identifies entire river network above minimum flow

### 4. Erosion and Terrain Modification

#### Droplet Erosion
- **Purpose**: Realistic terrain shaping through water flow
- **Algorithm**: `errodeDroples()`
- **Parameters**:
  - Water volume for each droplet
  - Dissolution rate for sediment pickup
  - Sediment capacity for deposition
  - Iteration count for terrain evolution

#### Process Flow
1. **Random Droplets**: Water droplets start at random locations
2. **Downhill Movement**: Droplets follow steepest descent path
3. **Erosion**: Remove material from steep sections
4. **Deposition**: Drop sediment in flatter areas
5. **Terrain Update**: Modify ground elevation based on erosion/deposition

### 5. Water Management Systems

#### Inflow/Outflow Control
- **Purpose**: Enable large-scale water projects (dams, canals)
- **Mechanisms**:
  - `init_outflow()`: Sets boundary conditions for drainage
  - `outflow_step()`: Advances outflow simulation
  - `extend_outflow()`: Expands drainage areas

## Integration with LandCraft

### Geographic Realism
- **Natural Rivers**: Form based on actual terrain elevation
- **Drainage Basins**: Realistic watershed boundaries
- **Water Bodies**: Lakes form naturally in depressions

### Strategic Impact
- **Transportation**: Rivers provide natural highways for boats
- **Resource Access**: Water bodies enable irrigation and power generation
- **Construction Challenges**: Rivers require bridges or detours for roads

### Large-scale Projects
- **Dam Construction**: Modify water flow for power generation
- **Canal Building**: Connect separate water systems
- **Irrigation**: Distribute water for agriculture
- **Flood Control**: Manage seasonal water level changes

## Advanced Features

### Conservation Laws
- **Water Volume**: Total water mass conserved during simulation
- **Sediment Transport**: Erosion/deposition maintains mass balance
- **Energy Minimization**: Water seeks lowest energy state

### Performance Optimization
- **Local Processing**: Only active areas require computation
- **Raster Methods**: Efficient row/column processing for large areas
- **Hexagonal Grids**: Optimal neighbor relationships for natural flow

### Visualization Integration
- **Real-time Updates**: Water levels update during gameplay
- **Flow Visualization**: River networks shown with width proportional to flow
- **Interactive Editing**: Direct terrain modification affects water flow

- **Cache Optimization**: Local processing maintains cache coherence

### Algorithm Complexity
- **Relaxation**: O(n) per iteration for n grid cells
- **River Tracing**: O(m) for m cells in drainage area
- **Erosion**: O(k×n) for k droplets over n iterations

This system enables realistic water simulation that directly impacts strategic gameplay decisions, making terrain modification and water management core aspects of the LandCraft experience.
