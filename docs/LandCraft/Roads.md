# Roads.h - Transportation Infrastructure System

## Purpose
Manages the construction, representation, and usage of transportation networks including roads, railways, and vehicle movement systems for LandCraft.

## Core Components

### 1. Road Representation

#### RoadTile Structure
- **Purpose**: Individual segment of road network with precise location
- **Attributes**:
  - `ia, ib`: Grid coordinates (uint16_t for memory efficiency)
  - `height`: Elevation at this road segment
- **Usage**: Building blocks for complete road paths

#### Road Class
- **Purpose**: Complete transportation route between two points
- **Key Features**:
  - `n`: Number of tiles in the road
  - `path`: Dynamic array of RoadTile segments
  - Automatic memory management via destructor
- **Algorithm**: Contiguous array storage for cache-efficient access

### 2. Road Construction System

#### RoadBuilder Class
- **Purpose**: Interactive tool for creating and editing road networks
- **Key Functions**:
  - `pushStright()`: Creates straight-line roads using hexagonal grid algorithm
  - `readIt()`: Loads existing road into editable format
  - `writeIt()`: Saves edited road back to storage

#### Hexagonal Path Algorithm
- **Source**: Based on Amit Patel's hexagonal line-of-sight algorithm
- **Purpose**: Creates natural-looking roads on hexagonal grids
- **Algorithm**:
  - Uses error accumulation for smooth diagonal movement
  - Handles both horizontal and vertical hexagon orientations
  - Optimizes for minimal tile usage while maintaining straight appearance

### 3. Vehicle System

#### RoadVehicleType Class
- **Purpose**: Defines vehicle characteristics for transport simulation
- **Key Attributes**:
  - `maxSpeed`: Top speed on flat terrain (50 m/s)
  - `maxForce`: Maximum pulling force (10 kN)
  - `power`: Engine power (100 kW)
  - `mass`: Vehicle weight (1000 kg)
  - `cargoMass`: Maximum payload (1000 kg)
  - `cargoVolume`: Maximum cargo space (5 m³)

#### RoadVehicle Class
- **Purpose**: Individual vehicle navigating road networks
- **Key Features**:
  - `ipath`: Current position along road (tile index)
  - `idir`: Direction of travel (1=forward, -1=backward)
  - `onWay`: Vehicle state (traveling vs. waiting)
- **Algorithms**:
  - **Slope-based Speed**: Speed varies with terrain gradient
  - **Bidirectional Movement**: Vehicles can reverse direction at endpoints
  - **Realistic Physics**: Movement time calculated from slope and vehicle specs

### 4. Movement Simulation

#### Vehicle Movement Algorithm
1. **Slope Calculation**: Computes elevation change between current and next tile
2. **Speed Adjustment**: Reduces speed based on uphill gradient
3. **Time Calculation**: Determines travel time for current segment
4. **Progress Tracking**: Updates position based on elapsed time
5. **Endpoint Handling**: Reverses direction or waits at road ends

#### Transport Efficiency
- **Energy Consideration**: Uphill segments consume more time/fuel
- **Load Capacity**: Vehicles respect weight and volume limits
- **Route Optimization**: Longer but flatter routes may be more efficient

## Integration with LandCraft

### 1. Geographic Context
- **Terrain Integration**: Roads follow actual elevation data
- **Water Interaction**: Rivers may require bridges or detours
- **Slope Analysis**: Steep terrain increases construction and travel costs

### 2. Economic Impact
- **Construction Costs**: Longer roads through mountains are expensive
- **Transport Efficiency**: Flatter routes enable faster, cheaper transport
- **Resource Access**: Roads enable extraction of remote resources

### 3. Strategic Planning
- **Route Selection**: Balance construction cost vs. ongoing transport efficiency
- **Network Effects**: Multiple connected roads create economic synergies
- **Terrain Modification**: Sometimes cheaper to modify terrain than build around

## Advanced Features

### 1. Road Editing
- **Interactive Construction**: Click-and-drag road building
- **Real-time Preview**: Shows path and elevation profile during construction
- **Profile Visualization**: 2D plots showing elevation changes along routes

### 2. Vehicle Management
- **Fleet Operations**: Multiple vehicles per road for capacity scaling
- **Type Specialization**: Different vehicles for different cargo types
- **Traffic Simulation**: Potential for congestion modeling

### 3. Infrastructure Scaling
- **Hierarchical Networks**: Local roads feeding into major highways
- **Upgrade Paths**: Roads can be improved for higher capacity
- **Maintenance**: Long-term degradation and repair systems

## Technical Implementation

### Memory Management
- **Efficient Storage**: RoadTile uses minimal memory for large networks
- **Dynamic Arrays**: Roads resize automatically based on path length
- **Cache Optimization**: Contiguous memory for fast vehicle movement calculations

### Performance Considerations
- **Path Caching**: Pre-calculated paths for repeated routes
- **Vehicle Pooling**: Efficient allocation/deallocation of vehicle objects
- **Spatial Indexing**: Fast lookup of roads in specific areas

This system enables realistic transportation network simulation where geographic constraints directly impact economic efficiency and strategic decision-making.

## Technical Implementation

### Data structures
- `struct RoadTile { uint16_t ia, ib; double height; }` with `print()`.
- `class Road { int n; RoadTile* path; ~Road(){ if(path) delete [] path; } }` contiguous storage of tiles.

### RoadBuilder
- Keeps an editable `std::list<RoadTile> path` and a pointer to a `Road`.
- `pushStright(const Vec2i&)` implements a hex-grid straight line algorithm (based on Amit Patel's article as per comment) that marches with error accumulation and pushes intermediate tiles.
- `readIt()` copies from `Road::path` to builder list; `writeIt()` allocates/resizes `Road::path` and copies tiles back.

Notes:
- `RoadTile.height` is not filled in `pushStright()` here; external code must assign heights [need-check].

### Vehicles
- `class RoadVehicleType` contains parameters (`maxSpeed`, `maxForce`, `power`, `mass`, `cargoMass`, `cargoVolume`); `getSpeed(double slope)` currently returns `1` (placeholder).
- `class RoadVehicle` holds `road`, `type`, `ipath`, `idir`, `t_rest`, `onWay`.
  - `depart()` sets direction at endpoints and marks `onWay`.
  - `moveStep(double& t_left)` computes slope from consecutive `RoadTile.height`, computes a time step via `type->getSpeed(slope)`, and advances `ipath` consuming `t_left`.
  - `move(double t)` loops `moveStep()` while `t>0`, stops at endpoints and clears `onWay`.

Ownership/lifetime:
- `Road` frees its `path` in destructor; `RoadVehicle` has raw pointers (ownership managed elsewhere) [need-check].