# PathFinder.h - Optimal Route Planning System

## Purpose
Provides intelligent pathfinding algorithms for infrastructure construction, considering terrain elevation, construction costs, and transportation efficiency in LandCraft's hexagonal grid system.

## Core Architecture

### 1. Pathfinding Foundation

#### PathFinder Class
- **Purpose**: Finds optimal routes between geographic points
- **Inherits**: Grid2DAlg - Provides hexagonal grid infrastructure
- **Key Data**:
  - `height`: Terrain elevation for slope calculations
  - `terrain_cost`: Additional costs (construction difficulty, materials)
  - `moveCosts`: Accumulated cost to reach each cell
  - `toBasin/toTile`: Back-pointers for route reconstruction

### 2. Cost Calculation System

#### Height-based Cost Function
- **Algorithm**: `heightCost(dh)`
- **Purpose**: Quantifies construction difficulty based on elevation changes
- **Formula**: `ch2 * dh² + chplus * dh` (uphill) or `ch2 * dh² - chminus * dh` (downhill)
- **Parameters**:
  - `ch2`: Quadratic term for steepness penalty
  - `chplus`: Additional uphill construction cost
  - `chminus`: Reduced downhill construction cost

#### Multi-factor Optimization
- **Elevation**: Steeper slopes cost exponentially more
- **Terrain Type**: Different surfaces have varying construction costs
- **Distance**: Longer routes accumulate linear costs
- **Directional Bias**: Uphill more expensive than downhill

### 3. Dijkstra-based Pathfinding

#### Wavefront Propagation
- **Algorithm**: Modified Dijkstra's algorithm for hexagonal grids
- **Process**:
  1. **Initialization**: Set infinite costs, zero at starting points
  2. **Wave Expansion**: Process cells in order of increasing cost
  3. **Neighbor Evaluation**: Calculate costs to adjacent hexagons
  4. **Cost Update**: Replace if new path is cheaper
  5. **Path Storage**: Maintain back-pointers for route reconstruction

#### Contour Processing
- **Purpose**: Efficient wavefront expansion using contour lists
- **Mechanism**:
  - `contour1/2`: Double-buffered active cell lists
  - `nContour`: Current wavefront size
  - `path_step()`: Advances wavefront by one iteration

### 4. Network Analysis

#### Connection Discovery
- **Purpose**: Identifies optimal connections between multiple destinations
- **Algorithm**: `findConnections()`
- **Process**:
  1. **Basin Identification**: Each destination creates cost basin
  2. **Boundary Analysis**: Examine edges between adjacent basins
  3. **Optimal Passes**: Find lowest-cost connections between basins
  4. **Path Storage**: Save best routes for later use

#### Symmetric Connection IDs
- **Purpose**: Avoid duplicate connection calculations
- **Algorithm**: `symetric_id(Vec2i)` creates unique IDs for unordered pairs
- **Benefit**: Efficient lookup and comparison of potential connections

### 5. Route Construction

#### Path Reconstruction
- **Method**: `track(Way, i1, i2)`
- **Algorithm**: Follow back-pointers from destination to source
- **Output**: Complete path as sequence of grid coordinates


## Integration with LandCraft

### 1. Road Construction
- Pre-planning: Calculate optimal routes before building
- Cost Estimation: Predict construction expenses for budget planning
- Alternative Routes: Compare multiple options for strategic decisions

### 2. Transportation Networks
- Hub Connections: Find best routes between resource centers
- Network Optimization: Minimize total construction cost for entire network
- Capacity Planning: Identify bottleneck segments for upgrades

### 3. Strategic Planning
- Resource Access: Plan roads to reach remote resources efficiently
- Economic Trade-offs: Balance construction cost vs. ongoing transport efficiency
- Future Expansion: Design networks that accommodate growth

## Advanced Features

### 1. Multi-destination Optimization
- Purpose: Find optimal network connecting multiple points
- Algorithm: Iterative connection of nearest neighbors
- Benefit: Minimizes total network construction cost

### 2. Dynamic Cost Factors
- Terrain Integration: Construction costs vary by surface type
- Seasonal Effects: Weather can modify construction difficulty
- Resource Availability: Local materials reduce construction costs

### 3. Real-time Adaptation
- Path Recalculation: Adjust routes when terrain changes
- Obstacle Avoidance: Reroute around new constructions
- Efficiency Updates: Recalculate when transport costs change
- **Memory Usage**: O(n) for cost arrays and back-pointers
- **Scalability**: Efficient for large maps (1000x1000+ cells)

### Hexagonal Grid Advantages
- **Natural Movement**: Six directions provide smooth directional changes
- **Distance Accuracy**: Hexagonal distance more accurate than rectangular
- **Aesthetic Appeal**: Roads appear more natural and organic

### Optimization Techniques
- **Early Termination**: Stop when all destinations reached
- **Priority Queues**: Efficient selection of next cell to process
- **Spatial Indexing**: Fast neighbor lookups in hexagonal grid

## Strategic Gameplay Impact

### 1. Geographic Constraints
- **Mountain Barriers**: Steep terrain forces longer, expensive routes
- **River Crossings**: Require bridges or longer detours
- **Valley Routes**: Natural corridors provide efficient pathways

### 2. Economic Decisions
- **Direct vs. Indirect**: Shorter steep routes vs. longer gentle ones
- **Investment Planning**: High initial cost vs. long-term savings
- **Risk Assessment**: Route redundancy vs. single-point failures

### 3. Long-term Planning
- **Network Effects**: Early routes influence future expansion
- **Scalability**: Design for increased traffic and capacity
- **Adaptability**: Routes that accommodate terrain modifications

This system transforms geographic analysis into strategic gameplay, where understanding terrain and optimizing routes becomes crucial for economic success in LandCraft.
