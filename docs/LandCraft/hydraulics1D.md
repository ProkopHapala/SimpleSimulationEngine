# hydraulics1D.h - 1D Hydraulic Simulation for River Profiles

## Purpose
Provides simplified 1D hydraulic simulation for analyzing river cross-sections and longitudinal profiles in LandCraft. Complements the 2D hydraulic system for detailed river analysis and water management planning.

## Core Components

### 1. 1D Terrain Generation

#### bisectNoise1D Function
- **Purpose**: Generates 1D terrain profiles using fractal subdivision
- **Algorithm**: Recursive midpoint displacement with random offsets
- **Parameters**:
  - `npow`: Power-of-two resolution (2^npow points)
  - `hs`: Output elevation array
  - `frndMin/Max`: Random displacement bounds
- **Process**:
  1. Start with endpoints at specified elevations
  2. Recursively subdivide intervals
  3. Add random offsets at each subdivision level
  4. Creates natural-looking terrain profiles

### 2. 1D Water Flow Simulation

#### Hydraulics1D Class
- **Purpose**: Simplified 1D water flow along terrain profile
- **Key Data**:
  - `ground`: Terrain elevation array
  - `water`: Water surface elevation array
  - `n`: Number of simulation points

#### Water Balance Algorithm
- **Method**: `step(rain, evapor)`
- **Purpose**: Simulates water redistribution along terrain profile
- **Physics**: 
  - Water flows from high to low elevation
  - Conservation of water volume
  - Gravity-driven flow with equilibrium states

#### Process Flow
1. **Input Calculation**: Apply rainfall and evaporation to each cell
2. **Elevation Analysis**: Calculate elevation differences between adjacent cells
3. **Water Redistribution**: Move water to achieve local equilibrium
4. **Volume Conservation**: Ensure total water mass remains constant

#### Equilibrium States
- **Uphill Flow**: Water accumulates in lower elevation cells
- **Downhill Flow**: Water spreads to achieve level surfaces
- **Partial Filling**: Cells fill to intermediate levels based on total water

### 3. Advanced Water Management

#### Deep Acceleration
- **Method**: `deepAccel(depth)`
- **Purpose**: Simulates rapid water level equalization over large areas
- **Algorithm**:
  1. Identify contiguous areas below target depth
  2. Calculate total water volume in affected region
  3. Redistribute water to achieve uniform level
  4. Update terrain and water arrays accordingly

#### Applications
- **Reservoir Filling**: Simulate dam construction and lake formation
- **Flood Analysis**: Model rapid water level changes
- **Canal Construction**: Predict water redistribution from major projects

### 4. Utility Functions

#### Memory Management
- **realloc(n_)**: Dynamically resize simulation arrays
- **clear()**: Reset all values to zero for new scenarios
- **Efficient Storage**: Contiguous arrays for cache performance
- **Manual Setup**: Direct assignment of terrain profiles
- **Random Generation**: Use bisectNoise1D for natural terrain
- **Import Data**: Load real-world elevation profiles

## Integration with LandCraft

### 1. River Profile Analysis
- **Cross-sections**: Analyze specific river segments
- **Longitudinal Profiles**: Study river gradient and flow patterns
- **Flood Modeling**: Predict water levels during extreme events

### 2. Infrastructure Planning
- **Dam Design**: Calculate reservoir capacity and impact
- **Canal Routing**: Determine water levels for navigation
- **Bridge Planning**: Establish clearance requirements

### 3. Environmental Impact
- **Ecosystem Changes**: Model habitat modification from water projects
- **Sediment Transport**: Track material movement along river profiles
- **Seasonal Variations**: Simulate wet/dry season water level changes

## Technical Characteristics

### Performance
- **Algorithm Complexity**: O(n) per simulation step
- **Memory Usage**: O(n) for terrain and water arrays
- **Real-time Capability**: Suitable for interactive analysis

### Accuracy vs. Speed
- **Simplified Physics**: Sacrifices 2D complexity for 1D speed
- **Valid Applications**: Longitudinal river profiles, cross-sections
- **Limitations**: Cannot model lateral flow or complex 2D phenomena

### Visualization Integration
- **Profile Plots**: 2D elevation vs. distance graphs
- **Real-time Updates**: Water levels update during simulation
- **Interactive Editing**: Direct modification of terrain and water levels

## Strategic Gameplay Applications

### 1. Water Resource Management
- **Hydroelectric Power**: Calculate potential energy from river gradients
- **Irrigation Systems**: Design water distribution networks
- **Flood Control**: Plan protective infrastructure

### 2. Transportation Planning
- **River Navigation**: Determine water depths for boat traffic
- **Port Location**: Identify suitable sites for water transport
- **Lock Systems**: Plan elevation changes for canal navigation

### 3. Economic Optimization
- **Resource Extraction**: Plan access to water-dependent resources
- **Settlement Location**: Choose sites with reliable water access
- **Agricultural Planning**: Design irrigation-dependent farming systems

## Comparison with 2D System

### 1D Advantages
- **Computational Speed**: Much faster than full 2D simulation
- **Profile Analysis**: Ideal for river longitudinal sections
- **Interactive Editing**: Real-time modification and analysis
- **Simplified Understanding**: Easier to interpret results

### 2D System Complement
- **Detailed Modeling**: 2D system handles complex flow patterns
- **Areal Analysis**: Full watershed simulation
- **Interactive Terrain**: Comprehensive terrain modification
- **Network Effects**: Multiple rivers and drainage systems

This 1D system provides the analytical depth needed for specific infrastructure projects while maintaining the speed required for interactive gameplay in LandCraft.
