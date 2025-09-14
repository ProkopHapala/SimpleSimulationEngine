# Ruler2DFast.h - 2D Grid Coordinate System and Spatial Operations

## Purpose
Provides fast and efficient 2D coordinate transformations, spatial indexing, and geometric operations for rectangular grid systems in LandCraft. Serves as the foundation for terrain analysis, construction planning, and spatial queries.

## Core Architecture

### 1. Coordinate System Management

#### Ruler2DFast Class
- **Purpose**: High-performance 2D grid coordinate system
- **Key Variables**:
  - `pos0`: Origin point of coordinate system (0,0 by default)
  - `step`: Grid spacing in x and y directions
  - `invStep`: Precomputed reciprocals for fast division
  - `n`: Grid dimensions (n.x, n.y)
  - `ntot`: Total number of grid cells

#### Index Transformations
- **i2ip(i)**: Converts 1D array index to 2D grid coordinates
- **ip2i(Vec2i)**: Converts 2D coordinates to 1D array index
- **Efficiency**: Uses bit operations for fast integer division

### 2. Spatial Transformations

#### World-to-Grid Mapping
- **x2i(x)**: Converts world x-coordinate to grid index
- **y2i(y)**: Converts world y-coordinate to grid index
- **i2x(ix)**: Converts grid index to world x-coordinate
- **i2y(iy)**: Converts grid index to world y-coordinate

#### Fractional Position Handling
- **x2id(x, dix)**: Returns grid index and fractional remainder
- **y2id(y, diy)**: Returns grid index and fractional remainder
- **Purpose**: Enables sub-grid accuracy for smooth positioning

#### Combined Operations
- **pos2index()**: Converts world position to grid index + fractional offset
- **index2pos()**: Converts grid index + fractional offset to world position

### 3. Spatial Queries

#### Overlap Detection
- **getOverlappingTiles()**: Finds all grid cells overlapping a circular area
- **Algorithm**: Efficient neighbor checking using grid spacing
- **Output**: Array of affected grid indices
- **Use Case**: Terrain modification, construction zones, influence areas

#### Line Segment Intersection
- **insertSegment()**: Rasterizes line segment into grid cells
- **Algorithm**: Bresenham's line algorithm variant
- **Output**: Array of all grid cells intersected by line
- **Use Case**: Road construction, wall placement, boundary definition

#### Triangle Rasterization
- **insertTriangle()**: Fills triangular area with grid cells
- **Algorithm**: Scan-line rasterization with edge sorting
- **Process**:
  1. Sort vertices by y-coordinate
  2. Calculate edge gradients
  3. Scan convert horizontal spans
  4. Output grid cells within triangle

### 4. Grid Geometry

#### Constant Definitions
- **Ruler2D_nEdges**: 4 edges for rectangular cells
- **Ruler2D_nVerts**: 4 vertices for rectangular cells
- **Purpose**: Standardized geometry calculations

#### Grid Properties
- **Uniform Spacing**: Regular grid simplifies calculations
- **Axis-aligned**: Efficient axis-parallel operations
- **Scalable**: Works with any grid resolution

## Integration with LandCraft

### 1. Terrain Analysis
- **Elevation Sampling**: Fast lookup of terrain height at any location
- **Slope Calculation**: Efficient gradient computation using grid neighbors
- **Area Analysis**: Quick summation over rectangular regions

### 2. Construction Planning
- **Building Placement**: Validate locations against terrain constraints
- **Area Calculations**: Determine construction material requirements
- **Overlap Detection**: Prevent building conflicts

### 3. Transportation Networks
- **Route Planning**: Convert world coordinates to grid paths
- **Distance Calculations**: Manhattan or Euclidean distance in grid units
- **Accessibility Analysis**: Check reachability between locations

## Advanced Features

### 1. Performance Optimization
- **Precomputed Values**: `invStep` eliminates division operations
- **Integer Math**: Avoids floating-point precision issues
- **Cache Efficiency**: Contiguous array access patterns

### 2. Flexibility
- **Configurable Origin**: Coordinate system can start anywhere
- **Variable Resolution**: Different step sizes for x and y directions
- **Arbitrary Size**: Supports any grid dimensions

### 3. Precision Control
- **Fractional Coordinates**: Sub-grid accuracy when needed
- **Boundary Handling**: Automatic clipping to grid bounds
- **Wraparound Support**: Optional periodic boundary conditions

## Technical Implementation

### Memory Efficiency
- **Minimal Storage**: Only stores essential parameters
- **No Dynamic Allocation**: All operations use provided arrays
- **Stack-based**: Suitable for real-time applications

### Algorithm Complexity
- **Index Conversion**: O(1) for all coordinate transformations
- **Line Rasterization**: O(n) for n cells intersected
- **Triangle Fill**: O(n²) for n×n triangle area

### Coordinate System Design
- **Right-handed**: Standard mathematical convention
- **Origin at Bottom-left**: Consistent with image coordinates
- **Positive Axes**: Right (x) and up (y) directions

## Comparison with Hexagonal Systems

### Rectangular Advantages
- **Simpler Math**: Straightforward index calculations
- **Axis-aligned**: Efficient for rectangular operations
- **Standard Algorithms**: Bresenham's line, scan conversion

### Hexagonal Trade-offs
- **Neighbor Count**: 4 vs 6 neighbors affects connectivity
- **Distance Metric**: Manhattan vs hexagonal distance
- **Aesthetic Quality**: Rectangular can appear more artificial

## Strategic Gameplay Applications

### 1. Resource Analysis
- **Extraction Zones**: Identify areas within resource deposits
- **Transport Corridors**: Calculate optimal paths for infrastructure
- **Influence Areas**: Determine regions affected by constructions

### 2. Construction Optimization
- **Material Estimation**: Calculate exact quantities needed
- **Labor Planning**: Determine work areas and access routes
- **Quality Control**: Verify construction meets specifications

### 3. Strategic Planning
- **Territory Analysis**: Evaluate control areas and boundaries
- **Defense Planning**: Calculate sight lines and coverage areas
- **Economic Zones**: Define regions for resource management

This system provides the geometric foundation for all spatial operations in LandCraft, enabling efficient calculation and analysis of terrain, construction, and transportation networks.

## Technical (from code)
Grounded in `cpp/common/maps/Ruler2DFast.h`:

### Class and fields
- `class Ruler2DFast` with fields `Vec2d pos0, step, invStep; Vec2i n; int ntot`.

### Indexing and transforms
- `setN(Vec2i)`, `i2ip(int)`, `ip2i(Vec2i)` for 1D/2D index conversion.
- `x2i(double)`, `y2i(double)`, `i2x(double)`, `i2y(double)` world/grid mapping.
- Fractional forms: `x2id`, `y2id`, and delta versions `x2idx`, `y2idy`.
- Combined: `pos2index(pos, dipos, ipos)` and `index2pos(ipos, dipos, pos)`.

### Spatial queries/rasterization
- `getOverlapingTiles(pos, r, results)` returns up to 4 neighboring tiles overlapped by a circle of radius `r` (in tile space), with an optional diagonal.
- `insertSegment(results, a, b)` rasterizes a segment with an integer stepping akin to Bresenham.
- `insertTriangle(results, a, b, c)` scan-converts a triangle by y-sorted edges and horizontal spans.

Constants
- `Ruler2D_nEdges=4`, `Ruler2D_nVerts=4`.

**Notes:**
- Some algorithms assume axis-aligned rectangular grid and uniform `step`.
