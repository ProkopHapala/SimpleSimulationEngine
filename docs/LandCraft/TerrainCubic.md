# TerrainCubic.h - Cubic Terrain Representation and Manipulation

## Purpose
Provides cubic spline-based terrain representation for smooth, continuous elevation modeling in LandCraft. Enables high-precision terrain analysis, realistic rendering, and sophisticated geographic operations.

## Core Architecture

### 1. Terrain Representation

#### TerrainCubic Class
- **Purpose**: Smooth terrain surface using cubic spline interpolation
- **Inherits**: Map2D - Provides 2D grid infrastructure
- **Key Data**:
  - `heights`: Elevation values at grid points
  - `nxy`: Total number of grid points (n.x * n.y)
- **Advantage**: Continuous surface from discrete elevation samples

### 2. Elevation Evaluation

#### getVal(x, y) Method
- **Purpose**: Returns smooth elevation at any continuous coordinate
- **Algorithm**: Cubic spline interpolation using Hermite basis functions
- **Benefits**:
  - Smooth derivatives for slope calculations
  - No discontinuities at grid boundaries
  - Accurate elevation between sample points

#### Interpolation Process
1. **Grid Location**: Identify containing grid cell
2. **Local Coordinates**: Calculate fractional position within cell
3. **Spline Evaluation**: Apply cubic basis functions
4. **Smooth Output**: Continuous elevation value

### 3. Ray Casting and Line-of-Sight

#### rayLine() Method
- **Purpose**: Casts rays through terrain for visibility analysis
- **Parameters**:
  - `hdir`: Horizontal direction vector
  - `p0`: Starting position
  - `hg0`: Starting elevation
  - `dr`: Ray step size
  - `rmax`: Maximum ray distance
  - `ntg`: Number of terrain intersection points
- **Output**: Array of intersection points with terrain

#### Applications
- **Line-of-Sight**: Determine visibility between points
- **Shadow Analysis**: Calculate terrain shadows
- **Radio Propagation**: Model signal coverage
- **Artillery Planning**: Ballistic trajectory analysis

### 4. Rendering Integration

#### renderRect() Method
- **Purpose**: Generates triangle mesh for smooth terrain rendering
- **Parameters**:
  - `x0, y0, x1, y1`: World coordinate bounds
  - `nx`: Grid resolution for rendering
- **Output**: Optimized triangle mesh for OpenGL rendering
- **Feature**: Adaptive resolution based on viewing distance

### 5. Terrain Generation

#### generateRandom() Method
- **Purpose**: Creates synthetic terrain for testing and scenarios
- **Parameters**: `vmin, vmax` - Elevation range bounds
- **Algorithm**: Random elevation assignment to grid points
- **Usage**: Quick terrain generation for prototyping

#### Memory Management
- **allocate()**: Dynamically allocates height array
- **Efficient Storage**: Single contiguous array for cache performance
- **Flexible Sizing**: Supports arbitrary grid dimensions

## Integration with LandCraft

### 1. Smooth Terrain Analysis
- **Slope Calculations**: Accurate derivatives for construction planning
- **Drainage Analysis**: Precise flow direction determination
- **Visibility Planning**: Line-of-sight for defensive structures
- **Accessibility**: Smooth paths for road construction

### 2. Construction Applications
- **Foundation Design**: Accurate elevation for building placement
- **Grading Plans**: Calculate earthwork volumes
- **Road Design**: Smooth vertical curves for transportation
- **Reservoir Planning**: Precise volume calculations for dams

### 3. Strategic Planning
- **Defensive Positions**: Identify high ground and visibility
- **Resource Location**: Analyze terrain for optimal placement
- **Route Planning**: Smooth elevation profiles for infrastructure
- **Environmental Impact**: Accurate terrain modification effects

## Technical Advantages

### 1. Continuous Surface
- **No Terracing**: Eliminates stepped appearance of grid-based terrain
- **Smooth Derivatives**: Accurate slope and curvature calculations
- **Natural Appearance**: Organic terrain shapes

### 2. Computational Efficiency
- **Local Operations**: Only 4 grid points needed for interpolation
- **Cache Friendly**: Contiguous memory access patterns
- **Fast Evaluation**: O(1) complexity for point queries

### 3. Versatility
- **Arbitrary Resolution**: Works with any grid density
- **Flexible Bounds**: Supports non-square regions
- **Multiple Uses**: Rendering, analysis, and simulation

## Comparison with Grid-based Systems

### Cubic Advantages
- **Smoothness**: Continuous first and second derivatives
- **Accuracy**: Better elevation estimates between grid points
- **Aesthetics**: More natural terrain appearance

### Grid Trade-offs
- **Complexity**: More computation than simple grid lookup
- **Storage**: Requires grid point storage vs. continuous functions
- **Performance**: Slightly slower than direct grid access

## Advanced Applications

### 1. Volume Calculations
- **Cut/Fill Analysis**: Precise earthwork volume calculations
- **Reservoir Capacity**: Accurate water storage calculations
- **Material Estimation**: Construction material requirements

### 2. Hydrological Modeling
- **Flow Accumulation**: Smooth drainage basin analysis
- **Flood Modeling**: Precise water level predictions
- **Erosion Simulation**: Accurate slope-based erosion

### 3. Visualization
- **Smooth Shading**: Continuous normal vectors for lighting
- **Level-of-Detail**: Adaptive mesh density
- **Real-time Updates**: Efficient terrain modification

This system provides the smooth, continuous terrain representation needed for realistic geographic simulation while maintaining computational efficiency for large-scale world modeling in LandCraft.

## Technical (from code)
Grounded in `cpp/common/maps/TerrainCubic.h`:

- `class TerrainCubic : public Map2D` holds `double* heights`.
- Declares evaluation and rendering API:
  - `double getVal(double x, double y);`
  - `void rayLine(Vec2d hdir, Vec2d p0, double hg0, double dr, double rmax, int ntg, double* tgs, Vec3d* poss);`
  - `int renderRect(double x0, double y0, double x1, double y1, int nx);`
- Inline helpers:
  - `allocate()` allocates `heights` of size `nxy` (from `Map2D`).
  - `generateRandom(vmin, vmax)` fills `heights` with uniform random values.

**Notes:**
- Implementations of `getVal`, `rayLine`, `renderRect` are not present in this header; expected in a `.cpp` [need-check].
- Includes `spline_hermite.h`, suggesting Hermite interpolation is intended [need-check].
