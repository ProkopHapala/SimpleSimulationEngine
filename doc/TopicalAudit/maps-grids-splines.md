---
type: TopicalAudit
title: Maps, Grids, Rulers & Splines (2D/3D Lattice Topology, Terrain, City/Castle Generation)
tags: [topic, cpp, grid, ruler, square-lattice, simplex, hexagonal, triangular, spline, hermite, barycentric, terrain, city-generation, castle, pathfinding, truss, fluid2d, nbody, tile, spatial-hash]
---

## Summary

This audit covers the `cpp/common/maps/` directory — the shared library for 2D/3D spatial indexing, lattice rulers (square, triangular/simplex, cubic), spline interpolation on grids, terrain representation and erosion, procedural city/castle generation, pathfinding, truss building, 2D fluid simulation, and N-body worlds. It complements `grids-rulers.md` (which focuses on spatial hashing and collision detection) by covering the **topological and interpolation** aspects of grid systems.

The core abstraction is the **Ruler** — a class that maps between continuous positions and discrete grid indices, with support for interpolation (barycentric on triangles, bicubic Hermite on squares), ray marching, and neighbor enumeration. Multiple lattice types are supported: **square/orthogonal** (`SquareRuler`, `Ruler2DFast`, `Map2D`, `GridMap2D`), **triangular/simplex/hexagonal** (`SimplexRuler`, `SimplexGrid`, `TerrainSimplex`), and **3D cubic** (`CubicRuler`, `GridShape`, `Grid3D`, `grids3D.h`).

## File Inventory

### Rulers & Grid Indexing

| File | Status | Description |
|------|--------|-------------|
| `GridIndex2D.h` | active | Base class: `n`, `ntot`, `i2ip()`, `ip2i()`, `wrap_index()`, `clamp_index()`, `fetchWraped()`. 2D index arithmetic. |
| `Map2D.h` | active | 2D grid stub with origin (`x0,y0`), step, `getIx()`, `getIy()`, `getIndex()`. Base for `TerrainCubic`. |
| `SquareRuler.h` | active | Square-lattice ruler with **bicubic Hermite spline** interpolation (`getValue_cubic`, `getDeriv` via `Spline_Hermite::val2D/dval2D`). Ray march stubs (empty). |
| `Ruler2DFast.h` | active | Fast 2D square ruler with per-axis step. `getOverlapingTiles()` (sphere→cells), `insertSegment()` (Bresenham-like), `insertTriangle()` (scanline up/down pass). |
| `SimplexRuler.h` | active | **Triangular/hexagonal** ruler. Barycentric coords via `toBaricentric`/`fromBaricentric`. `simplexIndex()` → (ia,ib,da,db) with triangle flip. `hexIndex()` for hex cell lookup. `getValue()` with barycentric interpolation. `getDeriv()` with analytic gradient. **Full ray marcher** (`rayStart`, `rayStep`, `rayCut`, `rayHorizonts`) — DDA across 3 edge families (a,b,c). |
| `SimplexGrid.h` | active | Hash-based spatial grid on triangular lattice. `SimplexField<NODE,TILE>` with hi/lo tiles. `simplexIndexBare()`, `nodePoint()`, `tilePoint()`. Uses `HashMap<OBJECT>` base. `raster_line()` stub (commented out). |
| `CubicRuler.h` | active | 3D cubic lattice ruler with arbitrary `Mat3d` basis. `pos2index()` → (iabc, dabc), `index2pos()`. 27-neighbor stencil (faces+edges+vertices). `xyz2id()` 64-bit key. |
| `Grid.h` | active | 3D `GridShape` with `Mat3d` cell vectors, `dCell`, `diCell`. `grid2cartesian()`, `cartesian2grid()`. `init(R, step, bPow2)`. `cut1D()` for line slices. Includes `Fourier.h` dependency. |
| `grids3D.h` | active | `CubeGridRuler` (bounded) and `CubeGridRulerUnbound` (64-bit interleaved key). See `grids-rulers.md` for details. |
| `grids2D.h` | active | `subdivideLoopGrid()` — Loop subdivision surface on square grid (PBC). Weights: `c_para=0.375`, `c_perp=0.125`, `c_on=0.625`, `c_off=0.0625`. |
| `Grid3D.h` | active | Templated 3D grid `Grid3D<TILE,OBJECT,NX,NY,NZ>` with compile-time sizes. `pos2box()`, sphere overlap insertion. |
| `GridUtils.h` | active | `coulombGrid_brute()` — brute-force Coulomb sum between two 3D charge density grids with rotation. |

### Tile Storage & Spatial Hashing

| File | Status | Description |
|------|--------|-------------|
| `GridMap2D.h` | active | Templated 2D grid map `GridMap2D<OBJECT,TILE>`. Specializations for `Segment2d` (Bresenham line rasterization) and `Triangle2d` (scanline triangle rasterization). |
| `HashMap2D.h` | active | 2D spatial hash with bucket indexing. See `spatial-hashing.md`. |
| `HashMap3D.h` | active | 3D unbounded spatial hash. See `grids-rulers.md`. |
| `TileBuffer2D.h` | active | Fixed-size 2D tile buffer `TileBuffer2D<OBJECT,NX,NY,M>`. Compile-time sizes, stack-allocated. Sphere overlap insertion. |
| `TileTree2D.h` | active | Two-level tile tree: `TileTree2D<OBJECT,POWER,NX,NY>` with `LeafTile2D` leaves. Sparse allocation (null tiles skipped). `TileTree2D_d` adds real-space coords. |
| `ArrayMap2D.h` | active | Fixed-size 2D array map `ArrayMap2D<OBJECT,NX,NY>`. Stack-allocated, compile-time dimensions. |
| `Buckets3D.h` | active | 3D bounded cell list with 27-neighbor stencil. See `grids-rulers.md`. |
| `AABBTree3D.h` | active | Dynamic AABB tree. See `collision-detection.md`. |
| `kBoxes.h` | active | K-d tree spatial partitioning. See `collision-detection.md`. |
| `DistanceHierarchy.h` | active | Distance hierarchy for nearest-neighbor queries. |
| `SphereTreeND.h` | active | N-dimensional sphere tree. |
| `SphereTreeND_2.h` | active | Alternative N-d sphere tree implementation. |
| `clustering3D.h` | active | `ClusterMap` (distance-based clustering) and `ClusterBox` (growing bounding box clusters). |
| `temp/MapTile2D.h` | experimental | Marching-squares-like contour rendering: `render_block()` with 16-case lookup table for triangulating quad cells based on corner values. OpenGL display list based. |
| `temp/GridMap2D_Line.h` | experimental | Older `GridMap2D` with `forward_list` storage, `insertLine()` (scanline), `makeStatic()` (freeze to array). Has bugs (undefined variables). |
| `temp/Map2D_grid.h` | experimental | Minimal grid map stub. |
| `temp/MapCell2D.h` | experimental | Single-value cell struct. |
| `temp/MapSite2D.h` | experimental | Minimal site struct. |
| `temp/GridMap2D.h` | experimental | Older grid map variant. |
| `temp/HashMap.h` | experimental | Generic hash map used by `SimplexGrid`. |

### Terrain

| File | Status | Description |
|------|--------|-------------|
| `TerrainCubic.h/.cpp` | active | Terrain on square grid with bicubic Hermite interpolation. `getVal(x,y)`, `rayLine()` (horizon raytracing), `renderRect()`. Extends `Map2D`. |
| `TerrainSimplex.h/.cpp` | active | Terrain on **triangular/simplex** grid. Same simplex indexing as `SimplexRuler`. Hydraulic erosion: `flow_errosion_step()`, `rain_and_evaporation()`, droplet erosion (`droplet_step()`, `errodeDroples()`). `genTerrainNoise()` with multi-octave noise. `raster_line()` for ray traversal. |
| `TerrainHydraulics.h/.cpp` | active | `HydraulicGrid2D` — water flow simulation on square grid. `relaxWaterRasterX/Y()` (water leveling), `gatherRain()`, `findAllRivers()`. River tracing via `PathFinder`-like contour propagation. `Lake` detection. Extends `Grid2DAlg`. |
| `TerrainRBF.h/.cpp` | active | Terrain via **Radial Basis Functions** (RBF). `TerrainRBF` extends `HashMap2D<RBF2D>`. `getVal()`, `insertRBF()`, `generateRandom()`. Spatial hashing for RBF lookup. |

### City & Castle Generation

| File | Status | Description |
|------|--------|-------------|
| `CityGeneration.h` | active | `QuadNode` — recursive quad-tree subdivision for city blocks. `split(n, jitter)`, `makeSubPoints()` (bilinear interp with jitter), `makeSubQuads()`, `inset()` (shrink quad for roads), `splitRecursive()` with random subdivision counts and insets. |
| `CityGeneration2.h` | active | `Quad2d`, `QuadSpliter<Fcond,Fleaf>` (template binary splitter with condition/leaf callbacks), `QuadNode2`, `splitOpen()` (strip-based subdivision with road widths). `drawQuadNode2Rec()` for visualization. |

### Pathfinding

| File | Status | Description |
|------|--------|-------------|
| `PathFinder.h` | active | `PathFinder` extends `Grid2DAlg`. Multi-source BFS/Dijkstra with height cost (`ch2*dh² + chplus*dh` or `chminus*dh`). `path_step()` contour propagation. `findConnections()` — watershed boundary pass detection. `track()` — backtrack path. `makePaths()`. Supports 4/6/8 neighbors. |
| `Grid2DAlgs.h/.cpp` | active | `Grid2DAlg` base with neighbor stencil. `initNeighsSquareN(4 or 8)`, `initNeighs_6()` (hex-like 6-neighbor on square grid), `initNeighsSquareMask()`. `bisecNoise()`, `bisecPaternNoise()` for terrain generation. |

### Truss Building (3D grid-based structures)

| File | Status | Description |
|------|--------|-------------|
| `TrussBuilder.h/.cpp` | active | Grid-based truss/structure builder. `GridNode` with (ix,iy,iz) and exist/fixed flags. `insertNode()`, `insertBond()`, `insertBox()` with mask. `removeNodesWithoutBond()`. `toSoftBody()` conversion. `toFile()`/`fromFile()`. 64-bit node IDs. Used by `test_TrussBuilder.cpp` and castle/structure generation. |

### Fluid & N-Body (grid-based simulation)

| File | Status | Description |
|------|--------|-------------|
| `Fluid2D.h/.cpp` | active | 2D stable fluid (Stam-style). `fluidStep_orig/simplified/minimal()`. Diffuse, advect, pressure projection. `visc`, `diff`, `dt`. Boundary conditions: zero, reflect, absorb, periodic. Extends `Grid2DAlg`. |
| `NBodyWorld2D.h/.cpp` | active | 2D N-body particle simulation with `HashMap2D` spatial hash. `Particle2D` with charge. `pairwiseForce()` (Lennard-Jones-like). String forces. Convergence tracking. |

### Spline Interpolation (shared math)

| File | Status | Description |
|------|--------|-------------|
| `spline_hermite.h` (in `cpp/common/math/`) | active | Cubic Hermite spline: `val()`, `dval()`, `ddval()`, `valdval()`, `valdd()`. 2D bicubic: `val2D()`, `dval2D()` (separable 1D Hermite along x then y, with central-difference derivatives). `basis()`, `dbasis()`, `ddbasis()` basis functions. `Sampler<T>` for streaming interpolation. Used by `SquareRuler`, `TerrainCubic`. |
| `SplineManager.h` (in `cpp/common/math/`) | active | Non-uniform spline manager with control points and optional explicit derivatives. `eval()`, `evalUniform()`, `insertPoint()`, `removePoint()`. Binary search for knot span. |

### Sphere Sampling (angular grids)

| File | Status | Description |
|------|--------|-------------|
| `SphereSampling.h` | active | Icosahedral and octahedral sphere mappings. See `radiosity-raytracing-scattering.md` for detailed coverage. |

## Lattice Types

### 1. Square / Orthogonal Lattice

**Basis vectors**: $\mathbf{a} = (s, 0)$, $\mathbf{b} = (0, s)$

**Index**: `i = ix + nx * iy` (row-major)

**Neighbors**: 4 (von Neumann) or 8 (Moore). `Grid2DAlgs.h` defines `SquareNeighs[8]` with distances (1, √2).

**Interpolation**: Bicubic Hermite spline via `Spline_Hermite::val2D()`. Uses 4×4 stencil (16 points) with central-difference derivatives. `SquareRuler::getValue_cubic()` and `getDeriv()`.

**Ray marching**: `SquareRuler::rayStart/rayStep` are **empty stubs** — not implemented. `Ruler2DFast` has `insertSegment()` (Bresenham) and `insertTriangle()` (scanline) for rasterization but no DDA ray marcher.

**Implementations**: `SquareRuler`, `Ruler2DFast`, `Map2D`, `GridMap2D`, `GridIndex2D`, `TileBuffer2D`, `TileTree2D`, `ArrayMap2D`, `TerrainCubic`, `TerrainHydraulics`, `Fluid2D`, `PathFinder`.

### 2. Triangular / Simplex / Hexagonal Lattice

**Basis vectors**: $\mathbf{a} = (1, 0)$, $\mathbf{b} = (0.5, 0.866025)$ (60° angle)

**Barycentric mapping**: `toBaricentric(p) = (p.x - p.y/√3, p.y * 2/√3)`, `fromBaricentric` is the inverse.

**Index**: `simplexIndex(x,y) → (ia, ib, da, db)` where (ia,ib) is the grid node and (da,db) are barycentric coords. Triangle flip when `da + db > 1.0` (upper triangle of the rhombus).

**Hex cell index**: `hexIndex()` — extends simplex index to 3 cases for hexagonal cells (center + 6 neighbors). Backported from JavaScript `HexGrid.js`.

**Neighbors**: 6 (hex-like) via `Grid2DAlgs::initNeighs_6()` — selects 6 of 8 square neighbors with distance 1.0 (approximation on square grid; true hex neighbors need simplex grid).

**Interpolation**: Barycentric linear interpolation via `trinagleInterp(dind, {h_a, h_b, h_c})`. Analytic gradient via `trinagleDeriv()`. `SimplexRuler::getValue()` and `getDeriv()`.

**Ray marching**: `SimplexRuler` has **full DDA ray marcher** across 3 edge families (a, b, c edges of the triangular grid). `rayStart()` initializes ray in barycentric coords, `rayStep()` advances to next edge crossing, returns `edgeKind` (0=a, 1=c, 2=b). `rayCut()` samples values along ray. `rayHorizonts()` finds horizon intersections (for terrain visibility).

**Implementations**: `SimplexRuler`, `SimplexGrid`, `TerrainSimplex`.

### 3. 3D Cubic Lattice

**Basis**: Arbitrary `Mat3d` (3×3 matrix) — supports non-orthogonal cells.

**Index**: `xyz2id(ix,iy,iz)` = 64-bit interleaved key. `CubicRuler::pos2index()` with `Mat3d invMat`.

**Neighbors**: 27 (face + edge + vertex) via `CubicRuler_neighs[27][3]`.

**Interpolation**: Not implemented for 3D cubic (no `getValue_cubic` equivalent). `GridShape` has `grid2cartesian`/`cartesian2grid` but no interpolation.

**Ray marching**: Not implemented for 3D cubic.

**Implementations**: `CubicRuler`, `GridShape` (`Grid.h`), `Grid3D`, `grids3D.h` (`CubeGridRuler`), `TrussBuilder`.

## Spline Interpolation on Grids

### Bicubic Hermite (Square Grid)

Used by `SquareRuler` and `TerrainCubic`. Implemented in `spline_hermite.h`.

**Method**: Separable application of 1D cubic Hermite along x then y. For each direction, uses 4 points (y₋₁, y₀, y₁, y₂) with central-difference derivatives: `dy₀ = (y₁ - y₋₁)/2`, `dy₁ = (y₂ - y₀)/2`.

**Formula**: `val2D(x, y, f00..f33)` = `val(y, val(x, row1), val(x, row2), ...)` — 4×1D Hermite along x producing 4 intermediate values, then 1×1D Hermite along y.

**Gradient**: `dval2D()` returns both `dfx` and `dfy` by computing derivatives along each axis.

**Stencil**: 4×4 = 16 grid points needed. `SquareRuler::fetchWraped()` handles periodic boundary.

**Status**: Complete and functional. Used for terrain height queries in `TerrainCubic`.

### Barycentric Linear (Triangular Grid)

Used by `SimplexRuler` and `TerrainSimplex`.

**Method**: Each grid rhombus splits into 2 triangles. Point falls in lower triangle (da+db < 1) or upper (da+db > 1). Interpolation: `F = h_a * λ_a + h_b * λ_b + h_c * λ_c` where λ are barycentric coordinates.

**Gradient**: `trinagleDeriv()` computes analytic gradient from the 3 vertex values. Constant per triangle (linear interpolation).

**Status**: Complete and functional. Lower order than bicubic (C⁰ vs C¹) but exact for piecewise-linear surfaces.

### Loop Subdivision (Square Grid)

`grids2D.h::subdivideLoopGrid()` — implements Loop subdivision surface scheme on a square grid with PBC.

**Weights**: `c_para = 0.375` (3/8), `c_perp = 0.125` (1/8), `c_on = 0.625` (5/8), `c_off = 0.0625` (1/16).

**Status**: Implemented but has debug `printf`s. Only works with PBC (periodic boundaries). Non-PBC version is commented out.

### RBF Interpolation (Scattered Data)

`TerrainRBF` uses radial basis functions placed at scattered points, with `HashMap2D` for spatial query of nearby RBFs. Not grid-based interpolation per se, but provides smooth terrain from control points.

### Spline Manager (Non-Uniform Knots)

`SplineManager.h` — general 1D spline with non-uniform knot vector, multiple control point channels, optional explicit derivatives. Not grid-based but used alongside grid interpolation for trajectory/curve evaluation.

## Terrain Representation & Erosion

### TerrainCubic (Square + Bicubic)

- Extends `Map2D` with `heights[]` array
- `getVal(x,y)` — bicubic Hermite interpolation
- `rayLine()` — horizon raytracing (find where ray intersects terrain surface)
- `renderRect()` — OpenGL mesh rendering of terrain patch
- Demo: `test_Horizont.cpp`

### TerrainSimplex (Triangular + Barycentric)

- Same simplex indexing as `SimplexRuler`
- `genTerrainNoise()` — multi-octave noise synthesis
- **Hydraulic erosion**: `flow_errosion_step()` — water flow erodes ground, deposits sediment
- **Droplet erosion**: `droplet_step()` — simulate single water droplet moving downhill, eroding and depositing
- `rain_and_evaporation()` — global water cycle
- `raster_line()` — ray traversal across triangular grid
- Demo: `test_SimplexGrid.cpp`

### TerrainHydraulics (Square + Water Flow)

- `HydraulicGrid2D` extends `Grid2DAlg`
- `ground[]`, `water[]` arrays
- `relaxWaterRasterX/Y()` — iterative water surface leveling (scanline relaxation)
- `gatherRain(minSinkFlow)` — find river sources from flow accumulation
- `findAllRivers(minFlow)` — trace rivers from sources to mouths
- `Lake` detection with `lakeMap[]`
- Used by `LandCraftWorld` for terrain generation and river simulation
- Demo: `test_TerrainHydraulics.cpp`

### TerrainRBF (Scattered RBF)

- `TerrainRBF` extends `HashMap2D<RBF2D>`
- Smooth terrain from placed RBF control points
- `generateRandom()` for testing
- `getVal()` — sum of RBF contributions within query radius

## City & Castle Generation

### CityGeneration.h — Recursive Quad Subdivision

`QuadNode` represents a quadrilateral plot with 4 corner indices into a shared point pool.

**Algorithm**:
1. Start with root quad (4 corner points)
2. `split(n, jitter)` — divide into n×a × n×b sub-quads with jittered split lines
3. `makeSubPoints()` — bilinear interpolation to place sub-grid points with jitter
4. `makeSubQuads()` — create child `QuadNode`s from sub-grid points
5. `inset(t00, t11)` — shrink quad inward (create road margins)
6. `splitRecursive(level, nmax, nmin, jitter, maxInset)` — recurse with random subdivision counts

**Features**: Random plot sizes, road margins via inset, jittered grid lines for organic feel.

**Status**: Functional. Demo in `test_CityGen.cpp`.

### CityGeneration2.h — Binary Splitting & Strip Subdivision

`QuadSpliter<Fcond,Fleaf>` — template-based binary quad splitter:
- `splitBinRec(quad, level, side)` — alternates split direction, calls `fcond()` to decide whether to split further, `fleaf()` for leaf quads
- Random split position in `[cmin, cmax]` (default 0.25–0.75)

`splitOpen()` — strip-based subdivision with road widths:
- Divides quad into n strips along one axis
- Each strip has road+sidewalk+plot layout (s, w parameters)
- `mask` controls which sub-plots to create (left/center/right)

**Status**: Functional. Demo in `test_CityGen.cpp`.

### CastleBuilder App

`cpp/apps/CastleBuilder/CastleWorld.h` — extends `TerrainSimplex`. Currently a stub (empty class body). Intended to combine triangular terrain with castle/structure placement.

### LandCraft App

`cpp/apps/LandCraft/LandCraftWorld.h` — full strategy game world:
- `HydraulicGrid2D hydro` — terrain + water
- `PathFinder pf` — road/path finding
- `RoadBuilder rb` — road construction
- `Economy` — factories, technologies
- `setNeighbors(4/6/8)` — switch between square/hex-like neighbor patterns

### LandTactics App

`cpp/apps/LandTactics/` — tactical combat on grid terrain. Uses `TerrainSimplex` for terrain, `PathFinder` for movement costs.

## Pathfinding

`PathFinder` extends `Grid2DAlg`:

- **Multi-source Dijkstra**: `pepare()` initializes frontier from `centers[]`. `path_step()` propagates frontier to neighbors with movement cost.
- **Height cost**: `heightCost(dh) = ch2 * dh² + (chplus or chminus) * dh` — asymmetric uphill/downhill cost.
- **Terrain cost**: Additional per-cell cost from `terrain_cost[]`.
- **Watershed boundaries**: `findConnections()` detects pass points between basins (minimum cost crossing).
- **Path tracking**: `track()` backtracks from endpoint to center via `toTile[]` pointers.
- **Neighbors**: Supports 4, 6 (hex-like), or 8 neighbors via `Grid2DAlgs`.

## Truss Building (Grid-Based Structures)

`TrussBuilder` — builds 3D truss structures on a cubic grid:

- `GridNode` with (ix, iy, iz) integer coordinates + exist/fixed flags
- `insertNode()`, `insertBond()`, `insertBox(mask)` — add nodes and bonds
- `removeNodesWithoutBond()` — cleanup
- `toSoftBody()` — convert to `SoftBody` for physics simulation
- `toFile()`/`fromFile()` — serialization
- 64-bit node IDs via `xyz2id()`
- Used for building structures (castles, bridges, spaceframes) on regular grids
- Demo: `test_TrussBuilder.cpp`

## 2D Fluid Simulation

`Fluid2D` — Stam-style stable fluids on square grid:

- `fluidStep_orig/simplified/minimal()` — different integration variants
- `diffuse()` — viscous diffusion (Gauss-Seidel relaxation)
- `advect()` — semi-Lagrangian advection
- `pressureBlur()` — pressure projection (Poisson solve via relaxation)
- Boundary conditions: zero, reflect, absorb, periodic
- Extends `Grid2DAlg` with velocity (`vx,vy`) and density (`dens`) fields

## N-Body World 2D

`NBodyWorld2D` — particle simulation with spatial hash:

- `HashMap2D` for neighbor search
- `Particle2D` with charge, position, velocity
- `pairwiseForce()` — Lennard-Jones-like interaction (soft repulsion + attraction)
- String forces between connected particles
- Convergence tracking for relaxation

## Parity Status

- **`SquareRuler` ↔ `SimplexRuler`**: Both are 2D rulers with `getValue()`, `getDeriv()`, `sampleLine()`, `rayCut()`, `rayHorizonts()`. Square uses bicubic Hermite (C¹, 16-point stencil); Simplex uses barycentric linear (C⁰, 3-point). Simplex has full ray marcher; Square has empty stubs. No formal parity test.
- **`TerrainCubic` ↔ `TerrainSimplex`**: Same terrain features (noise gen, erosion, raytracing) on different grids. TerrainSimplex has more complete erosion (droplet, flow, rain). TerrainCubic has horizon raytracing. No cross-grid parity.
- **`CityGeneration.h` ↔ `CityGeneration2.h`**: Different subdivision strategies (recursive n×m vs binary split). No shared interface. No parity.
- **`PathFinder` ↔ `TerrainHydraulics`**: Both use contour propagation on `Grid2DAlg`. PathFinder finds least-cost paths; TerrainHydraulics finds water flow paths. Similar algorithm structure, different cost functions.
- **`GridMap2D<Segment2d>` ↔ `Ruler2DFast::insertSegment`**: Both rasterize line segments onto square grid. `GridMap2D` inserts into tile storage; `Ruler2DFast` returns index list. Same Bresenham-like algorithm.
- **`GridMap2D<Triangle2d>` ↔ `Ruler2DFast::insertTriangle`**: Both rasterize triangles. Same scanline algorithm (up-pass + down-pass).
- **`TrussBuilder` ↔ `CubicRuler`**: Both use 3D integer grid indices. `TrussBuilder` uses `int_fast16_t` coords with `unordered_map`; `CubicRuler` uses `Mat3d` basis with 64-bit IDs. No shared interface.

## Open Issues

- **`SquareRuler::rayStart()`/`rayStep()` empty** — no DDA ray marcher for square grid (SimplexRuler has one). TerrainCubic works around this with `rayLine()` using direct spline evaluation.
- **`GridIndex2D::fetchWraped()` bug** — line 18: `return hs[n.x*ip.y+ip.x]` uses unwrapped `ip` instead of wrapped `ia,ib`. Should be `hs[n.x*ib+ia]`.
- **`SquareRuler::getValue_cubic()` bug** — all 4 rows of the 4×4 stencil pass the same values (`h00,h01,h02,h03` repeated 4 times). Should be `h00..h03` for row 0, `h10..h13` for row 1, etc. The `getDeriv()` has the same bug.
- **`grids2D.h::subdivideLoopGrid()` has debug `printf`s** — should be gated by verbosity flag.
- **`Grid2DAlgs::initNeighs_6()` is approximate** — uses 6 of 8 square neighbors with distance 1.0, but true hex neighbors on square grid have mixed distances (1.0 and √2). Only correct on actual triangular/simplex grid.
- **`temp/GridMap2D_Line.h::insertLine()` has undefined variables** — `xa, ya, xb, yb, ixa, iya` used but not defined (should be `ax, ay, bx, by, ix, iy`). Broken code.
- **`temp/MapTile2D.h::getIndex()` bug** — `(iy<<pow) + iy` should be `(iy<<pow) + ix`.
- **`CityGeneration.h` has no road network generation** — only plot subdivision. `RoadBrancher` class is commented out stub.
- **`CastleWorld` is empty stub** — extends `TerrainSimplex` but adds nothing.
- **No hexagonal grid terrain** — `TerrainSimplex` uses triangular grid, not true hexagonal cells. `hexIndex()` exists in `SimplexRuler` but no terrain class uses hex cells.
- **No 3D grid interpolation** — `CubicRuler` and `GridShape` lack `getValue()`/`getDeriv()`. Only 2D rulers have interpolation.
- **No 3D ray marcher** — `CubicRuler` has no ray marching. Only `SimplexRuler` (2D) has full DDA.
- **`SimplexGrid::raster_line()` commented out** — triangular grid line rasterization not implemented.
- **`TerrainRBF` not integrated with terrain pipeline** — standalone, not used by `LandCraft` or `LandTactics`.
- **`Fluid2D` not integrated with terrain** — standalone fluid sim, no coupling with `TerrainHydraulics`.
- **`NBodyWorld2D` tightly coupled to `Particle2D`** — TODO comment says should be generalized to any `PointBody2D` subclass.

## Related Audits

- **`grids-rulers.md`** — Spatial partitioning for collision detection and neighbor search. Covers `CubeGridRuler`, `HashMap3D`, `Buckets3D`, `HashMap2D` from the collision-detection perspective. This audit covers the same files from the terrain/interpolation/city-generation perspective.
- **`spatial-hashing.md`** — Hash-based spatial indexing (`HashMap64`, `HashMap3D`, `HashMap2D`, `BroadSpaceMapHash`). `SimplexGrid` and `TerrainRBF` use hash maps from this infrastructure.
- **`collision-detection.md`** — `AABBTree3D`, `kBoxes`, `BroadSpaceMapHash`. `Ruler2DFast::insertSegment/Triangle` and `GridMap2D` specializations are used for 2D collision rasterization.
- **`parallel-particle-cell.md`** — Grid infrastructure for N-body and MD. `NBodyWorld2D` and `Buckets3D` are related. GPU cell-list acceleration.
- **`radiosity-raytracing-scattering.md`** — `SphereSampling.h` (icosahedral/octahedral sphere maps) lives in `cpp/common/maps/` but covers angular sampling for radiation transport.
- **`soft-body-truss-dynamics.md`** — `TrussBuilder` converts to `SoftBody` for truss physics simulation.
- **`fluid-dynamics.md`** — `Fluid2D` (Stam stable fluids) and `TerrainHydraulics` (water flow on terrain) are 2D fluid solvers in this directory.
- **`land-tactics.md`** — `LandTactics` and `LandCraft` apps use `TerrainSimplex`, `TerrainHydraulics`, `PathFinder`, and city generation for strategy game terrain.
