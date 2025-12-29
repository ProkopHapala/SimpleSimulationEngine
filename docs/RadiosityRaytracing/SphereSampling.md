# SphereSampling: Overview

## Background
The module `cpp/common/maps/SphereSampling.h` provides a seamless parameterization of the sphere via an icosahedral “diamond” charting. A scalar field on S² (e.g., displacement/bump map) is sampled over 10 rhombi (each rhombus = two triangles) derived from the icosahedron. For any unit vector `p=(x,y,z)`, the code:
- Selects the rhombus sector using azimuth and plane tests (`icosa_edge_planes[]`).
- Selects upper vs. lower triangle via a second plane test.
- Computes barycentric coordinates w.r.t. the triangle’s three spherical vertices (`icosa_polar_verts[]`).
- Converts them into diamond coordinates `(fa, fb)∈[0,1]^2` with controlled flips so that adjacent faces share the same edge orientation.

This enforces continuity across seams (no cracks) and consistent orientation of shared edges between neighbors.

## Icosahedral vs. Octahedral mapping (which is used?)
- __Icosahedral mapping (used, finished for forward mapping/rendering)__
  - Implemented by `sampleIcosa2quads()`, `diTri2cartes()`, and rendered via `drawIcosaMap()`.
  - Provides seamless charts with consistent edge orientation; this is what the demo uses to generate and display heightfields.
- __Octahedral mapping (experimental/incomplete)__
  - Only `sampleOctahedron()` exists, which identifies octants and returns `|x|,|y|` as provisional coordinates.
  - No full seamless parameterization or reconstruction helpers; not used by the renderer.

## Files and Roles
- `cpp/common/maps/SphereSampling.h`
  - Data for icosahedral sector planes (`icosa_edge_planes[]`) and sphere vertices (`icosa_polar_verts[]`).
  - Mapping helpers S² ⇄ rhombus diamond: `sampleIcosa2quads()`, `diTri2cartes()`, `icosa2cartes()`, `sampleTri()`.
  - Minimal octahedral helper `sampleOctahedron()` (octant id only).
- `cpp/common_SDL/SDL2OGL/DrawSphereMap.h`
  - Rendering/tessellation of the icosa-diamond charts: `drawIcosaMap()`, `drawDiTri()`, `drawDiTri_seam()` (and wire variants). Seam routines apply index flips to ensure edge orientation matches between neighboring faces.
- `cpp/sketches_SDL/3D/test_SphereSampling.cpp`
  - Demo app generating a spherical heightfield (noise/craters) and drawing it via `drawIcosaMap()`.

## Functions (API summary)
- `sampleIcosa2quads(Vec3d p, uint8_t& iface, double& a, double& b)` — Map unit vector `p` to rhombus index `iface∈[0..9]` and diamond coords `(a,b)`; uses plane tests and barycentric conversion. Applies flips on the upper triangle so shared edges retain the same parametric direction across faces.
- `getIcosaFace(Vec3d p, uint8_t& iface)` — Determine the icosahedron face (triangle) index using plane tests only; no `(a,b)` returned.
- `template<typename T> diTri2cartes(T fa, T fb, const Vec3T<T>& a, const Vec3T<T>& b, const Vec3T<T>& c, const Vec3T<T>& d, Vec3T<T>& p)` — Reconstruct spherical point `p` from diamond coords `(fa,fb)` and the rhombus corners `(a,b,c,d)`. Piecewise formula preserves continuity across the diagonal.
- `icosa2cartes(Vec2i ns, int iface, float fa, float fb, Vec3d& p)` — Convenience wrapper over `diTri2cartes()` for the given rhombus index, using built-in vertex tables.
- `sampleTri(Vec3d p, Vec3i tri, Vec3d* verts, Vec3d& c)` — Compute (normalized) barycentric weights of `p` w.r.t. triangle `tri` with vertices `verts`.
- `sampleOctahedron(Vec3d p, uint8_t& iface, double& a, double& b)` — Identify octant and return `|x|,|y|` as coords (a minimal octahedral mapping helper; not a full seamless scheme).

Rendering helpers (`cpp/common_SDL/SDL2OGL/DrawSphereMap.h`):
- `drawIcosaMap(Vec2i n, float* heights, float hscale)` — Render the full sphere by tessellating the 10 rhombi; calls seam routines that reverse index order as needed to ensure continuity/orientation.
- `drawDiTri(...)`, `drawDiTri_seam(...)` (+ wire variants) — Per-rhombus tessellation and seam stitching with explicit index “views”.

## How continuity and edge orientation are enforced
- __Sector selection by planes__ — `icosa_edge_planes[]` contains separating plane normals. `getIcosaFace()` and `sampleIcosa2quads()` use dot-products to decide which belt (0..4) and cap (top/bottom) a point belongs to, robustly separating regions.
- __Barycentric on spherical triangles__ — `cs.fromLinearSolution(va,vb,vc,p)` gives weights normalized to sum 1 for the chosen triangle.
- __Diamond coordinate flip__ — In `sampleIcosa2quads()`, for the upper triangle: `a=1-cs.b`, `b=1-cs.a`; for the lower triangle: `a=cs.a`, `b=cs.b`. This makes UV orientation consistent along shared edges.
- __Seam index views__ — `drawDiTri_seam()` passes views like `{n.a,(n.b-1)}`/`{-1,n.b-1}` so that scan directions match on both sides of an edge, guaranteeing equal samples on both sides at identical positions.

## Tutorial
### Evaluate a value at any unit vector p
Suppose each rhombus has an `n={na,nb}` grid and the full buffer `heights` is laid out consecutively per rhombus (block size `na*nb`).
1) Map to diamond coords:
   - `sampleIcosa2quads(p, iface, fa, fb)`.
2) Locate cell and triangle:
   - Convert `(fa,fb)` to integer cell indices `(ia,ib)` and fractional offsets.
   - Decide lower vs. upper cell triangle (e.g., `fa<=fb` is lower; `fa>fb` is upper) mirroring `diTri2cartes()` split.
3) Interpolate:
   - Fetch the 3 grid node values of that triangle from `heights + iface*(na*nb)` and barycentrically interpolate.

Note: The repo does not (yet) expose `evalIcosaGrid(...)`; the above steps mirror how rendering works in `drawDiTri()`.

### Indexing and interpolation details (evaluation side)
- __Per-face block layout__: For grid `n={na,nb}`, each rhombus (face10 index `f∈[0..9]`) has a contiguous block of `nab=na*nb` samples. The base offset is `baseOff = f*nab`. Global buffer length is `total = 10*nab`.
- __Local indices in a cell__: For integer cell indices `(ia,ib)` and fractional offsets `(u,v)` from `(fa,fb)`:
  - `i00 = ia*nb + ib`
  - `i10 = (ia+1)*nb + ib`
  - `i01 = ia*nb + (ib+1)`
  - `i11 = (ia+1)*nb + (ib+1)`
  These must satisfy `0 ≤ i00,i10,i01,i11 < nab` and `0 ≤ baseOff + i11 < total`.
- __Bounds bug (fixed)__: Do not use bitwise OR to validate indices. Use min/max: `imin = min(i00,i10,i01,i11)`, `imax = max(...)`, then assert `imin≥0 && imax<nab`.
- __Triangle split and weights__: Interpolate within the cell consistent with `diTri2cartes()`:
  - If `u ≤ v` (upper diagonal triangle): `h = (1−v)·h00 + (v−u)·h01 + u·h11`.
  - Else (lower): `h = (1−u)·h00 + (u−v)·h10 + v·h11`.
  This ensures evaluation matches how geometry is tessellated.

### Rendering primitives: drawDiTri()
Function: `drawDiTri(Vec2i n, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d, float* hs, double hscale)` in `cpp/common_SDL/SDL2OGL/DrawSphereMap.h`.

- __Purpose__: Render one rhombus (two triangles a–b–c and a–b–d in diamond coordinates) as `GL_TRIANGLE_STRIP`s.
- __Scan pattern__: For each strip `ia = 0..na-2`, it emits a strip over `ib = 0..nb-1` using two rows: `(ia,*)` and `(ia+1,*)`.
- __Positioning__: For each vertex, it converts `(fa,fb)` where `fa=ia/na or (ia+1)/na`, `fb=ib/nb` via `diTri2cartes(fa,fb,a,b,c,d)` to a point on the sphere. Optional normalization (`bNormalize`) and radial relief (`bRelief`, `hscale`) are applied.
- __Sampling__: Heights are fetched as `hs[ia*nb+ib]` and `hs[(ia+1)*nb+ib]`. Colors are derived by `heightColor(h)`.

### Seam stitching and wrapping: drawDiTri_seam()
Function: `drawDiTri_seam(int n, int n2, const Vec3f& a, const Vec3f& b, const Vec3f& c, float* hs, float* hs2, const Vec2i& view, const Vec2i& view2, float hn, float hscale)`.

- __Purpose__: Stitch the common edge between two neighboring rhombi so that vertices and samples align even if the two faces have opposite scan directions along that edge.
- __Edge geometry__: The common edge is parameterized by `a→b` (one face) and advanced toward `c` by `dc=(c−a)/n2`. The routine emits a single strip of `(n+1)` steps along the shared edge.
- __Index views__: The helper `index(i, view) = i*view.x + view.y` remaps 1D indices to read a row/column possibly reversed and with an offset.
  - Example views used by `drawIcosaMap()` for side seams:
    - Left/right edges: `view={n.a, n.b−1}` vs. `view2={−1, n.b−1}`.
    - Top/bottom seams: `view={−1, n.a*n.b−1}` vs. `view2={−n.a, n.a*(n.a−1)}`.
  These flip the readout direction so the two edges traverse samples in the same geometric order.
- __Pole corner value (`hn`)__: At the very last step (`i==n`), one side lacks a sample (corner beyond the last column/row). The code uses an explicit corner height `hn` (usually the first element of the neighbor block) to close the strip without cracks.
- __Wrapping behavior__: There is no circular modulo inside a block. “Wrapping” across faces is achieved by:
  1) Choosing the neighboring face’s buffer pointer (`hs2` vs `hs` or with `+5*nab` for cap faces), and
  2) Passing `view/view2` so that the last column of one face is paired with the first column of its neighbor in the same geometric orientation.

See `drawIcosaMap()` for concrete arguments per seam; it invokes four seams per belt index `i` to connect i↔i+1 and top↔bottom halves.

### Accumulate (inverse projection) at p
To add a value `v` back to the grid:
1) `sampleIcosa2quads(p, iface, fa, fb)`.
2) Choose the triangle within the cell (as in evaluation).
3) Compute barycentric weights and add `v*weight` to the three involved nodes (optionally also accumulate weights for later normalization).

This inverse function is not present yet; implement alongside `evalIcosaGrid` to complete the pipeline.

## Notes & Pitfalls
- Normalize `p` before mapping; clamp indices when on borders; choose a deterministic rule on the diagonal (`fa==fb`).
- Pole/corner stitching uses a single corner value (see `hn` in `drawDiTri_seam()`); mirror this behavior in evaluation to avoid small mismatches.
- Do not change constants in `icosa_edge_planes[]`/`icosa_polar_verts[]` without re-derivation.
- Validate per-face local indices (`0..nab-1`) instead of using bitwise operations; check global bounds `baseOff+i < 10*nab`.

## Status
- __Icosahedral mapping__: Complete for forward mapping and rendering with continuous seams (`sampleIcosa2quads()`, `diTri2cartes()`, `drawIcosaMap()`).
- __Octahedral mapping__: Only a minimal helper (`sampleOctahedron()`), not a full seamless mapping with reconstruction.
- __High-level helpers__: `evalIcosaGrid(...)` and `accumIcosaGrid(...)` are not present yet; they can be added directly on top of `sampleIcosa2quads()` and the diamond split in `diTri2cartes()`.
