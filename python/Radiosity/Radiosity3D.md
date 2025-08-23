# Radiosity3D: From 3D Polygons to Hex-Elements and Radiosity

---

## Goal
Build a minimal, fast 3D radiosity pipeline by reusing the existing 2D hex rasterizer (`PolygonRasterization.py`) and radiosity core (`Radiosity.py`). Input is either a Wavefront `.obj` mesh or Python lists of 3D points (vertices) with polygon faces. Each polygon is assumed planar.

## Pipeline Overview
1) Geometry input (OBJ or Python lists)
2) Per-polygon local frame (orthonormal basis)
3) Project polygon to 2D (UV plane)
4) Hex rasterization in UV (reused `HexRasterizerNP`)
5) Map rasterized elements back to 3D (center, area, normal)
6) Assemble element arrays (eps_front/back, rho, heat_in)
7) Reuse radiosity solver from `Radiosity.py` (3D-aware)
8) Visualize in 3D (scatter + optional faces and normals)

## Key Design Points
- Planarity: For each face, take any three non-collinear vertices to build the plane normal `N = normalize((p1-p0) x (p2-p0))`. If degenerate (area~0), skip or raise.
- Local frame: Choose `U` as the normalized projection of the longest edge onto the plane, `V = N x U`. This makes the UV mapping isometric within the plane (unit Jacobian), so areas in UV equal areas in 3D.
- Projection/mapping: For any 3D point `P`, the UV is `u = dot(P-O,U)`, `v = dot(P-O,V)`; back-map `P = O + u*U + v*V`.
- Rasterizer reuse: `HexRasterizerNP(polygon_verts_2d, hex_size)` returns a dict of elements with fields: `area` (2D), `com` (2D), and `geoms` (list of clipped 2D polygons). We create one 3D element per raster element with:
  - `center3 = O + com[0]*U + com[1]*V`
  - `area3 = area2d` (isometry)
  - `normal3 = N`
- Radiosity core reuse: We construct `elements` dict compatible with `Radiosity.py` functions `compute_view_factors`, `solve_radiosity_system`, `calculate_temperatures`:
  - Required keys used by solver: `center (n,3)`, `normal (n,3)`, `area (n,)`, `eps_front (n,)`, `eps_back (n,)`, `rho_front (n,)`, `rho_back (n,)`, `heat_in (n,)`
  - Compatibility: `compute_view_factors` obtains `n = len(elements['length'])`. Provide `length = np.ones(n)` as a harmless placeholder.
- Two-sided surfaces, single DOF per element: Same treatment as in 2D code. We keep one unknown radiosity per element, with per-side epsilons applied on reception side in the transport operator.

## Data Model
- Input mesh:
  - `V`: array (M,3) float64
  - `F`: list of faces, each face is list[int] indices into `V` (>=3 vertices, planar)
- Optional face properties: default emissivities `eps_front`, `eps_back`; optional per-face heating density `P_face` [power per area]. For test, we assign one hot face with uniform `P_face_hot`, others zero.
- Elements (assembled):
  - `center`: (N,3), `normal`: (N,3), `area`: (N,)
  - `eps_front`, `eps_back`, `rho_front=1-eps_front`, `rho_back=1-eps_back`
  - `heat_in`: per element total power = `P_density_element * area`
  - `length`: placeholder ones(N) for compatibility

## Numerical/Debug Rules (repo conventions)
- Vectorized NumPy; avoid Python loops in hot paths
- Float64 arrays; explicit shapes
- Fail loudly with assertions on planarity and frames
- Verbosity flag to gate diagnostics (counts, min/max, etc.)

## Visualization
- Matplotlib 3D:
  - Scatter the element centers colored by T^4 (or T)
  - Optional: draw source mesh faces (wireframe), and normals via `quiver`

## Minimal Test Plan
- Use a unit cube with 6 quad faces
- Rasterize each face with hex size h (e.g., 0.25)
- Make one face hot: `P_face_hot = 100.0` W/m^2; others 0
- eps_front=eps_back=0.8
- Solve radiosity and visualize

## Future Extensions
- OBJ loader (minimal): parse `v`, `f` supporting polygons; convert to 0-based indices
- Visibility/occlusion: not included; current transport is fully geometric Lambertian without shadowing
- Per-face materials; temperature-to-color transfer functions; save figures
- Speed-ups: BVH for near-field subsets, block-sparse assembly, Numba
