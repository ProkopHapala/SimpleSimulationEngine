# Radiosity — Radiative Heat Transfer Solver

Steady-state radiative heat transfer between surface elements using the radiosity method. 2D and 3D versions.

## Physics

Lambertian diffuse emission + reflection, Stefan-Boltzmann law (`P = A·ε·σ·T⁴`), view factors `F_ij`. Two-sided surfaces with per-side emissivity (`eps_front`, `eps_back`), single heat source per element. Linear system:

```
(I - W)·B = heat_in / area
```

solved via `numpy.linalg.solve`. Temperature recovered from `T⁴`.

## Files

- **`Radiosity.py`** — 2D solver: discretize line segments into elements, compute view factors (NumPy vectorized), solve, visualize. Supports two-sided surfaces with per-side emissivities.
- **`Radiosity3D.py`** — 3D solver: rasterize polygon faces into hex elements (reuses `PolygonRasterization.py`), compute 3D view factors, solve, visualize in 3D. Input: OBJ-style (V, F) mesh.
- **`TriangleOcclusionRaytracer.py`** — GPU occlusion: PyOpenCL ray-triangle intersection to build N×N binary occlusion matrix. Uses shared kernel `cpp/common_resources/cl/radiosity.cl` (`occlusion_matrix` kernel). One work-item per element i, loops over all j and all triangles. Brute-force O(N²·T).
- **`PolygonRasterization.py`** — Hex rasterizer for 2D polygons (used by 3D solver to discretize faces).
- **`OctahedralSphereMaping.py`** — Octahedral sphere mapping utility.

## Pipeline

1. Discretize surfaces into elements (line segments in 2D, hex tiles per face in 3D)
2. Compute geometry-only view factor matrix `F` (N×N) — midpoint inverse-square kernel with Lambert cosines
3. Compute occlusion matrix on GPU (ray-triangle, binary visibility)
4. Mask: `F *= (1 - occ)` to zero out blocked pairs
5. Solve `(I - W)·B = P` for radiosity `B`
6. Recover `T⁴` from `B` via thin-sheet energy balance

## Key Design Points

- Two-sided surfaces: single DOF per element (N unknowns, not 2N), per-side emissivities applied on reception side
- View factor: `F_ij ≈ (|n_i·r̂_ij| |n_j·r̂_ji| / (π r²)) * A_j`
- Occlusion: brute-force, no acceleration structure yet (BVH needed for large N)
- Vectorized NumPy for view factors; GPU for occlusion

## Run

```bash
# 2D
python Radiosity.py --labels temp --show-matrix

# 3D
python Radiosity3D.py --hex 0.35 --epsFront 0.8 --epsBack 0.8 --hotFace 0 --P_hot 100

# Occlusion only
python TriangleOcclusionRaytracer.py --geom three_quads --elements
```

## Shared Infrastructure

- Occlusion kernel: `cpp/common_resources/cl/radiosity.cl` — shared with pyScatter
- Hex rasterizer: `PolygonRasterization.py` — reused by pyScatter's `SimulationPipeline.py`
- Radiosity solver functions reused by `pyScatter/SimulationPipeline.py`
