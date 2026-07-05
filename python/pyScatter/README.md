# pyScatter — X-ray / Neutron Radiation Scattering Simulation

Monte Carlo and deterministic radiation transport (X-ray, neutron, gamma) through spacecraft geometry composed of sparse thin materials (slender tubes, thin metallic sheets).

## Physics

- **Attenuation**: Beer-Lambert law `I = exp(-Σ μ·L)` along ray through intersected primitives
- **Scattering**: Monte Carlo path sampling Source→Scatterer→Detector with transmission `T = exp(-Σ Σ_t·thickness)` per segment (design doc, not fully implemented)
- **Radiosity coupling**: reuses Radiosity solver for thermal radiative heat transfer on same geometry

## Files

- **`xray_sim.py`** — Deterministic X-ray projection: source→detector rays, accumulate optical depth through all intersected primitives (triangles + tubes). Spatial clustering + local memory caching + cone-based broad-phase culling for GPU acceleration.
- **`SimulationPipeline.py`** — Full radiosity pipeline for spacecraft geometry: imports Radiosity modules (`rasterize_face_to_elems`, `assemble_elements`), generates sparse visibility channels on GPU (`compute_sparse_channels`), solves sparse fixed-point radiosity (`solve_radiosity_sparse`, `calculate_temperatures_sparse`). No dense N×N matrices.
- **`SpacecraftGeometry.py`** — Procedural "Christmas-tree" spacecraft geometry generator (vertices + triangle faces).
- **`AttenuationTest.py`** — Standalone attenuation test.
- **`run_all.py`** — Runner script for all sub-simulations.
- **`cl/cl/Xray.cl`** — GPU kernel: `xray_simulation` — per-pixel ray traversal with cluster-based local memory caching, cone culling, tube + triangle intersection, Beer-Lambert attenuation.
- **`cl/SpacecraftSimulation.cl`** — GPU kernels: `compute_occlusion` (brute-force N², legacy), `compute_occlusion_tiled` (tiled cluster occlusion, legacy dense output), `compute_sparse_channels` (sparse top-K visible channels with Xray-style cluster culling + local memory), `compute_attenuation` (per-detector-pixel ray-triangle attenuation).
- **`cl/cl/Scatter.cl`** — Skeleton kernel for Monte Carlo scattering (incomplete).
- **`MonteCarloScattere.md`** — Design document: MC path sampling, transmission, scattering cross-sections, GPU thread layout, BVH acceleration strategy.

## Three Approaches

1. **Deterministic attenuation** (`xray_sim.py` + `Xray.cl`): Source→Detector, sum `μ·L` through all hits. Accelerated with spatial clustering + local memory + cone culling. One thread per detector pixel.

2. **Sparse radiosity channels** (`SimulationPipeline.py` + `SpacecraftSimulation.cl::compute_sparse_channels`): Per-receiver top-K visible channels with geometric view-factor weights. Xray-style cluster broad-phase culling + local-memory triangle cache. Output is sparse `(idx, weight)` arrays of shape `(N, K)`, not dense N×N. Occlusion tested via ray-triangle against clustered triangles. Scales as O(N·K·C) where K=channels, C=clusters.

3. **Monte Carlo scattering** (`MonteCarloScattere.md`): Source→Scatterer→Detector with transmission. Design stage only — `Scatter.cl` is a skeleton.

## Key Design Points

- Primitives: tubes (slender cylinders, 1D edges with radius) + triangles (thin sheets)
- GPU acceleration: spatial clustering of primitives into bounding spheres, cone-based broad-phase culling per workgroup, cooperative local-memory loading of active clusters
- BVH identified as needed for scalability (currently brute-force or simple clustering)
- Reuses Radiosity modules for surface discretization, view factors, and linear solve

## Run

```bash
export PYOPENCL_CTX='1'
python xray_sim.py --width 800 --height 600 --tubes 200 --tris 200
python SimulationPipeline.py
python AttenuationTest.py
```

## End-to-End Test Results

Tested with real spacecraft geometry from `tests_bash/Orbital/ship_ICF_marksman_2.obj` (1271 vertices, 4268 edges, **zero faces** — wireframe/truss model, ~200×400×1100 meters).

### X-ray attenuation: ✅ Works on full ship

- Loaded 4268 edges as tube primitives (1m radius, μ=0.5)
- 67 clusters, 400×300 detector
- Kernel: 0.006s on Intel GPU
- Output: `OUT-xray_ship.png` — spacecraft silhouette with attenuation
- **Scales because it's M×K** (detector pixels × clusters), not N²

### Sparse radiosity channels: ✅ Works on subsets, designed for full ship

Replaced brute-force N² dense occlusion + dense view-factor matrices with **sparse top-K visible channels**:

- **Kernel**: `compute_sparse_channels` — per-receiver top-K candidates by geometric view-factor weight, occlusion-tested via Xray-style cluster culling + local-memory triangle cache
- **Output**: `(idx, weight)` arrays of shape `(N, K)` — sparse, no N×N matrix ever materialized
- **Memory**: O(N·K) instead of O(N²)
- **Row chunking**: Kernel launched in bounded `row_chunk=256` batches to avoid GPU watchdog timeouts on large geometries
- **Solver**: Sparse fixed-point (Jacobi) iteration, not dense `np.linalg.solve`

#### Test: `high_res_mesh.obj` (84 faces)

- 120 elements, 84 triangles, 3 clusters, kmax=32
- Sparse channels: 2759 / 3840 nnz
- Max T⁴ = 68.59 (vs 68.82 with previous dense method — close agreement)
- No dense matrix built

#### Test: Ship subset (1000 edges → 2000 strip triangles)

- 2000 elements, 2000 triangles, 63 clusters, kmax=32
- Sparse channels: 63752 / 64000 nnz
- Max T⁴ = 63.13
- Runtime: 16.73s
- Memory: 64K sparse slots instead of 4M dense matrix entries

### Also tested with `high_res_mesh.obj` (84 faces) — legacy dense path

- Radiosity (dense): ✅ 120 elements, 72% pairs occluded, max T⁴=68.82
- Attenuation: ✅ 84 triangles, but source/detector positions hardcoded — need geometry-aware placement

### What's missing for full spacecraft support

1. **Full-ship sparse run** — row chunking is in place but not yet tested on full 4268-edge ship
2. **Morton/Z-order spatial sorting** — current cluster sort is simple z+x; Morton code would improve cluster locality
3. **GPU sparse solver** — current `solve_radiosity_sparse` is CPU gather; should move to GPU for large N
4. **Geometry-aware source/detector placement** — attenuation test uses hardcoded positions that may miss the geometry
5. **Monte Carlo scattering** — `Scatter.cl` is a skeleton, not functional

## Shared Infrastructure

- Radiosity solver: imports `Radiosity3D.py` functions (`rasterize_face_to_elems`, `assemble_elements`); sparse solve done locally in `SimulationPipeline.py`
- Hex rasterizer: imports `PolygonRasterization.py` via Radiosity
- Ray-triangle intersection: `SpacecraftSimulation.cl` duplicates the Möller-Trumbore pattern from `radiosity.cl` — should be unified
- Sparse channel kernel reuses `Xray.cl` patterns: cluster spheres, cooperative local-memory load, broad-phase culling
