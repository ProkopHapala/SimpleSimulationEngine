---
type: TopicalAudit
title: Parallel Particle & Cell Methods
tags: [topic, cpp, opencl, python, nbody, cell-list, grid-ruler, neighbor-search, particle-in-cell]
---

## Summary

Particle simulation infrastructure for N-body problems, molecular dynamics, and truss systems. Core abstraction is `CubeGridRuler` for uniform grid cell indexing (bounded and unbounded variants). Cell-based neighbor search via spatial hash maps. OpenCL parallelization for N-body gravity and truss dynamics. Python PyOpenCL wrappers for GPU particle simulation.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/maps/grids3D.h` | active | `CubeGridRuler` — bounded 3D grid with `pos2box()`, `icell()`, `ixyz2i()`, `overlap_Sphere()` (2×2×2 neighborhood), `overlap_BBox()`. `CubeGridRulerUnbound` — unbounded using 64-bit interleaved key (`NBIT3D=21`, `OFFSET3D`). `Grids3D::insert_SphereOfInfluence()` — templated sphere insertion with edge/corner checks. |
| C++ | `cpp/common/engine/broadPhaseCollision.h` | active | `BroadSpaceMapHash` — uses `CubeGridRuler` + `std::unordered_multimap<int,int>` for cell→object mapping. Paint-based deduplication. See `collision-detection.md`. |
| C++ | `cpp/common/maps/HashMap2D.h` | active | 2D hash map for spatial neighbor search. See `collision-detection.md`. |
| C++ | `cpp/apps/temp/NBodyWorld/NBodyWorld2D.h` / `.cpp` | active | 2D N-body world with collision detection, force assembly, particle update |
| C++ | `cpp/sketches_SDL/2D/test_NBodyColHashMap.cpp` | active | N-body collision using `HashMap2D`: 9-bucket stencil (self + 8 neighbors), pairwise force assembly, boundary reflection |
| C++ | `cpp/common/OCL/OCL_Orb.h` | active | OpenCL truss/particle dynamics: `run_getTrussForces()`, `run_projective_dynamics()`. Buffer management for points/velocities/forces/neighs/params. See `projective-dynamics.md`. |
| OpenCL | `cpp/common_resources/cl/spacecraft.cl` | active | Kernels: `getTrussForces` (neighbor-based force assembly), `updateJacobi_neighs` (per-neighbor constraint projection), `PD_perdictor`/`PD_corrector` (leap-frog time step) |
| Python+OpenCL | `python/GLCL/NBody_glcl.py` | active | Simple O(N²) N-body gravity on GPU. PyOpenCL + PyQt5 OpenGL. Inline kernel source. No cell list — brute force all-pairs. |
| Python+OpenCL | `python/GLCL2/OCLsystem.py` | active | `OCLSystem` — generic OpenCL manager: context, queue, program loading, buffer dict, kernel caching, `BakedKernelCall` for pre-baked kernel execution |
| Python+OpenCL | `python/pyMolecular/OCL/OpenCLBase.py` | active | `OpenCLBase` — base class for OpenCL apps: `load_program()`, `extract_kernel_headers()`, `create_buffer()`, `toGPU()`/`fromGPU()`, `generate_kernel_args()` (auto-parse kernel headers for arg binding). Used by `MolecularDynamics`. |
| Python+OpenCL | `python/pyMolecular/OCL/MolecularDynamics.py` | active | `MolecularDynamics(OpenCLBase)` — multi-system MD with neighbor lists, MMFF force field, grid force field. `allocate_cl_buffers()`, `setup_kernels()`, `run_getMMFFf4()`, `run_runMD()`. Batched per-system execution. |
| Python+OpenCL | `python/GLCL2/GLCLBrowser.py` | active | `GLCLBrowser` — GUI for loading/executing simulation scripts with OpenCL+OpenGL. Script-driven kernel execution, parameter editing, buffer sync CL→GL. |
| Java | `java/Common/CellSort.java` / `CellSort2D.java` | active | Cell-sort collision grids. See `collision-detection.md`. |

## Sub-topics

### Grid Ruler (Cell Indexing)

`CubeGridRuler` (bounded):
- `setup(pmin, pmax, step)` — defines grid extent and resolution
- `pos2box(pos, ipos, dpos)` — position → cell index + fractional position within cell
- `icell(pos)` — position → linear cell index `i = ix + nx*(iy + ny*iz)`
- `overlap_Sphere(pos, r, icells)` — returns up to 8 cells (2×2×2 neighborhood) overlapping sphere
- `overlap_BBox(p0, p1, icells)` — all cells overlapping bounding box
- `ixyz2i_wrap(ip)` — periodic boundary wrapping

`CubeGridRulerUnbound` (unbounded):
- 64-bit interleaved key: `NBIT3D = 21` bits per axis, `OFFSET3D = 2^20`
- `ixyz2long(ip)` — interleaved key: `(ip.x+OFF) + ((ip.y+OFF) + ((ip.z+OFF)<<NBIT)<<NBIT)`
- Range: ±2^20 ≈ ±10^6 (with `step=1`)
- No boundary checking needed

### N-Body Force Assembly

CPU (`test_NBodyColHashMap.cpp`):
- `HashMap2D` cell grid → 9-bucket stencil (self + 8 neighbors)
- `assembleForces()`: onside (within bucket) + offside (neighbor buckets) pairwise
- Boundary: velocity reflection at domain edges

GPU (`NBody_glcl.py`):
- O(N²) brute-force kernel: each particle sums force from all others
- Softening: `dist_sq + 0.1` to avoid singularity
- No cell list — not scalable beyond ~10K particles

GPU (`spacecraft.cl` `getTrussForces`):
- Neighbor-list based: each particle iterates `nNeighMax` neighbors
- `neighs[i*nNeighMax+j]` → neighbor index, `-1` = no neighbor
- `params[j]` = `{l0, kPress, kPull, damping}` per neighbor
- `inv_dt2` scales inertial term

### Molecular Dynamics (Multi-System)

`MolecularDynamics(OpenCLBase)`:
- Multiple molecular systems simulated in parallel (batched)
- Per-system buffers: `apos`, `aforce`, `avel`, `neighs`, `neighCell`, `REQs`, `apars`, `bLs`, `bKs`
- Kernels: `getMMFFf4` (bonded forces), `getNonBond` (non-bonded), `updateAtomsMMFFf4` (integration), `runMD` (full step)
- Grid force field support: `getNonBond_GridFF_Bspline` (texture and non-texture variants)
- Work sizes: `sz_na = (roundup(natoms, nloc), nSystems)`, `nloc=32`

## Parity Status

- **CPU `CubeGridRuler` ↔ GPU cell methods**: CPU has full grid ruler; GPU uses neighbor lists (precomputed on CPU). No direct parity.
- **CPU N-body (`test_NBodyColHashMap`) ↔ GPU N-body (`NBody_glcl.py`)**: Different algorithms (cell list vs brute force). No parity test.
- **CPU truss forces (`evalTrussForce`) ↔ GPU `getTrussForces`**: Same neighbor-list formula. See `projective-dynamics.md` for parity details.
- **`OpenCLBase` ↔ `OCLSystem`**: Two separate OpenCL wrapper classes with overlapping functionality. `OpenCLBase` has kernel header parsing + auto-arg generation. `OCLSystem` has kernel caching + `BakedKernelCall`. No convergence.

## Open Issues

- `NBody_glcl.py` uses O(N²) brute force — no cell list acceleration on GPU
- `CubeGridRuler::overlap_Line()` and `overlap_Triangle()` return `-1` (not implemented)
- `CubeGridRulerUnbound` has no `overlap_Sphere()` — only `overlap_BBox()`
- Two parallel OpenCL wrapper hierarchies (`OpenCLBase` vs `OCLSystem`) — should converge
- `BroadSpaceMapHash` uses VLA `int inds[nIndTmpMax]` on stack — overflow risk for large objects
- No GPU-side cell list construction — all neighbor lists built on CPU then uploaded
- `MolecularDynamics` assumes all systems use same MMFF parameters (`self.mmff_list = [mmff] * nSystems`)
- No periodic boundary support in GPU truss kernels (only in `CubeGridRuler::ixyz2i_wrap`)

## Related Audits

- **`continuum-mechanics-impact.md`** — Eulerian compressible multi-material solver with level set. Uses fixed Cartesian grids (same `CubeGridRuler` infrastructure). Contains a comprehensive **Algorithm Review: Heterogeneous Material Simulation Methods** covering PIC, FLIP, MPM, SPH and their implementation status.
- **`fluid-dynamics.md`** — Overview of all fluid dynamics implementations. References this audit for grid infrastructure.
- **`aerodynamics-hydrodynamics.md`** — Potential flow / vortex methods. Uses particle-based vortex filaments that could benefit from cell-list neighbor search.
- **`soft-body-truss-dynamics.md`** — Truss/mass-spring dynamics. GPU truss kernels (`getTrussForces`, `updateJacobi_neighs`) use the neighbor-list infrastructure documented here.
- **`collision-detection.md`** — Broad-phase collision detection using `BroadSpaceMapHash` (built on `CubeGridRuler`) and `HashMap2D`.
- **`projective-dynamics.md`** — GPU projective dynamics solver using `OCL_Orb` and neighbor-list kernels documented here.
- **`spatial-hashing.md`** / **`grids-rulers.md`** — Detailed documentation of spatial hash maps and grid rulers.
