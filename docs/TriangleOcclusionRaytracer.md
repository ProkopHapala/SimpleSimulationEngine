---
description: GPU-accelerated occlusion prototype (COG-based) using PyOpenCL
---

# TriangleOcclusionRaytracer (plan)

Goal: Compute an element-to-element occlusion matrix on GPU using triangle obstacles. Initial version tests visibility only between centers of geometry (COGs) of surface elements.

Key constraints (from repo rules):
- Minimal deps: NumPy + PyOpenCL only; small, testable steps; fail loudly.
- Compute separated from plotting; keep kernels simple first.
- Prefer float64 on CPU, float32 on GPU buffers.

## Data model
- Elements (receivers/emitters): use rasterized elements from `python/Radiosity/Radiosity3D.py` (`rasterize_face_to_elems`), then pack:
  - `points`: float32 (n,4) with xyz = element center; w = face_id (as float for compactness).
- Obstacles: triangles from the base mesh (triangulate each polygon face by a fan around v0):
  - `tris`: float32 (ntris*3,4): sequence of A,B,C; each has xyz vertex; w = surface index (face_id) to allow skipping self.
- Output:
  - `occ`: float32 (n,n): 1.0 if any triangle blocks the segment i→j (simple binary occlusion), else 0.0.

## Kernel (OpenCL)
- One work-item per source element i.
- For each target j over [0..n):
  - Build ray direction d = normalize(pj - pi), local basis (hX,hY) using `getSomeUp()` and cross.
  - For every triangle t:
    - Skip if triangle's surface index equals either endpoint's surface index (avoid self-occlusion).
    - Test if the ray (origin at pi, direction d) intersects triangle in projected 2D via `rayInTriangle(a-pi, b-pi, c-pi, hX, hY)`.
    - Early-exit on first hit.
- Write `occ[i*n + j]` to 0/1.

Notes:
- Start with naive O(n*n*ntris). For small n this is OK. Later add local-memory tiling and broad-phase culling.
- Reuse helper routines from `cpp/common_resources/cl/radiosity.cl` (`getSomeUp`, `originInTriangle`, `rayInTriangle`) adapted inline.

## Host pipeline (Python)
1. Build geometry
   - Use `unit_cube()` as test mesh; triangulate faces (fan).
   - Rasterize faces to elements: `rasterize_face_to_elems(..., hex_size=--size or --hex)`, collect centers and face_ids.
2. Pack GPU buffers (float32):
   - points: (n,4), tris: (ntris*3,4), occ: (n*n,)
3. Compile kernel (inline string) using `OpenCLBase.load_program` from a generated temporary `.cl` file or direct string build.
4. Set `kernel_params` and map buffers via `buffer_dict` names aligned with kernel header.
5. Launch with global size (roundUpGlobalSize(n),), local size `(nloc,)`.
6. Download `occ`, print basic stats; optionally validate symmetry and self-entries.

## Validation
- For unit cube with two visible faces, expect many pairs visible (occ=0) and none occluded unless triangles from other faces are added.
- Sanity:
  - `occ[i*n+i] == 0`.
  - Symmetry not guaranteed by loop order; optionally enforce `occ[j,i]=occ[i,j]` post hoc for diagnostics.

## Next steps (after COG prototype)
- Triangle subdivision per element (use clipped sub-triangles from rasterizer) to compute fractional occlusion.
- Two-stage acceleration (grouped points vs grouped obstacles) as sketched in `radiosity.cl` comments.
- Integrate occlusion into radiosity view factor computation (masking or attenuation).
- Add broad-phase rejection using distance-to-ray bounds and AABBs to reduce triangle checks.

## Hierarchical occlusion: broad-phase and narrow-phase

We will keep both a simple brute-force reference kernel and an optimized staged approach in the shared library `cpp/common_resources/cl/radiosity.cl`.

like [Broad-Phase Collision Detection](https://developer.nvidia.com/gpugems/gpugems3/part-v-physics-simulation/chapter-32-broad-phase-collision-detection-cuda), check it out if there are usefull ideas for GPU-accelerated collision detection.

- __Brute-force (reference)__
  - Kernel: `occlusion_matrix(points, tris, occ)` now lives in `cpp/common_resources/cl/radiosity.cl`.
  - Tests every i→j ray against all triangles; used for debugging/validation and small cases.

- __Broad-phase (group-level)__
  - Kernel: `makeOcclusionMatrix(ns, points, obstacles, distMin, occluders, noccs, Rrange)`.
  - Purpose: for each pair of point-groups (i,j), quickly classify occlusion and collect only likely occluders:
    - If any obstacle is closer than `R_full` to the i→j ray: full occlusion (no need for narrow-phase).
    - If all obstacles are farther than `R_safe`: no occlusion.
    - Otherwise: partial/uncertain — record indices of nearby obstacle-groups into `occluders[(i,j), 0:nOccMax]` and count in `noccs[(i,j)]`.
  - Data model:
    - `points`: float4 (x,y,z,R) for group centers with radius R (grouped elements).
    - `obstacles`: triangles representing higher-level “plates” or aggregate surfaces; very few globally → small `nOccMax` suffices.
    - `Rrange = {R_full, R_safe}` thresholds for distance-to-ray tests; `R_full < R_safe`.
  - Output aids hierarchy: `distMin[(i,j)]` (closest obstacle distance to ray), compact `occluders` list for the next pass.
  - Note: alternative broad-phase shapes like disks defined by point(x,y,z), normal(u) and radius r could be used later; even elliptical disks if u encodes orientation (e.g., quaternion). Not implemented now.

- __Narrow-phase (fine, per-element)__
  - For pairs (i,j) marked partial by broad-phase, test exact occlusion using only the triangles grouped under listed occluder IDs.
  - Can compute binary blocking or fractional coverage using the element’s sub-division triangles from rasterization.

Consolidation: the Python prototype `TriangleOcclusionRaytracer.py` now builds from the shared `radiosity.cl` so we can drop the separate `python/Radiosity/kernels/triangle_occlusion.cl` after validation.

### Narrow-phase kernel design (per macro pair)

- __Goal__: For a macro pair (i,j), determine fractional occlusion over all sub-element pairs (ik in i, jk in j) using only listed occluders from broad-phase.
- __Assumptions__:
  - `MAX_SUB_I`, `MAX_SUB_J` ≤ 32 sub-elements per macro element can fit in local memory.
  - `nOccMax` occluders per (i,j) from broad-phase; cannot load all their triangles at once.
  - One work-group processes one (i,j) pair; threads stride over `ni*nj` sub-pairs.
- __Local memory__:
  - Cache sub-points of i and j: `LI[MAX_SUB_I]`, `LJ[MAX_SUB_J]`.
  - Tile occluder triangles to local: `TRILOC[3*TILE_TRIS]`.
- __Parallelization__:
  - Map thread-local index to sub-pair: `kk -> ik=kk/nj, jk=kk%nj`.
  - For each occluder group, cooperatively load triangles by tiles to `TRILOC`, then test hits for current sub-pair.
- __Result__:
  - Accumulate `blocked` hits across threads and write fraction `occPairs[pairId] = blocked/(ni*nj)`.

#### Pseudo-code
```c
// one work-group per (i,j)
load LI[0..ni), LJ[0..nj) to local
blocked = 0
for kk in thread-stride over 0..ni*nj:
  (ik,jk) = divmod(kk,nj)
  P1=LI[ik], P2=LJ[jk]; setup ray basis
  hit=0
  for each occluder group o in listed occluders:
    for base in range(0, nverts(o), 3*TILE_TRIS):
      cooperative load TRILOC[base:base+3*TILE_TRIS]
      for v in 0..nverts_tile step 3:
        if ray intersects tri(TRILOC[v:v+3]) then hit=1; break
      if hit break
  blocked += hit
reduce blocked over work-group, write frac to occPairs[pairId]
```
