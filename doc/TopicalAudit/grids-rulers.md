---
type: TopicalAudit
title: Grids, Rulers & Spatial Hashing
tags: [topic, cpp, java, grid, ruler, hash-map, spatial-partition, cell-indexing, buckets, neighbor-search]
---

## Summary

Spatial partitioning infrastructure for neighbor search and collision detection. `CubeGridRuler` (bounded) and `CubeGridRulerUnbound` (unbounded, 64-bit interleaved key) for cell indexing. `HashMap3D` combines ruler + `HashMap64` for unbounded spatial hashing. `Buckets3D` extends `Buckets` with `CubeGridRuler` for bounded cell lists with 27-neighbor stencil. `HashMap2D` for 2D spatial hashing. Java cell-sort grids for 2D collision detection.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/maps/grids3D.h` | active | `CubeGridRuler` (bounded): `pos2box()`, `icell()`, `overlap_Sphere()` (2×2×2), `overlap_BBox()`, `ixyz2i()`, `ixyz2i_wrap()` (PBC). `CubeGridRulerUnbound`: 64-bit interleaved key, `pos2long()`, `overlap_BBox()`. `Grids3D::insert_SphereOfInfluence()` template. See `parallel-particle-cell.md`. |
| C++ | `cpp/common/maps/HashMap3D.h` | active | `HashMap3D(HashMap64)`: unbounded spatial hash. `CubeGridRulerUnbound` + `HashMap64`. `insert(o, Box)`, `collide(o, Box, boxes, colPairs)` — BBox overlap + paint-based dedup. `nMaxCell=256`, `nMaxFound=4096`. |
| C++ | `cpp/common/maps/HashMap64.h` | active | `HashMap64`: 64-bit key hash map with open addressing. `insert()`, `getAllInBoxOnce()` (paint-based dedup). |
| C++ | `cpp/common/maps/HashMapT64.h` | active | `HashMapT64`: templated 64-bit hash map for very big scenes |
| C++ | `cpp/common/maps/HashMap2D.h` | active | 2D hash map: bucket indexing, insert/remove/query in rectangular region. See `collision-detection.md`. |
| C++ | `cpp/common/maps/Buckets.h` | active | `Buckets`: cell-based bucket storage. `resizeCells()`, `resizeObjs()`, `updateCells()`, `getInCell()`. `obj2cell` mapping, `cell2obj` linked list. |
| C++ | `cpp/common/maps/Buckets3D.h` | active | `Buckets3D(Buckets, CubeGridRuler)`: bounded 3D cell list. `pointsToCells()` — assign points to cells. `getNeighbors()` — 27-neighbor stencil (3×3×3 minus self). `getForwardNeighbors()` — 13 neighbors (half stencil, avoids double-counting). PBC wrapping via `ixyz2i_wrap()`. |
| C++ | `cpp/common/maps/kBoxes.h` | active | K-d tree spatial partitioning. See `collision-detection.md`. |
| C++ | `cpp/common/maps/AABBTree3D.h` | active | Dynamic AABB tree. See `collision-detection.md`. |
| C++ | `cpp/common/engine/broadPhaseCollision.h` | active | `BroadSpaceMapHash`: `CubeGridRuler` + `std::unordered_multimap`. See `collision-detection.md`. |
| Java | `java/Common/CellSort.java` | active | Cell-sort 1D grid for collision detection |
| Java | `java/Common/CellSort2D.java` | active | Cell-sort 2D grid for collision detection |
| Java | `java/Common/CellSort2D_T.java` | active | Templated cell-sort 2D grid |

## Sub-topics

### CubeGridRuler (Bounded)

- `setup(pmin, pmax, step)` — defines grid extent
- `pos2box(pos, ipos, dpos)` — position → cell index + fractional position
- `icell(pos)` → `ixyz2i(ip)` = `ix + nx*(iy + ny*iz)` (row-major)
- `overlap_Sphere(pos, r, icells)` — up to 8 cells (2×2×2 neighborhood), edge/corner checks via `dr² < r²`
- `overlap_BBox(p0, p1, icells)` — all cells in BBox range
- `ixyz2i_wrap(ip)` — periodic boundary wrapping via `wrap_index_fast()`
- `overlap_Line()` and `overlap_Triangle()` — **not implemented** (return -1)

### CubeGridRulerUnbound

- 64-bit interleaved key: `NBIT3D = 21` bits per axis, `OFFSET3D = 2^20`
- `ixyz2long(ip)` = `(ip.x+OFF) + ((ip.y+OFF) + ((ip.z+OFF)<<NBIT)<<NBIT)`
- Range: ±2^20 ≈ ±10^6 (with `step=1`)
- No boundary checking — unbounded domain
- Only `overlap_BBox()` — no `overlap_Sphere()`

### HashMap3D (Unbounded Spatial Hash)

- `CubeGridRulerUnbound` + `HashMap64` (open addressing)
- `insert(o, Box)`: compute overlapping cells → insert object into each cell
- `collide(o, Box, boxes, colPairs)`: find objects in overlapping cells → paint-based dedup → BBox overlap test → collision pairs
- `isOld[]` array for paint-based dedup (reset after each query)
- Limits: `nMaxCell=256`, `nMaxFound=4096` (exits on overflow)

### Buckets3D (Bounded Cell List)

- `Buckets` + `CubeGridRuler` (bounded)
- `pointsToCells(np, ps)`: assign each point to its cell, update cell lists
- `getNeighbors(ip, neighs)`: 27-neighbor stencil (3×3×3 minus self), PBC wrapping
- `getForwardNeighbors(ip, neighs)`: 13 neighbors (half stencil) — avoids double-counting for symmetric interactions
- `maxInBucket` tracking, `neighs_in` buffer auto-resize (`maxInBucket * 14`)

### HashMap2D

- 2D spatial hash with bucket indexing
- Insert/remove/query in rectangular region
- Used by `test_NBodyColHashMap.cpp` for 2D N-body collision

## Parity Status

- **`CubeGridRuler` ↔ `CubeGridRulerUnbound`**: Same cell indexing concept, different key types (int vs uint64). Bounded has `overlap_Sphere()`, unbounded does not.
- **`HashMap3D` ↔ `BroadSpaceMapHash`**: Both use `CubeGridRuler` + hash map. `HashMap3D` uses custom `HashMap64`, `BroadSpaceMapHash` uses `std::unordered_multimap`. Different dedup strategies (paint vs `findInCells_paint`).
- **`Buckets3D` ↔ `HashMap3D`**: `Buckets3D` is bounded with PBC, `HashMap3D` is unbounded. `Buckets3D` uses cell lists (array-based), `HashMap3D` uses hash map.
- **C++ `Buckets3D` ↔ Java `CellSort2D`**: Both are cell-list approaches. C++ is 3D with PBC, Java is 2D. No formal parity.

## Open Issues

- `CubeGridRuler::overlap_Line()` and `overlap_Triangle()` not implemented (return -1)
- `CubeGridRulerUnbound` lacks `overlap_Sphere()` — only BBox
- `HashMap3D` exits on overflow (`nMaxCell`, `nMaxFound`) — should resize
- `Buckets3D::getNeighbors()` uses 27 hardcoded `getInCell()` calls — could use loop
- `Buckets3D::getForwardNeighbors()` has 13 hardcoded calls — mirror of `getNeighbors()`
- No GPU implementation of spatial hashing — all CPU-side
- `HashMap64` not fully audited — open addressing collision resolution unknown
- `BroadSpaceMapHash` uses VLA `int inds[nIndTmpMax]` on stack — overflow risk
- No 2D equivalent of `Buckets3D` (only `HashMap2D` which is hash-based, not cell-list)
- `CubeGridRuler::ipcell()` has operator precedence bug: `(int)(pos.x-pos0.x)*invStep` — cast applies only to subtraction, not multiplication
