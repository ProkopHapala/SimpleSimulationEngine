---
type: TopicalAudit
title: Spatial Hashing
tags: [topic, cpp, java, hash-map, spatial-hash, collision, neighbor-search, open-addressing, paint-dedup]
---

## Summary

Spatial hashing for broad-phase collision detection and neighbor search. `HashMap64` (open-addressing 64-bit key hash map), `HashMap3D` (unbounded 3D spatial hash with `CubeGridRulerUnbound`), `HashMap2D` (2D bucket-based hash), Java `CellSort` variants. Paint-based deduplication for multi-cell object queries. Used by collision detection, N-body simulation, and molecular dynamics.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/maps/HashMap64.h` | active | Open-addressing hash map with 64-bit keys. `insert(key, val)`, `getAllInBoxOnce(key, out, isOld)` — paint-based dedup. `store[]` array for values. |
| C++ | `cpp/common/maps/HashMapT64.h` | active | Templated 64-bit hash map for very big scenes |
| C++ | `cpp/common/maps/HashMap3D.h` | active | `HashMap3D(HashMap64)`: `CubeGridRulerUnbound` + `HashMap64`. `insert(o, Box)`, `collide(o, Box, boxes, colPairs)`. See `grids-rulers.md`. |
| C++ | `cpp/common/maps/HashMap2D.h` | active | 2D hash map: bucket indexing, insert/remove, query in rectangular region. See `collision-detection.md` and `grids-rulers.md`. |
| C++ | `cpp/common/engine/broadPhaseCollision.h` | active | `BroadSpaceMapHash`: `CubeGridRuler` + `std::unordered_multimap`. `findInCells_paint()` dedup. See `collision-detection.md`. |
| Java | `java/Common/CellSort.java` | active | 1D cell-sort grid for collision detection |
| Java | `java/Common/CellSort2D.java` | active | 2D cell-sort grid |
| Java | `java/Common/CellSort2D_T.java` | active | Templated 2D cell-sort grid |

## Sub-topics

### HashMap64 (Open Addressing)

- 64-bit key → hash → open addressing probe sequence
- `insert(key, val)`: stores value in `store[]` at probed slot
- `getAllInBoxOnce(key, out, isOld)`: finds all values for key, uses `isOld[]` paint array to avoid duplicates
- Paint reset: `isOld[o] = false` after each object is processed
- Collision resolution: linear probing (assumed — not fully audited)

### HashMap3D (3D Spatial Hash)

- Combines `CubeGridRulerUnbound` (position → 64-bit cell key) with `HashMap64` (key → objects)
- `insert(o, Box)`: compute BBox-overlapping cells → insert object into each cell
- `collide(o, Box, boxes, colPairs)`:
  1. Find overlapping cells via `ruler.overlap_BBox()`
  2. Gather objects from cells via `hmap.getAllInBoxOnce()` (paint dedup)
  3. Test BBox overlap for each candidate
  4. Output collision pairs
- Limits: `nMaxCell=256`, `nMaxFound=4096` (exits on overflow)

### HashMap2D (2D Spatial Hash)

- Bucket-based: 2D grid of buckets, each containing list of objects
- `insert(o, x, y)`: hash (x,y) → bucket → append
- `remove(o, x, y)`: find and remove from bucket
- `queryRect(x0, y0, x1, y1, out)`: iterate buckets in rectangle, collect objects
- Used by `test_NBodyColHashMap.cpp` for 2D N-body collision

### Paint-Based Deduplication

Problem: same object may appear in multiple cells (e.g., large BBox overlaps several cells).
Solution: `isOld[]` boolean array:
- Before query: all `isOld[i] = true` (or false, depending on convention)
- When object found: check `isOld[o]` — if already found, skip; else mark and add
- After query: reset marks for found objects
- Avoids `std::unordered_set` overhead — O(1) with array

### BroadSpaceMapHash (std::unordered_multimap)

- Uses `std::unordered_multimap<int,int>` instead of custom `HashMap64`
- `findInCells_paint()`: same paint-based dedup as `HashMap3D`
- `insert()` for point, sphere, box, line, triangle shapes
- `overlap()` for sphere, box, line, triangle queries
- Bounded grid (`CubeGridRuler` with `int` keys)

## Parity Status

- **`HashMap64` ↔ `std::unordered_multimap`**: Custom open-addressing vs STL. `HashMap64` has paint dedup built in; `BroadSpaceMapHash` adds it externally. No formal benchmark.
- **`HashMap3D` ↔ `HashMap2D`**: 3D unbounded vs 2D bounded. Different APIs. No parity test.
- **C++ `HashMap2D` ↔ Java `CellSort2D`**: Both 2D spatial partitioning. C++ uses hash buckets, Java uses cell-sort. No formal parity.

## Open Issues

- `HashMap64` collision resolution not fully audited — may be linear probing
- `HashMap3D` exits on overflow (`exit(0)`) — should resize or return error
- `HashMap3D::collide()` has `bDEBUG` exits — should be conditional compilation
- No resize/rehash support in `HashMap64` — fixed capacity
- `BroadSpaceMapHash` uses VLA `int inds[nIndTmpMax]` on stack — overflow risk
- No GPU spatial hash implementation
- Java `CellSort` variants not fully audited
- `HashMap2D` bucket count fixed — no dynamic resize
- Paint array (`isOld[]`) requires `setBodyBuff()` call — easy to forget
