---
type: TopicalAudit
title: Collision Detection
tags: [topic, cpp, java, collision, broad-phase, sweep-prune, aabb-tree, hash-grid]
---

## Summary

Broad-phase collision detection using three strategies: (1) spatial hash grid (`BroadSpaceMapHash`, `HashMap2D`), (2) sweep-and-prune (`SAPbuff`, `SweepAndPrune`), (3) K-d tree / AABB tree (`KBoxes`, `AABBTree3D`, `AABBTreeBin`). Java has cell-sort-based collision grids. Narrow-phase is delegated to object virtual methods. Documentation covers design trade-offs for each strategy.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/engine/broadPhaseCollision.h` | active | `BroadSpaceMapHash` — uniform grid hash using `std::unordered_multimap<int,int>` + `CubeGridRuler`. Insert/overlap for sphere, box, line, triangle. Paint-based deduplication for multi-cell objects. `sweepAndPrune_simple()` — single-axis SAP. `SAPbuff` — multi-direction SAP with `std::sort`. `SweepAndPrune` — `std::map`-based SAP with `lower_bound`/`upper_bound`. `CrossSAP` — templated cross-set SAP. |
| C++ | `cpp/common/maps/kBoxes.h` | active | `KBoxes` — K-d tree spatial partitioning with branch dissolve/regenerate. `collideSelf()` — intra-tree collision pair generation. Sweep-and-prune integration within branches. `collideBranches()` / `collideBranchesNew()` — inter-branch collision with overlap selection. |
| C++ | `cpp/common/maps/AABBTree3D.h` | active | `AABBNode3D` — templated node with fixed-size leaf arrays. `AABBNode6` — 6-branch node layout. `AABBTreeBin` — Box2D/Bullet-style dynamic AABB tree with `AABBNodeBin` (parent/child1/child2/height). `evalInsertCost()` — SAH (surface area heuristic) cost. Incomplete: `insertLeaf`, `deleteNode` not fully implemented. |
| C++ | `cpp/common/maps/HashMap2D.h` | active | `HashMap2D` — 2D uniform grid hash map. `UHALF` (unsigned short) cell indices with `MAP_OFFSET` for negative coords. `getObjectsInRect()`, `getBucketsInRect()`. Used by NBody collision demo. |
| C++ | `cpp/common/dataStructures/HashMat.h` | active | Hash matrix data structure (see spatial-hashing.md) |
| C++ | `cpp/sketches_SDL/2D/test_NBodyColHashMap.cpp` | active | N-body collision simulation using `HashMap2D`. `assembleForces()` — onside (within bucket) + offside (8 neighbor buckets) pairwise interaction. Boundary reflection. |
| Java | `java/Common/CellSort.java` | active | Cell-sort collision grid (Java) |
| Java | `java/Common/CellSort2D.java` | active | 2D cell-sort collision grid |
| Java | `java/Common/CellSort2D_T.java` | active | Templated 2D cell-sort |
| Java | `java/Tests/Test_CollisionGrids.java` | active | Test harness for collision grid variants |
| Doc | `cpp/common/engine/Notes/BroadPhaseCollision.md` | doc | Design notes: shared cell2object array, multi-cell object handling, projectile-to-object SAP, direction-separated SAP, AABB tree references |
| Doc | `cpp/common/engine/Notes/AABBTree.md` | doc | AABB tree design: leaf/branch distinction, array storage, cache considerations, balancing by tree rotation |
| Doc | `cpp/common/engine/Notes/BoxSweepAndPrune.md` | doc | Box-sweep-prune with const-length/const-number boxes to reduce SAP O(m²) worst case |

## Sub-topics

### Spatial Hash Grid

- `BroadSpaceMapHash`: `CubeGridRuler` computes cell index from position. `std::unordered_multimap<int,int>` maps cell→object. Paint-based deduplication (`isNew[]` flag) for objects spanning multiple cells.
- `HashMap2D`: 2D version with `UHALF` (16-bit) indices, `MAP_OFFSET=0x7FFF` for negative coordinate support. Bucket key: `(iy << 16) + ix`.
- NBody demo: 9-bucket stencil (self + 8 neighbors) for pairwise force assembly.

### Sweep and Prune (SAP)

Three variants:
1. `sweepAndPrune_simple()`: Sort by axis, scan forward, break when `bj.a.x > bi.b.x`. Check remaining axes.
2. `SAPbuff`: Multi-direction SAP. Each direction has sorted `SAPitem` list. `addSphere()` projects onto each direction. `collideSelfObjects()` scans within one direction. `collideCrossObjects()` for two-set collision.
3. `SweepAndPrune`: `std::map<double, AOO>` ordered by projection. `lower_bound`/`upper_bound` range query with `maxSize` padding.

Direction-separated SAP for projectiles: group by minimal velocity projection direction, sort each group, sweep. Amortize sorting via insertion sort (temporal coherence).

### K-d Tree (KBoxes)

- Branches with `KBox` (span + range `{i0, n}` into permuted body array)
- `update()`: dissolve high-cost branches, pick random pivots, reinsert free bodies
- `collideSelf()`: intra-branch + inter-branch collision pair generation
- Optional sweep-and-prune within branches (`bSweep` flag)

### Dynamic AABB Tree (AABBTreeBin)

- Box2D/Bullet-style: `AABBNodeBin` with parent/child1/child2/height
- Surface area heuristic for insertion cost
- AVL-style balancing by tree rotation
- **Incomplete**: `insertLeaf()`, `deleteNode()` not fully implemented

## Parity Status

- **C++ `HashMap2D` ↔ Java `CellSort2D`**: Both implement uniform grid collision. Different approaches (hash map vs cell sort). No formal parity test.
- **C++ `SAPbuff` ↔ Java collision grids**: No direct correspondence. No parity test.
- **`KBoxes` ↔ `AABBTreeBin`**: Both are tree-based spatial structures but different algorithms. No comparison.

## Open Issues

- `AABBTreeBin` incomplete — `insertLeaf()`, `deleteNode()` not implemented, `insertSubs()` commented out
- `AABBNode3D::insert()` always returns `false` — not functional
- `BroadSpaceMapHash` uses `std::unordered_multimap` — performance may be poor vs custom hash map
- `KBoxes::update()` uses `new[]`/`delete[]` for `freeBodies` — should use preallocated buffer
- `SAPbuff::sort()` uses `std::sort` (quicksort) — should use `std::stable_sort` for temporal coherence
- No narrow-phase collision detection (shape-shape intersection) — only broad-phase pair generation
- No GPU/OpenCL collision detection
- `collideBranchesNew()` uses VLAs (`int isel[Bi.n]`) — not portable, stack overflow risk
- Java collision grid tests not audited in detail
