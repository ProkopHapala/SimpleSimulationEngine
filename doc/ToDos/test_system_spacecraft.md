# Spacecraft Simulator — Concrete Test System Design

**Derived from**: `test_system_for_engine.md` (general framework) + `SpacecraftSimulator.todo.md` (specific tasks)
**Goal**: Map each spacecraft simulator TODO item to concrete test strategies using the 4-layer interaction model.

---

## Core Design Principles

1. **Lua-driven scenarios** — tests use the existing Lua scripting API (`Node`, `Girder`, `Ring`, `Slider`, `Thruster`, etc.) to build spacecraft shapes. This tests the real construction pipeline end-to-end, not abstract arrays.
2. **Topology correctness first** — the key question is: did the mesh builder produce the right connectivity? Tests must dump the edge list / connectivity graph and judge it against expectations.
3. **NaN initialization** — fill all padding/ghost/unused slots with `NaN`. After any topology operation, assert `isfinite` on all active slots. Any accidental read of an empty slot poisons downstream math and makes the bug visible immediately (see `agentic_debugging_principles` §1.1).
4. **Robustness testing** — feed invalid inputs (duplicate points, out-of-range indices, zero-length girders, degenerate triangles) and assert the builder either handles them gracefully or fails loud (never silent corruption).
5. **Visual artifacts for human review** — generate `.svg` (2D wireframe) and `.obj` (3D mesh) for every test case so the user can quickly verify shape correctness without running a GUI.

---

## Test Infrastructure

### Lua Scenario Runner

Tests are Lua scripts that build spacecraft via the existing API. The test harness loads the script, runs `BuildCraft_blocks` (mesh generation), and inspects the result.

```lua
-- tests/lua/scenarios/test_girder_simple.lua
require( "data/lua/utils" )
n1 = Node( origin )
n2 = Node( { 10.0, 0.0, 0.0 } )
g1 = Girder( n1, n2, xvec, 5, 2, { 1.0, 1.0 }, "Steel" )
```

```cpp
// Test harness loads Lua script, runs mesh generation, extracts topology
struct MeshTestResult {
    // Raw topology
    int nVerts, nEdges, nTris, nChunks;
    std::vector<Vec3d> vertPositions;
    std::vector<Vec2i> edgeList;      // (a, b) pairs
    std::vector<Vec3i> triList;       // (a, b, c) triples

    // Quality metrics
    double bboxMin[3], bboxMax[3];
    double minEdgeLen, maxEdgeLen;
    bool anyNaN;                      // any vert/edge/tri slot is NaN/Inf
    bool anyDuplicateVerts;           // two verts within R_snap distance
    bool anyInvalidIndex;             // edge/tri references out-of-range or deleted vert
    bool manifold;                    // every edge has ≤2 faces
    bool noDegenerate;                // no zero-length edges, no zero-area tris
};

MeshTestResult runLuaMeshTest(const char* luaScript, const char* debugDir);
```

### Topology Verification

After mesh generation, the harness:

1. **Dumps edge list as text** — human-readable list of all edges `(a, b)` with vertex positions, so the LLM or human can judge if the connectivity is correct:
   ```
   Edge 0: (0) -> (1)   [0,0,0] -> [2,0,0]   len=2.0
   Edge 1: (0) -> (2)   [0,0,0] -> [0,2,0]   len=2.0
   Edge 2: (1) -> (3)   [2,0,0] -> [2,2,0]   len=2.0
   ...
   ```

2. **Dumps connectivity graph** — adjacency list per vertex:
   ```
   Vert 0 [0,0,0]: neighbors = {1, 2, 4}
   Vert 1 [2,0,0]: neighbors = {0, 3, 5}
   ...
   ```

3. **NaN check** — all unused/padding slots initialized to NaN before mesh generation. After generation, assert active slots are finite. Unused slots must still be NaN (not silently overwritten).

4. **Duplicate point detection** — any two vertices within `R_snapVert` distance flagged as potential welding failure.

5. **Index validity** — all edge/tri indices must reference valid (non-deleted) vertices.

### Artifact Generation

For every test case, generate:

- **`.obj` file** — 3D wireframe + faces, viewable in MeshViewer, Blender, or browser-based WebGL viewer
- **`.svg` file** — 2D projected wireframe (for simple geometries, faster to review than 3D)
- **`topology.txt`** — edge list + connectivity graph (for LLM/human judgment)
- **`metrics.json`** — vertex/edge/tri counts, bbox, quality metrics (for regression comparison)

All artifacts go to `debug/<test_name>/`.

### SVG Export for Topology Verification

The `writeSVG` / `writeSVGMulti` functions in `cpp/apps/MeshViewer/MeshFileFormats.h` provide headless 2D SVG export from any `CMesh`. This is the fastest way to visually verify mesh topology without launching a GUI.

**Key features:**
- **Tight bounding box** — canvas shrinks to geometry extent, no wasted space
- **Configurable projection** — drop X, Y, or Z axis (orthographic)
- **Faces** with opacity (0–1) and optional backface culling
- **Edges** — from explicit edge list, or derived from triangles if `nedge == 0`
- **Points** (small circles, off by default)
- **Vertex/edge numbers** — `<text>` labels for topology debugging (off by default, continuous numbering across multi-mesh)
- **Face normals** — lines from face centroids (off by default)
- **Multi-mesh** — `writeSVGMulti()` composes multiple `CMesh` objects into one SVG with combined bbox and continuous numbering

**Usage in test harness:**
```cpp
#include "MeshFileFormats.h"

// After mesh generation, export topology debug view
SVGExportOpts opts;
opts.showFaces = false;       // wireframe only for topology
opts.showEdges = true;
opts.showVertNums = true;     // label vertex indices
opts.showEdgeNums = true;     // label edge indices
opts.projAxis = 2;            // XY plane (drop Z)
writeSVG("debug/test_girder/topology.svg", *mesh, opts);

// Compare two meshes side by side (e.g., C++ vs JS port)
const CMesh* meshes[2] = { meshCPP, meshJS };
SVGExportOpts optsM;
optsM.showFaces = false;
optsM.showEdges = true;
writeSVGMulti("debug/test_parabola/compare.svg", meshes, 2, optsM);

// Face normals for mesh quality check
SVGExportOpts optsN;
optsN.showFaces = true;
optsN.faceOpacity = 0.3;
optsN.showNormals = true;
optsN.normalLength = 0.15;
writeSVG("debug/test_parabola/normals.svg", *mesh, optsN);
```

**Headless batch export:** `test_svg_export <file.obj> [output_dir] [file2.obj]` loads any `.obj` and exports 7 SVG variants (default, edges-only, transparent+bothfaces, opaque+side, all+front, numbers, normals, multi-mesh).

**Interactive viewer:** `meshViewer` (in `cpp/apps/MeshViewer/`) provides:
- 3D rendering with solid/wireframe/points toggles (W/S/P keys)
- Billboarded vertex/edge/face number labels (7/8/9 keys)
- 3D face normals display
- SVG export with live option toggles (X to export, 1–6 for faces/opacity/backface/edges/axis/points, N/M for vert/edge numbers, B for normals)
- Folder browser for loading `.obj`, `.obj+`, `.npz` files

**When to use SVG vs OBJ:**
- **SVG** — fast 2D topology check (edge connectivity, vertex numbering, face normals direction). Best for flat or single-projection geometries. No GUI needed.
- **OBJ** — full 3D inspection in MeshViewer/Blender. Needed for complex 3D shapes where a single projection is ambiguous.

---

## Per-TODO Test Design

### P0: Bond Breaking (Structural Damage)

**Lua scenario**: Build a simple girder between two nodes, apply extreme load.

```lua
-- tests/lua/scenarios/test_bond_break.lua
require( "data/lua/utils" )
n1 = Node( origin )
n2 = Node( { 10.0, 0.0, 0.0 } )
g1 = Girder( n1, n2, xvec, 10, 2, { 1.0, 1.0 }, "Steel" )
-- Set some bonds weak (via material with low Spull)
```

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | Run simulation N steps with extreme force. Assert: (1) bonds with strain > maxStrain are removed from `neighs`/`neighBs`/`neighB2s`, (2) back-neighbor consistency maintained for surviving bonds, (3) no NaN in active slots after break, (4) NaN still present in unused slots (not silently overwritten) | ✅ |
| **L0** | Topology diff: dump edge list before and after break. Assert only expected bonds are removed. | ✅ |
| **L0** | Conservation: before break `|ΔE| < tol`. After break, energy may jump but must be finite. | ✅ |
| **L0** | Robustness: break all bonds simultaneously (extreme load). Assert no crash, no NaN, all neighbor lists empty. | ✅ |
| **L2** | `.obj` snapshot before/after break. `.svg` wireframe showing broken bonds in red. `topology.txt` with edge list diff. | plot |
| **L3** | `./test_BondBreak --review --t=breakStep` — see structure at moment of break | GUI |

### P0: Self-Collision Detection

**Lua scenario**: Two girders connected by one weak bond. Break bond, push fragments together.

```lua
-- tests/lua/scenarios/test_self_collision.lua
require( "data/lua/utils" )
n1 = Node( origin )
n2 = Node( { 5.0, 0.0, 0.0 } )
n3 = Node( { 10.0, 0.0, 0.0 } )
g1 = Girder( n1, n2, xvec, 5, 2, { 1.0, 1.0 }, "WeakMat" )  -- breaks easily
g2 = Girder( n2, n3, xvec, 5, 2, { 1.0, 1.0 }, "Steel" )
-- After break, apply opposing velocities to fragments
```

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | After bond break + collision: assert (1) `pointBBs`/`edgeBBs` detected contact, (2) impulse applied, (3) no interpenetration (min distance > 0 between fragments), (4) momentum conserved `|ΔP| < tol` | ✅ |
| **L0** | Topology: verify collision pairs are in neighbor list with `l0 < 0` sentinel (collision flag). Verify no bonded pairs are also flagged as collision. | ✅ |
| **L0** | Robustness: fragments with identical positions (degenerate). Assert no NaN, no infinite impulse. | ✅ |
| **L2** | `.obj` at collision moment showing contact points. `topology.txt` with collision neighbor list. | plot |
| **L3** | `./test_SelfCollision --review --t=collisionStep` | GUI |

### P1: Projectile Damage Model

**Lua scenario**: Cube truss, fire projectile ray.

```lua
-- tests/lua/scenarios/test_projectile.lua
require( "data/lua/utils" )
n = {}
for i=0,7 do n[i] = Node( { (i&1)*2.0, ((i>>1)&1)*2.0, ((i>>2)&1)*2.0 } ) end
-- build cube edges as girders
g1 = Girder( n[0], n[1], xvec, 3, 2, { 0.5, 0.5 }, "Steel" )
g2 = Girder( n[1], n[3], xvec, 3, 2, { 0.5, 0.5 }, "Steel" )
-- ... etc for all 12 cube edges
```

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | Cast ray from (10,1,1) toward (-1,1,1). Assert: (1) `TriangleRayTracer::getOcclusion()` returns hit, (2) nearest bond to hit point identified, (3) bond breaks, (4) impulse applied to nodes within radius | ✅ |
| **L0** | Projectile cloud (100 rays, spread): assert hit probability ≈ cross-sectional area / total area (within 10% statistical tolerance) | ✅ |
| **L0** | Robustness: ray parallel to truss (no hit). Ray through deleted vertex. Assert no crash, no NaN. | ✅ |
| **L2** | `.obj` with hit points as red spheres, ray paths as lines. `topology.txt` showing which bonds broke. | plot |
| **L3** | `./test_ProjectileDamage --review` — see projectile, hit, damage | GUI |

### P1: Slider Path Initialization

**Lua scenario**: Ring with girders and sliders.

```lua
-- tests/lua/scenarios/test_slider_paths.lua
require( "data/lua/utils" )
n1 = Node( { 0.0, 0.0, 0.0 } )
n2 = Node( { 10.0, 0.0, 0.0 } )
g1 = Girder( n1, n2, xvec, 8, 2, { 1.0, 1.0 }, "Steel" )
r1 = Ring( { 5.0, 0.0, 2.0 }, yvec, xvec, 16, 3.0, { 0.4, 0.4 }, "Steel" )
-- sliders connecting ring to girder
```

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | After `BuildCraft_blocks`: dump each slider's path as vertex index list. Assert: (1) path has expected length (= nseg or nseg×4 depending on method), (2) all path indices reference valid vertices on the rail component, (3) path is ordered (consecutive along rail) | ✅ |
| **L0** | `findNearestPoint()`: for known position on rail, assert returns correct `t` parameter. Test edge cases: position exactly at path vertex, position between vertices, position off-path. | ✅ |
| **L0** | B-spline parity: C++ `Path::getPos(t)` vs JS `bsplineInterpolate()` for same control points. Tolerance `1e-6`. | ✅ |
| **L0** | Robustness: path on deleted component. Path with 0 vertices. Path with 1 vertex. Assert no crash, no NaN. | ✅ |
| **L2** | `.obj` with slider paths as colored line segments (one color per slider). `.svg` 2D projection showing ring + paths. `topology.txt` with path vertex lists. | plot |
| **L3** | `./test_SliderPaths --review` — see ring with sliders, drag to test | GUI |

### P1: Actuator / Control System

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | PID step response: target = 90°. Assert: (1) all sliders reach target within tolerance, (2) no slider exceeds `maxSpeed`/`forceMax`, (3) overshoot < 10%, (4) settling < 100 steps | ✅ |
| **L0** | Input replay: record keyboard → replay → assert same slider positions at each step (deterministic) | ✅ |
| **L0** | Robustness: target beyond physical limit (slider at end of path). Assert clamped, no NaN, no crash. | ✅ |
| **L2** | Plot: slider position vs time + target line. Plot: control force vs time. | plot |
| **L3** | `./test_Actuator --review` — drive with keyboard | GUI |

### P1: Parabolic Mesh Generators (C++ Port from JS)

**Lua scenario** (or direct C++ API call):

```lua
-- tests/lua/scenarios/test_parabola.lua
require( "data/lua/utils" )
-- ParabolaSheet(R0=1, R=3, L=2, nU=16, nV=8)
-- (requires new Lua binding for parabolic generators after port)
```

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | Generate same mesh in C++ and JS. Assert: (1) same vertex count, (2) same edge count, (3) same tri count, (4) vertex positions match within `1e-6` | ✅ |
| **L0** | Topology: dump edge list from both. Assert identical connectivity graph (same edges, same vertex mapping). | ✅ |
| **L0** | Mesh quality: no degenerate tris (zero area), min angle > 15°, edge length ratio < 5, manifold (every edge ≤2 faces) | ✅ |
| **L0** | Robustness: `R0 = 0` (degenerate), `R0 = R` (zero-height sheet), `nU = 1` (minimal). Assert no crash, no NaN, no zero-area tris (or fail loud if degenerate input should be rejected). | ✅ |
| **L2** | `.obj` from C++ and JS side by side. `.svg` wireframe overlay (C++ blue, JS red). `topology.txt` with edge lists. | plot |
| **L3** | `./test_ParabolaGen --review` — parameter sliders | GUI |

### P2: MHD Thruster Integration

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | Simple truss + magnetic nozzle. Assert: (1) MHD force finite and positive (thrust direction), (2) force on correct attachment nodes, (3) truss deforms along thrust axis | ✅ |
| **L0** | Conservation: `|P| = F_thrust × t` (rigid case, no internal loss) | ✅ |
| **L0** | Robustness: zero magnetic field. Assert zero thrust, no NaN. | ✅ |
| **L2** | Plot: thrust vs time. Heatmap: field intensity. `.obj` deformation snapshot. | plot |
| **L3** | `./test_MHDThruster --review` — field lines + structural response | GUI |

### P2: Radiosity Integration

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | 2 parallel plates. C++ vs Python `Radiosity3D.py`. Assert temperatures match within `1e-3` | ✅ |
| **L0** | Energy conservation: total radiated = total absorbed. View factor: `Σ_j F_ij = 1` | ✅ |
| **L0** | Robustness: single plate (no view factor). Degenerate triangle (zero area). Assert no NaN. | ✅ |
| **L2** | Heatmap on surface. C++ vs Python comparison plot. | plot |
| **L3** | `./test_Radiosity --review` — heat distribution on spacecraft | GUI |

### P2: Telescopic Damper Collision Avoidance

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | Compress damper to 50%. Assert: (1) no inner tube vertex penetrates outer tube wall, (2) collision response prevents interpenetration, (3) force finite | ✅ |
| **L0** | Compress to 90% (extreme). Assert: collision stops further compression, no NaN, no bond break. | ✅ |
| **L0** | Robustness: fully compressed (0% extension). Assert no crash, no NaN. | ✅ |
| **L2** | `.obj` at max compression. Plot: compression distance vs time. | plot |
| **L3** | `./test_TelescopicDamper --review --t=maxCompression` | GUI |

### P2: Port Python OpenCL Truss Solvers to C++

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | Same truss problem, same solver (Jacobi, GS, VBD). C++ vs Python. Assert positions match within `rtol=1e-5` (float32) / `rtol=1e-8` (float64) | ✅ |
| **L0** | Convergence: same iteration count to reach tolerance. Graph coloring: same color count. | ✅ |
| **L0** | Robustness: single-vertex truss (no neighbors). Disconnected truss (two islands). Assert no NaN, no crash. | ✅ |
| **L2** | Plot: convergence curve C++ vs Python, all solvers. | plot |

### P3: SpaceCraft2Truss Integration Gaps

**Lua scenario**: Build spacecraft with all component types.

```lua
-- tests/lua/scenarios/test_full_spacecraft.lua
require( "data/lua/utils" )
n1 = Node( origin )
n2 = Node( { 10.0, 0.0, 0.0 } )
n3 = Node( { 5.0, 10.0, 0.0 } )
g1 = Girder( n1, n2, xvec, 5, 2, { 1.0, 1.0 }, "Steel" )
g2 = Girder( n2, n3, xvec, 5, 2, { 1.0, 1.0 }, "Steel" )
Shield( g1, 0.0, 1.0, g2, 0.0, 1.0 )
Radiator( g1, 0.0, 0.5, g2, 0.0, 0.5 )
Thruster( { 0.0, 0.0, 0.0 }, yvec, { 5.0, 300.0, 20.0 }, "ICF_Ebeam_magNozzle" )
tanks( 2, 0.0, 5.0, 2.0, 50.0, 100.0 )
```

| Layer | Test | Auto? |
|-------|------|-------|
| **L0** | After `toTruss()`: dump edge list + connectivity graph. Assert: (1) every component has truss nodes, (2) shields → triangulated faces, (3) radiators → faces, (4) thruster → node with force, (5) tanks → nodes with mass. Verify by checking topology, not just counts. | ✅ |
| **L0** | No orphan components: every `ShipComponent*` has truss representation. Check by walking `SpaceCraft` vectors and matching to truss nodes/edges/faces. | ✅ |
| **L0** | Robustness: component with no material. Component with zero-size. Assert fail-loud (crash with message) or graceful skip with warning. | ✅ |
| **L2** | `.obj` with components color-coded by type. `topology.txt` with full edge list + component mapping. | plot |
| **L3** | `./test_SpaceCraft2Truss --review` — see full spacecraft in truss | GUI |

---

## Test Execution Matrix

| Test | L0 Auto | L1 LLM | L2 Artifacts | L3 GUI | Ref Data |
|------|---------|--------|--------------|--------|----------|
| Bond Breaking | ✅ | — | .obj, .svg, topology.txt | ✅ | ✅ |
| Self-Collision | ✅ | — | .obj, topology.txt | ✅ | ✅ |
| Projectile Damage | ✅ | hit pattern | .obj, topology.txt | ✅ | ✅ |
| Slider Paths | ✅ | — | .obj, .svg, topology.txt | ✅ | ✅ |
| Actuator Control | ✅ | — | response plot | ✅ | ✅ |
| Parabolic Mesh | ✅ | — | .obj, .svg, topology.txt | ✅ | ✅ |
| MHD Thruster | ✅ | field sanity | thrust plot, heatmap | ✅ | — |
| Radiosity | ✅ | — | heatmap, comparison | ✅ | ✅ |
| Damper Collision | ✅ | — | .obj, compression plot | ✅ | — |
| Solver Parity | ✅ | — | convergence plot | — | ✅ |
| Component→Truss | ✅ | completeness | color .obj, topology.txt | ✅ | — |

## Priority Implementation Order

1. **Lua scenario runner** — loads Lua script, runs `BuildCraft_blocks`, extracts topology (edge list, connectivity graph, positions)
2. **Topology dump + NaN check** — `topology.txt`, `metrics.json`, NaN initialization and verification
3. **Artifact generators** — `.obj` exporter (3D mesh), `.svg` exporter (2D wireframe projection)
4. **First test: girder simple** — validate the runner with a basic girder, judge edge list correctness
5. **Bond breaking test** (P0) — first dynamics test
6. **Self-collision test** (P0) — builds on bond breaking
7. **Slider path test** (P1) — validates path initialization
8. **Parabolic mesh parity test** (P1) — validates JS→C++ port
9. **`run_all_tests.sh`** — aggregate runner once 3+ tests exist
10. **GUI review mode** — add `--review` to test apps
