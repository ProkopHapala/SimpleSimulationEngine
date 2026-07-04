---
type: TopicalAudit
title: Solid Modeling & CSG
tags: [topic, cpp, glsl, javascript, sdf, ray-marching, csg, primitives, mesh-builder, selection]
---

## Summary

Solid modeling via Signed Distance Functions (SDF) and Constructive Solid Geometry (CSG). C++ SDF primitives (`SDfuncs.h`) for mesh selection and editing. GLSL SDF library (IQ-inspired) for ray-marching rendering with CSG operations (union, subtract, intersect, blend). JavaScript shader generator for WebGL SDF ray-marching with scene tree. Mesh builder with SDF-based vertex/edge selection.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common/geometry/SDfuncs.h` | active | SDF functors: `SDF_point2` (distance squared), `SDF_Sphere` (`|p-c|-r`), `SDF_AABB` (box, signed), `SDF_Cylinder` (capped/uncapped, from two endpoints). `operator()` interface for generic use. |
| C++ | `cpp/common/geometry/MeshBuilder2.h` | active | `selectVertsBySDF(sdf, threshold)`, `selectEdgesBySDF(sdf, threshold)` — SDF-based mesh selection. Uses C++20 ranges. |
| C++ | `cpp/common/geometry/Selection.h` | active | `Selection` class: add/remove/toggle/contains, `selectByPredicate()` template with ranges. `SelectionBanks` for multiple selections. |
| GLSL | `cpp/common_resources/shaders/ShaderToy/SDF_RayMarchingPrimitivesCommented.glslf` | active | IQ's SDF primitives with extensive comments: `sdPlane`, `sdSphere`, `sdBox`, `sdEllipsoid`, `udRoundBox`, `sdTorus`, `sdHexPrism`, `sdCapsule`, `sdTriPrism`, `sdCylinder`, `sdCone`, `sdConeSection`, `sdPryamid4`. Non-Euclidean variants (`length6`, `length8`). CSG ops: `opS` (subtract), `opU` (union, with material), `opI` (intersect), `opBlend` (smooth min). Domain ops: repetition, displacement. |
| GLSL | `js/GLSL_solid_modeling/Primitives.glslf` | active | SDF primitive functions for WebGL ray-marching |
| JS | `js/GLSL_solid_modeling/GLSLscreen.js` | active | WebGL screen-space ray-marching renderer |
| JS | `js/GLSL_solid_modeling/ListOfPrimitives.html` | active | Interactive primitive list demo |
| Doc | `docs/SolidModeling/SolidModeling_js.md` | doc | Design document: JS SDF scene tree, `SdfNode` base class with `compileSDF()` and `compileTrace()` modes, `ShaderGenerator` class, CSG limitations in analytical raytracing |

## Sub-topics

### C++ SDF Functors

`SDfuncs.h` provides `operator()` interface:
- `SDF_Sphere`: `return (p - center).norm() - radius;`
- `SDF_AABB`: `q = abs(p-center) - halfSpan; return max(q.maxComponent(), 0) + length(max(q, 0))`
- `SDF_Cylinder`: from two endpoints + radius. Capped: distance to cap circle. Uncapped: distance to axis line segment.
- Used by `MeshBuilder2::selectVertsBySDF()` and `selectEdgesBySDF()`

### GLSL SDF Primitives (Ray-Marching)

Based on Inigo Quilez's distance functions:
- **Basic**: sphere, box, ellipsoid, plane, torus, capsule, cylinder, cone, hex prism, tri prism, pyramid
- **Non-Euclidean**: `length6`, `length8` for squarish cross-sections (`sdTorus82`, `sdTorus88`, `sdCylinder6`)
- **CSG operations**:
  - `opU(d1, d2)`: union — `min(d1.x, d2.x)` with material code in `.y`
  - `opS(d1, d2)`: subtraction — `max(-d2, d1)`
  - `opI(d1, d2)`: intersection — `max(d1, d2)`
  - `opBlend(d1, d2)`: smooth union — polynomial `smin()` with blend radius `k`
- **Domain operations**: repetition (`mod(p, period)`), displacement

### JavaScript SDF Scene Tree

From `SolidModeling_js.md`:
- `SdfNode` base class with `compile(pName)` → GLSL expression string
- `Sphere extends SdfNode`: `compileSDF(p)` → `length(p) - radius`, `compileTrace(ro, rd)` → `iSphere(ro, rd, radius)`
- `ShaderGenerator`: traverses tree, collects GLSL functions, builds `map()` function, adds ray-march loop + lighting
- **Analytical CSG limitation**: union is trivial in SDF (`min`), but extremely hard in analytical raytracing (requires interval logic)

### SDF-Based Mesh Selection

`MeshBuilder2`:
- `selectVertsBySDF(sdf, threshold)`: selects vertices where `sdf(pos) < threshold`
- `selectEdgesBySDF(sdf, threshold)`: selects edges where both endpoints satisfy `sdf < threshold`
- Uses C++20 ranges and `std::function<double(const Vec3d&)>`

## Parity Status

- **C++ `SDfuncs.h` ↔ GLSL SDF primitives**: C++ has sphere, AABB, cylinder only. GLSL has 15+ primitives. No formal parity — different use cases (selection vs rendering).
- **GLSL `SDF_RayMarchingPrimitivesCommented.glslf` ↔ JS `Primitives.glslf`**: Both based on IQ's functions. Coverage differences unknown.
- **C++ SDF ↔ JS `SolidModeling_js.md` design**: C++ uses functor objects; JS uses class hierarchy with `compile()`. Different paradigms.

## Open Issues

- C++ SDF library has only 3 primitives (sphere, AABB, cylinder) — missing cone, torus, capsule, etc.
- No CSG operations in C++ SDF (no union, subtract, intersect for functors)
- `SDF_Cylinder` capped mode: distance to cap is approximate (uses endpoint distance, not true capped cylinder SDF)
- `SolidModeling_js.md` is design document — implementation status of JS scene tree unclear
- No SDF-to-mesh conversion (marching cubes or similar)
- `MeshBuilder2::selectEdgesBySDF` requires both endpoints inside — no partial edge selection
- `Selection::remove()` marks as `-1` in vec but doesn't compact — may cause issues
- No GPU-side SDF evaluation in C++ (only GLSL fragment shader)
