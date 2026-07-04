# geometry

Geometric primitives, mesh data structures, and procedural mesh generation. The core of the engine's 3D modeling and simulation pipeline.

- **raytrace.h** — Ray-shape intersection primitives (ray-sphere, ray-plane, ray-triangle Möller–Trumbore, ray-AABB slab method, ray-cylinder/cone). Foundation for mouse picking, raytracing, and collision queries. Uses t_inf=1e+300 as "no hit" sentinel.
- **geom2D.h** — 2D computational geometry: point-in-shape tests, line/segment intersection, triangle area, convex hull. All templated on Vec2T<T> for float/double/int. Companion to geom3D.h.
- **geom3D.h** — 3D computational geometry: distance queries (point-to-triangle with edge clipping), intersection tests (ray-triangle, sphere-sphere, AABB-AABB), analytic geometry (circle from 3 points, plane from 3 points). Companion to geom2D.h.
- **Solids.h** — Static compile-time CMesh definitions of platonic solids (tetrahedron, octahedron, cube, icosahedron, dodecahedron). Zero allocation — CMesh points to static arrays. Some solids include plane cross-sections for ConstructionBlock insertion points.
- **CMesh.h** — Non-owning C-style mesh view (raw pointers, no methods). Bridge between Builder2's std::vector storage and C/OpenCL APIs. Used for GPU upload and file I/O.
- **SDfuncs.h** — Signed distance functions as functor objects (SDF_Sphere, SDF_AABB, SDF_Cylinder). Composable — pass any SDF to selection code. Functor pattern enables inlining.
- **UVfuncs.h** — Parametric surface generators mapping (u,v)→Vec3f (sphere, cone, cylinder, paraboloid, bicubic patch). Normal computation via finite-difference cross product. For procedural mesh creation.
- **MeshBuilder.h** — Lightweight render-ready mesh builder (Mesh::Builder). Stores Vec3f positions/normals/UVs + indices in flat vectors. Rendering abstraction without rendering dependencies — same code produces OpenGL VBO, OpenCL buffer, or .obj. Predecessor to MeshBuilder2.
- **MeshBuilder2.h** — The main mesh editing system (Builder2). Flat arrays of VertT (64-byte union: static pos/nor/uv or dynamic flags/uid/neighbors), Quat4i edges/tris/chunks, strips for polygon storage. Soft-remove with alive flags + compactDead() GC. Topology ops: bevel, extrude, bridge, edge loop sorting. Inherits SelectionBanks for multi-selection editing.
- **MeshBuilder2.cpp** — Implementation of Builder2: compactDead(), neighbor rebuild, OBJ I/O, truss export, all topology operations, selection by SDF/ray/box/sphere.
- **Mesh.h** — Deprecated OMesh (half-edge mesh with std::vector + soft-removal). Contains OBJ import/export, MeshEdge with crease angle, Polygon with triangulation cache. Predecessor to MeshBuilder2 — kept for OBJ parser compatibility.
- **Mesh.cpp** — Implementation file for Mesh.h (includes arrayAlgs.h and raytrace.h).
- **Selection.h** — Dual-storage selection (ordered vector + hash map) for interactive editing. remove() marks -1 (preserves indices). SelectionBanks manages multiple named selections simultaneously.
- **ConstructionBlock.h** — ZomeTool-like construction set: build trusses from platonic solid connectors. cubeIndex3D() packs 4D cube element IDs into 8 bits. Tracks beam connectivity for truss simulation.
- **IntersectionCurve.h** — Traces 1D curve from intersection of N implicit surfaces F_i(p)=C_i. Gradient descent relaxation + curve-direction stepping with curvature-adaptive step size. Same principle as NEB/string methods in molecular dynamics.
- **Voronoi.h** — Fortune's algorithm for 2D Voronoi diagrams (beach-line BST + priority queue). O(N log N). Based on Ivan Kuckovic's implementation. For terrain partitioning, procedural layout, nearest-neighbor.
- **Voronoi.cpp** — Implementation of Voronoi.h.
- **voronoi.h** — Alternative Voronoi implementation (Steven Fortune's original C structure). Uses pool allocator (Freelist/Freenode) and half-edge data structure. Coexists with Voronoi.h for different tradeoffs.
- **voronoi.cpp** — Implementation of voronoi.h.
- **Convex2d.cpp** — Implementation of Convex2d (convex polygon cutting by line). Header is in math/Convex2d.h. Used for 2D clipping operations.
