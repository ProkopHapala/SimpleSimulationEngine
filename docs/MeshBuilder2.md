# `Mesh::Builder2` - The Unified Mesh Toolkit

## MeshBuilder2 Documentation

The `MeshBuilder2` class is a versatile tool designed for generating and manipulating 3D meshes within the SimpleSimulationEngine. Unlike typical polygon mesh modelers primarily focused on rendering, `MeshBuilder2` is specifically tailored for physics simulations, emphasizing volumetric 3D meshes (trusses) for volumetric 3D element simulations and truss simulations. It focuses more on points connected by edges (mass points and distance constraints - tubes, sticks) rather than solely on polygons. However, it is now also being adapted for radiosity and light-scattering simulations, which will involve more extensive use of surfaces (polygons).

It consolidates functionality that might otherwise be spread across multiple classes, providing a unified API for:
-   **Low-level primitive creation**: Adding individual vertices, edges, and triangles.
-   **High-level shape generation**: Creating complex objects like boxes, girders, and tubes with single function calls.
-   **Topological editing**: Performing operations like extruding faces, bridging polygons, and creating faces from edge loops.
-   **Interactive selection**: Picking and selecting vertices, edges, and faces for user-driven editing.

## 2. Core Data Structures

The builder stores all geometric and topological data in a set of `std::vector` containers. Understanding these is key to using the class effectively.

-   **`verts`**: `std::vector<Vert>`: The fundamental list of all vertices in the mesh. Each `Vert` object stores a `pos` (position), `nor` (normal), and `uv` (texture coordinate).

-   **`edges`**: `std::vector<Quat4i>`: A list of all edges. Each edge is a `Quat4i` storing `{v1, v2, type, type2}`, where `v1` and `v2` are indices into the `verts` vector.

-   **`tris`**: `std::vector<Quat4i>`: A list of all triangles. Each triangle is a `Quat4i` storing `{v1, v2, v3, face_id}`, where `v1`, `v2`, `v3` are vertex indices. The `face_id` links the triangle back to a logical face (a "chunk").

-  This is the system for defining complex primitives like polygons:
   -  **`chunks`** :`std::vector<Quat4i>`:  A metadata vector. Each `Quat4i` chunk entry stores `{i0, i0+n, n, type}`, where `i0` is the start index in the `strips` vector, `n` is the number of elements, and `type` is an enum (`ChunkType::face`, `ChunkType::edgestrip`, etc.)
   -  **`strips`** : `std::vector<int>`:  A flat vector of integers that stores sequences of vertex or edge indices.
    -   Example: A square face chunk would point to a section in `strips` containing four vertex indices followed by four edge indices. This flexible system allows for polygons of any size.

-   **`vert2edge`** `std::unordered_map<uint64_t, int>`: An optional hash map for accelerating edge lookups. It maps a key derived from two vertex indices to the corresponding edge index. This is much faster than a brute-force search.

### 2.1. The Chunk/Strip System

A key feature of `MeshBuilder2` is its ability to store complex, variable-length primitives like polygons or line strips. This is achieved through a flexible two-part system composed of the `chunks` and `strips` vectors.

*   **`strips`** (`std::vector<int>`): This is a flat, contiguous buffer that stores the raw data for all complex primitives. It is simply a long list of integer indices. These indices can refer to vertices, edges, or other primitives depending on the context defined by the chunk.

*   **`chunks`** (`std::vector<Quat4i>`): This is the metadata or descriptor vector. Each `Quat4i` entry in this vector defines a single "chunk" and describes a specific segment within the `strips` vector. The `Quat4i` for a chunk is interpreted as `{ i0, i_end, n, type }`:
    -   `x` (`i0`): The starting index of the chunk's data within the `strips` vector.
    -   `y` (`i_end`): The end index of the first part of the chunk's data (e.g., the vertex list for a face).
    -   `z` (`n`): The number of elements in the primitive (e.g., the number of vertices in a polygon).
    -   `w` (`type`): The type of the chunk, defined by the `ChunkType` enum (`face`, `edgestrip`, `trianglestrip`).

**Example: A Face Chunk (`ChunkType::face`)**

This is the most common and important use case. When a face chunk is created by `polygonChunk()`, its data is laid out in a specific format within the `strips` vector: `n` vertex indices followed by `n` edge indices.

`strips`: `[..., v0, v1, v2, v3, e0, e1, e2, e3, ...]`

The corresponding `chunk` descriptor would point to this data. For a 4-vertex face (a quad) starting at index `i0` in `strips`, the chunk would be: `{ i0, i0+4, 4, ChunkType::face }`.

**Accessing Chunk Data**

-   `getChunkStrip(ich)`: This is the low-level accessor. It returns a raw `int*` pointer to the beginning of the chunk's data in `strips` (i.e., to `strips[i0]`). The calling code must know the chunk's type to correctly interpret the subsequent data.

-   `loadChunk(ich, iedges, iverts)`: This is a high-level helper function specifically designed for `ChunkType::face`. It understands the `[vertices][edges]` layout and automatically populates the provided `iverts` and `iedges` arrays with the correct indices from the `strips` vector. This is the preferred way to read face data for operations like extrusion.

## 3. Key Concepts and Workflows

The `MeshBuilder2` provides a comprehensive set of tools for modifying mesh topology and geometry. These functions enable operations similar to those found in advanced 3D modeling software, adapted for the specific needs of both rendering and physics simulations with polygon surfaces (e.g. in Radiosity and Light-Scattering simulations) and finite element simulations (e.g. dynamics of trusses, mass-spring systems, etc).

The `MeshBuilder2` provides a complete suite of tools for interactive mesh editing, driven by user input.

1.  **Set Mode**: The user selects a `SelectionMode` (`vert`, `edge`, or `face`) using the `I`, `O`, `P` keys in the demo app.
2.  **Pick**: The user clicks, and `pickSelect()` is called. This uses ray-tracing (`pickVertex`, `pickEdge`, `pickTriangle`) to find the element under the cursor. The picked element's index is stored.
3.  **Select**: The picked index is added to the `selection` vector and `selset` set. The user can select multiple elements if `bAdditiveSelect` is true. A selection box can also be used with `selectRect()`.
4.  **Operate**: The user presses a key to perform an action on the selection.
    -   **`F` key -> `selectionToFace()`**: If a loop of edges is selected, this function sorts them and creates a new polygon face.
    -   **`E` key -> `extrudeFace(ipick, ...)`**: If a face is selected, this function extrudes it. It gets the face's vertices from its chunk, creates new vertices offset along the normal, creates a new face, and calls `bridge_quads` to connect the old and new faces, forming the sides of the extrusion.

### 3.1. Snapping and Finding Vertices

Many functions, especially the `snap...Face` family, rely on finding existing vertices rather than creating new ones.

-   **`findVert(pos, R_snapVert)`**: This method searches for a vertex within a given radius (`R_snapVert`) of a specified position. This is crucial for "welding" new geometry to an existing mesh.
-   **`findOrAddEdges(verts)`**: This method checks if an edge between two vertices already exists (using `findEdgeByVerts`). If not, it creates a new one. This prevents duplicate edges.

### 3.2. Face Connector Generation (`snap...Face`)

These functions are the bridge between the high-level `ConstructionBlock` and the low-level mesh. They are responsible for creating the physical connector geometry on the face of a block.

-   `snapBoxFace(p0, rot, La, Lb)`: Creates a simple, flat quad face, snapping to existing vertices.
-   `frustrumFace(p0, rot, La, Lb, h, Lbh, Lah)`: Creates a frustum-shaped face.
-   `snapFrustrumFace(p0, rot, La, Lb, h, Lbh, Lah, bFace)`: Creates a frustum-shaped face, snapping to existing vertices.
-   `prismFace(p0, rot, La, Lb, h, Lbh)`: Creates a prism-shaped face.
-   `snapPrismFace(p0, rot, La, Lb, h, Lbh, bFace)`: Creates a prism-shaped face, snapping to existing vertices.
-   `quad(q, face_type, edge_type)`: Creates a quadrilateral face from four vertex indices.
-   `bridgeFacingPolygons(p0, p1, ch1, ch2, nseg, stickTypes, maks)`: Bridges two facing polygons between two points.
-   `bridgeFacingPolygons(nrod, edges, points, nseg, chs, stickTypes, maks)`: Bridges multiple facing polygons using edge and point data.
-   `facingNodes(cmesh, nnod, points, out_chs, nplane, planes, planeVs)`: Identifies nodes facing a certain direction within a CMesh. It uses `findVert` to locate the four corner vertices of the face on the existing block's mesh and then uses `findOrAddEdges` and `polygon` to create the face chunk.
-   `snapPrismFace(p0, rot, ...)`**: Creates a prism-shaped connector that provides two new faces for connections. It finds two vertices on the block and creates four new vertices to form the prism. It then generates two quad face chunks.
-   `snapFrustrumFace(p0, rot, ...)`**: Creates a more complex frustum-shaped connector that provides five new faces. It finds two vertices on the block and creates six new vertices to form the frustum, generating five quad face chunks.

The `bFace` parameter in these functions controls whether the connecting quad polygons are actually generated and stored as chunks. This is essential for the `ConstructionBlock` workflow, which needs these chunks to connect girders.

## 4 MeshBuilder2 Functions

### 4.1 Primitive Creation
-   `vert(pos, ...)`: Adds a vertex.
-   `edge(v1, v2, ...)`: Adds an edge between two vertices.
-   `tri(v1, v2, v3, ...)`: Adds a triangle.
-   `addVerts(n, ps)`: Adds multiple vertices from an array of positions.
-   `addEdges(n, iedges, types, types2, iv0)`: Adds multiple edges from arrays of vertex index pairs and types.
-   `addFaces(n, nVerts, iverts, bPolygonToTris, iv0)`: Adds multiple faces (polygons) from arrays of vertex counts and vertex indices.
-   `addCMesh(cmesh, bFaces, p0, sc, rot, edge_type)`: Adds a complete CMesh object to the builder, with optional transformations.
-   `add_box(p0, p1, bEdges, bFaces)`: Adds a rectangular box.
-   `girder1(p0, p1, up, n, width, stickTypes, bCaps)`: Creates a 3D truss-like girder structure between two points.
-   `triangle_strip(p0, p1, up, n, width, stickType, bCaps)`: Generates a flat, zig-zagging triangle strip.
-   `chunk(...)`: Adds a metadata chunk pointing to a strip of primitives.
-   `quad(q, face_type, edge_type)`: Creates a quadrilateral face from four vertex indices.
-   `addPointCross(const Vec3d& p, double d)`: Adds a 3D cross shape at a given position `p` for debugging visualization.
-   `addArrow(const Vec3d& p1, const Vec3d& p2, double d)`: Adds an arrow from `p1` to `p2` for debugging visualization.
-   `addSphere(const Vec3d& p, double r, int n)`: Adds a sphere at a given position `p` with radius `r` and `n` subdivisions.
-   `addCone(const Vec3d& p, double r, double h, int n)`: Adds a cone at a given position `p` with radius `r`, height `h`, and `n` subdivisions.
-   `addCylinder(const Vec3d& p, double r, double h, int n)`: Adds a cylinder at a given position `p` with radius `r`, height `h`, and `n` subdivisions.

### 4.2  Complex Shape Generation
-   `box(p, ls, rot)`: Creates a box with 8 vertices and 12 edges.
-   `ring(...)`: Creates a circular loop of vertices and edges.
-   `tube(...)`: Creates a tube by connecting a series of rings.
-   `plate(...)`: Creates a subdivided flat surface between four corner points.
-   `plate_quad(ip00, ip01, ip10, ip11, typs, n, fillType)`: Creates a subdivided plate from four vertex indices.
-   `plateOnGriders(ns, prange1, prange2, byN, offs, span1, span2, stickTypes)`: Creates a plate structure on top of existing girders.
-   `girder1_caps(ip0, ip1, kind)`: Creates a girder structure between two existing vertex indices, adding caps.
-   `girder1(ip0, ip1, up, n, width, stickTypes)`: Creates a girder structure between two existing vertex indices.
-   `wheel(p0, p1, ax, n, wh, stickTypes)`: Creates a wheel-shaped truss structure between two points.
-   `ngon(p0, p1, ax, n, stickType)`: Creates an N-sided polygon (ngon) structure.
-   `rope(p0, p1, nseg, ropeType, anchorType, Rcolapse, r)`: Creates a rope-like structure between two points with specified segments and types.
-   `ropes(nv, vs, ne, nseg, ends, ropeType, anchorType, Rcolapse, r)`: Creates multiple rope-like structures between arrays of vertices and ends.
-   `vstrip(p0, p1, n, et)`: Creates a vertical strip of edges.
-   `fstrip(ip0, ip1, n, ft, et)`: Creates a strip of faces between two vertex indices.
-   `panel(...)`: Creates a subdivided panel with truss-like internal structure.

### 4.3 Selection and Picking
These functions allow for precise selection of mesh elements (vertices, edges, faces) for subsequent operations.
-   `pickVertex(ray0, hRay, R)`: Finds the closest vertex to a ray.
-   `pickEdge(ro, rh, Rmax)`: Finds the closest edge to a ray.
-   `pickTriangle(ro, rh, bReturnFace)`: Finds the triangle intersected by a ray, with an option to return the face.
-   `pickEdgeSelect(ro, rh, Rmax)`: Selects an edge by ray intersection.
-   `findClosestVert(p0, i0, n)`: Finds the closest vertex to a point within a specified range.
-   `closestInSelection(p0, Rmax, n, sel)`: Finds the closest vertex to a point within the current selection.
-   `pickSelect(ro, rh, Rmax)`: Main picking function that dispatches to the above based on `selection_mode`.
-   `selectRect(p0, p1, rot)`: Selects all elements within a screen-space rectangle.
-   `selectRectEdge(p0, p1, rot)`: Selects edges within a screen-space rectangle.
-   `selectRectVert(p0, p1, rot)`: Selects vertices within a screen-space rectangle.
-   `select_in_sphere(p0, r, i0, imax)`: Selects vertices within a specified sphere.
-   `select_in_cylinder(p0, fw, r, l, i0, imax)`: Selects vertices within a specified cylinder.
-   `select_in_box(p0, fw, up, Lmin, Lmax)`: Selects vertices within a specified bounding box.
-   `clearSelection()`: Clears the current selection.
-   `toggleSelSet(i)`: Toggles the selection status of an element.
-   `makeSelectrionUnique()`: Ensures that selected elements are unique.
-   `selectVertsAlongLine(p0, p1, r, bSort)`: Selects vertices along a line segment.
-   `selectVertsAlongPolyline(r, bSort)`: Selects vertices along a polyline.

### 4.4 Geometric Transformations and Utilities
These functions modify the positions of mesh elements.
-   `move_verts(indices, shift)`: Translates a specified set of vertices by a given shift vector.
-   `scale_verts(indices, p, s)`: Scales a specified set of vertices relative to a pivot point `p` by scale factors `s`.
-   `rotate_verts(indices, p, rot)`: Rotates a specified set of vertices around a pivot point `p` by a given rotation matrix `rot`.
-   `getCOG(n, ivs)`: Calculates the center of gravity for a set of vertices.
-   `getChunkNormal(ich)`: Calculates the normal vector of a specified chunk. (redudant?)
-   `polygonNormal(ich)`: Calculates the normal vector of a polygon chunk.
-   `findVert(p0, Rmax, ...)`: Finds the closest vertex to a point within a given radius.
-   `facingNodes(cmesh, nnod, points, out_chs, nplane, planes, planeVs)`: Identifies nodes facing a certain direction within a CMesh.

### 4.5 Topological Utilities
-   `findEdgeByVerts(verts)`: Finds the index of an edge connecting two vertices. Can use a fast map (`_map`) or slow brute-force (`_brute`).
-   `sortEdgeLoop(n, iedges, iverts)`: Sorts an array of edge indices into a contiguous loop, optionally returning vertex indices.
-   `sortPotentialEdgeLoop(n, edges, iverts)`: Sorts a potential edge loop.
-   `findEdgeByVerts_brute(verts)`: Finds an edge by its two vertex indices using a brute-force search.
-   `findEdgeByVerts_map(verts)`: Finds an edge by its two vertex indices using a map for faster lookup.
-   `findOrAddEdges(verts, t, t2)`: Finds an edge by its vertices, or adds it if it doesn't exist.
-   `buildVerts2Edge()`: Builds `vert2edge` map map for fast lookup of edges connected to vertices.
-   `build_edgesOfVerts()`: Builds adjacency list mapping vertices to connected edges.
-   `loadNeighbours(iv, ivs, ies, n)`: Gets neighboring vertices and edges for vertex `iv`.
-   `vertNormalByEdges(iv, bNormalizeEach)`: Computes vertex normal from connected edges.
-   `sortVertEdgesByNormal(p, nor, n, ies)`: Sorts edges around vertex normal.
-   `alling_polygons(n, ivs1, ivs2, ipiv)`: Aligns two sets of polygons based on their vertex indices.
-   `loadChunk(ich, iedges, iverts)`: Loads data from a specified chunk into edge and vertex arrays.
-   `polygonChunk(n, iedges, ivs, bPolygonToTris)`: Creates a polygon chunk from edge and vertex indices.
-   `polygonToTris(i)`: Converts a polygon chunk to triangles.
-   `findMostFacingNormal(hray, nch, chs, cosMin, bTwoSide)`: Finds the chunk with the normal most facing a given ray from a list of chunks.
-   `findMostFacingNormal(hray, chrange, cosMin, bTwoSide)`: Finds the chunk with the normal most facing a given ray within a chunk range.
-   `clear()`: Clears all mesh data (vertices, edges, triangles, chunks, etc.).

### 4.6 Topological Operations
These functions modify the connectivity and structure of the mesh.
-   `extrudeFace(ich, L, stickTypes, maks)`: Extrudes a face chunk outwards by distance `L`.
-   `extrudeVertLoop(n, iverts, d, bEdges, bFace, bTris, bSort)`: Extrudes a loop of vertices.
-   `bevel_vert(iv, L, h, ies, nor0)`: Bevels a single vertex with given width and height.
-   `bevel(ne, ies, L, h, nseg)`: Bevels multiple edges with given width and height.
-   `bridge_quads(q1, q2, nseg, stickTypes, mask, bAlling)`: Creates a bridge between two quads.
-   `selectionToFace()`: Creates a new polygon face from the currently selected edge loop.
-   `polygon(n, iedges)`: Creates a polygon face from an array of edge indices.
-   `plateBetweenVertStrips(n, ivs1, ivs2, nsub)`: Creates a plate (subdivided surface) between two strips of vertices.
-   `plateBetweenEdges(nsub, r, bSort)`: Creates a plate between selected edges.
-   `conect_vertex(iv, stickType, n, iverts)`: Connects a vertex to a set of other vertices.
-   `conected_vertex(p, stickType, n, iverts)`: Connects a new vertex at `p` to a set of other vertices.
-   `make_anchor_point(p, stickType, Rcolapse, r, fw, l, i0, n)`: Creates an anchor point, potentially connecting to nearby vertices.
-   `make_anchor_points(nv, vs, ivrts, anchorType, Rcolapse, r, fw, l)`: Creates multiple anchor points.
-   `bondsBetweenVertRanges(v1s, v2s, Rmax, et)`: Creates bonds (edges) between vertices in two specified ranges.
-   `duplicateBlock(iblock)`: Duplicates an existing block of geometry (vertices, edges, triangles) and returns the index of the new block.
-   `snapBoxFace(p0, rot, La, Lb)`: Creates a simple, flat quad face, snapping to existing vertices.
-   `frustrumFace(p0, rot, La, Lb, h, Lbh, Lah)`: Creates a frustum-shaped face.
-   `snapFrustrumFace(p0, rot, La, Lb, h, Lbh, Lah, bFace)`: Creates a frustum-shaped face, snapping to existing vertices.
-   `prismFace(p0, rot, La, Lb, h, Lbh)`: Creates a prism-shaped face.
-   `snapPrismFace(p0, rot, La, Lb, h, Lbh, bFace)`: Creates a prism-shaped face, snapping to existing vertices.
-   `bridgeFacingPolygons(p0, p1, ch1, ch2, nseg, stickTypes, maks)`: Bridges two facing polygons between two points.
-   `bridgeFacingPolygons(nrod, edges, points, nseg, chs, stickTypes, maks)`: Bridges multiple facing polygons using edge and point data.

### 4.7 Debugging and I/O
-   `printSizes()`: Prints the number of elements in each data vector.
-   `printVerts()`, `printEdges()`, `printSelection()`: Print detailed information for debugging.
-   `printSelectedVerts()`: Prints detailed information about selected vertices.
-   `printChunkRange(ich, ich2)`: Prints information about chunks within a specified range.
-   `checkAllPointsConnected(bExit, bPrint)`: Checks if all points in the mesh are connected, with options to exit on failure or print details.
-   `export_pos(ps, i0, i1)`: Exports vertex positions to a Vec3d  array.
-   `export_pos(ps, i0, i1)`: Exports vertex positions to a float4 array.
-   `export_edges(eds, i0, i1)`: Exports edge data to an array.
-   `export_tris(tri, i0, i1)`: Exports triangle data to an array.
-   `write_obj(fname, mask)`: Exports mesh data to a Wavefront OBJ file with specified components.
-   `read_obj(fname, mask)`: Imports mesh data from a Wavefront OBJ file with specified components.

## 6. ToDo and Future Work

This section outlines potential future enhancements and additions to the `MeshBuilder2` class, particularly focusing on comprehensive mesh editing tools similar to those found in advanced 3D modeling software like Blender's BMesh, while keeping in mind the dual purpose of the engine for both volumetric physics simulations and surface-based rendering (radiosity/light-scattering).

### General Mesh Editing Operations

-   **Bevel (Vertices, Edges, Faces)**: Implement functions to bevel selected vertices, edges, or faces, creating new geometry with rounded or chamfered corners. This is crucial for creating more realistic and complex shapes for rendering.
-   **Extrude Along Normal/Path**: Enhance the existing extrude functionality to allow extrusion along a specified normal or a custom path, providing more control over the shape of the extruded geometry.
-   **Dissolve (Vertices, Edges, Faces)**: Implement functions to dissolve selected vertices, edges, or faces, removing them while attempting to preserve the surrounding geometry. This can simplify meshes and remove unnecessary detail.
-   **Subdivide (Edges, Faces)**: Add functionality to subdivide edges and faces, increasing mesh density and allowing for finer detail. This is particularly useful for smooth shading and high-resolution rendering.
-   **Merge/Collapse (Vertices, Edges)**: Implement tools to merge or collapse selected vertices or edges, reducing polygon count and simplifying mesh topology.
-   **Bridge Edge Loops**: Extend `bridge_quads` to handle more general edge loops, allowing for seamless connections between arbitrary open boundaries.
-   **Inset Faces**: Create new faces inset from existing ones, useful for creating borders or paneling effects.
-   **Loop Cut and Slide**: Implement a tool to insert new edge loops across a mesh, with the ability to slide them along the surface.
-   **Knife Tool**: A freehand cutting tool to divide faces and edges.

### SpaceCraft Building

#### Build Plates (shields, radiators)

- we want to create trinagular panel (plate, sail) in between two line segments (aka girders) starting at common point (corner)
- we select points of the mesh which are closer than `r` from each line segment using `select_in_cylinder` 
- we sort these vertexes according to the distance from conter, and take the minimum of number of vertexes in each selection (we need same number of points along both girders)
- we create trinagular grid defined by the to vertex-strips, i.e. connecting the vertexes in the same order as they diverge from corner using rope, the number of rope segments should be equal to index of the endpoint in the sorted selection
- we create edges connecting vertexes of the ropes to create triangular grid
- we also fill the trinagles between these edges by tris  

#### Build Magnetic Nozzle

- Magnetic nozzle is basically fan of triangular grids deformed by UV-function to parabolic shape (see `DrawUV.h`)
- we first typically create fan of girders (typically 6-8)



### Physics Simulation Specific Tools

-   **Volumetric Meshing Tools**: Develop more advanced tools for generating and manipulating volumetric meshes (e.g., tetrahedral meshes) for finite element analysis, beyond the current truss-based approach.
-   **Constraint-Based Editing**: Integrate mesh editing with physical constraints, allowing users to intuitively modify geometry while respecting physical properties or simulation requirements.
-   **Topology Optimization**: Tools to automatically optimize mesh topology for better simulation performance or accuracy.

### Radiosity and Light-Scattering Specific Tools

-   **UV Unwrapping and Editing**: Essential tools for preparing meshes for texture mapping and lightmap generation.
-   **Normal Editing**: Functions to manipulate vertex normals for improved shading and lighting calculations.
-   **Decimation/Simplification**: Algorithms to reduce polygon count while preserving visual fidelity, important for real-time rendering and performance.
-   **Mesh Repair Tools**: Functions to automatically detect and fix common mesh issues like non-manifold geometry, flipped normals, and holes, which are critical for accurate radiosity calculations.