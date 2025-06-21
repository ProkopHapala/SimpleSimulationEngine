# `Mesh::Builder2` - The Unified Mesh Toolkit

## 1. Overview

The `Mesh::Builder2` class is a powerful and versatile toolkit for the procedural generation, manipulation, and editing of 3D geometric data within the SimpleSimulationEngine. It is designed to be a central hub for creating everything from simple primitive shapes to complex, high-resolution truss structures.

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

*   **`chunks`** (`std::vector<Quat4i>`): This is the metadata or descriptor vector. Each `Quat4i` entry in this vector defines a single "chunk" and describes a specific segment within the `strips` vector.

The `Quat4i` for a chunk is interpreted as `{ i0, i_end, n, type }`:
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

### 3.1. Snapping and Finding Vertices

Many functions, especially the `snap...Face` family, rely on finding existing vertices rather than creating new ones.

-   **`findVert(pos, R_snapVert)`**: This method searches for a vertex within a given radius (`R_snapVert`) of a specified position. This is crucial for "welding" new geometry to an existing mesh.
-   **`findOrAddEdges(verts)`**: This method checks if an edge between two vertices already exists (using `findEdgeByVerts`). If not, it creates a new one. This prevents duplicate edges.

### 3.2. Face Connector Generation (`snap...Face`)

These functions are the bridge between the high-level `ConstructionBlock` and the low-level mesh. They are responsible for creating the physical connector geometry on the face of a block.

-   `snapBoxFace(p0, rot, La, Lb)`: Creates a simple, flat quad face. It uses `findVert` to locate the four corner vertices of the face on the existing block's mesh and then uses `findOrAddEdges` and `polygon` to create the face chunk.
-   `snapPrismFace(p0, rot, ...)`**: Creates a prism-shaped connector that provides two new faces for connections. It finds two vertices on the block and creates four new vertices to form the prism. It then generates two quad face chunks.
-   `snapFrustrumFace(p0, rot, ...)`**: Creates a more complex frustum-shaped connector that provides five new faces. It finds two vertices on the block and creates six new vertices to form the frustum, generating five quad face chunks.

The `bFace` parameter in these functions controls whether the connecting quad polygons are actually generated and stored as chunks. This is essential for the `ConstructionBlock` workflow, which needs these chunks to connect girders.

### 3.3. Interactive Editing

The builder provides a complete suite of tools for interactive mesh editing, driven by user input.

1.  **Set Mode**: The user selects a `SelectionMode` (`vert`, `edge`, or `face`) using the `I`, `O`, `P` keys in the demo app.
2.  **Pick**: The user clicks, and `pickSelect()` is called. This uses ray-tracing (`pickVertex`, `pickEdge`, `pickTriangle`) to find the element under the cursor. The picked element's index is stored.
3.  **Select**: The picked index is added to the `selection` vector and `selset` set. The user can select multiple elements if `bAdditiveSelect` is true. A selection box can also be used with `selectRect()`.
4.  **Operate**: The user presses a key to perform an action on the selection.
    -   **`F` key -> `selectionToFace()`**: If a loop of edges is selected, this function sorts them and creates a new polygon face.
    -   **`E` key -> `extrudeFace(ipick, ...)`**: If a face is selected, this function extrudes it. It gets the face's vertices from its chunk, creates new vertices offset along the normal, creates a new face, and calls `bridge_quads` to connect the old and new faces, forming the sides of the extrusion.

## 4. Method Reference (Grouped by Functionality)

### Primitive Creation
-   `vert(pos, ...)`: Adds a vertex.
-   `edge(v1, v2, ...)`: Adds an edge between two vertices.
-   `tri(v1, v2, v3, ...)`: Adds a triangle.
-   `chunk(...)`: Adds a metadata chunk pointing to a strip of primitives.
-   `quad(...)`: A helper to create a quad from four vertices (as two triangles).

### Complex Shape Generation
-   `box(p, ls, rot)`: Creates a box with 8 vertices and 12 edges.
-   `rope(ip0, ip1, ...)`: Creates a sequence of connected edges between two vertices.
-   `ring(...)`: Creates a circular loop of vertices and edges.
-   `tube(...)`: Creates a tube by connecting a series of rings.
-   `plate(...)`: Creates a subdivided flat surface between four corner points.
-   `girder1(...)`: Creates a complex, lightweight girder structure.
-   `wheel(...)`: Creates a wheel-shaped truss structure.
-   `panel(...)`: Creates a subdivided panel with truss-like internal structure.

### Mesh Editing and Topology
-   `extrudeFace(ich, L, ...)`: Extrudes a face chunk outwards by distance `L`.
-   `extrudeVertLoop(...)`: A lower-level function to extrude a loop of vertices.
-   `bridge_quads(q1, q2, ...)`: Creates a segmented bridge/tunnel between two quads.
-   `move_verts(...)`: Translates a specified set of vertices by a given shift vector.
-   `scale_verts(...)`: Scales a specified set of vertices relative to a pivot point `p` by scale factors `s`.
-   `rotate_verts(...)`: Rotates a specified set of vertices around a pivot point `p` by a given rotation matrix `rot`.
-   `duplicateBlock(...)`: Duplicates an existing block of geometry (vertices, edges, triangles) and returns the index of the new block.
-   `selectionToFace()`: Creates a new face from the currently selected edge loop.
-   `polygon(n, iedges)`: Creates a polygon face from an array of edge indices.

### Selection & Picking
-   `pickVertex(ray0, hRay, R)`: Finds the closest vertex to a ray.
-   `pickEdge(ro, rh, Rmax)`: Finds the closest edge to a ray.
-   `pickTriangle(ro, rh, ...)`: Finds the triangle intersected by a ray.
-   `pickSelect(ro, rh, Rmax)`: Main picking function that dispatches to the above based on `selection_mode`.
-   `selectRect(p0, p1, rot)`: Selects all elements within a screen-space rectangle.
-   `clearSelection()`: Clears the current selection.

### Topological Utilities
-   `findVert(p0, Rmax, ...)`: Finds the closest vertex to a point within a given radius.
-   `findEdgeByVerts(verts)`: Finds the index of an edge connecting two vertices. Can use a fast map (`_map`) or slow brute-force (`_brute`).
-   `sortEdgeLoop(n, iedges, ...)`: Sorts an array of edge indices into a contiguous loop.
-   `buildVerts2Edge()`: Populates the `vert2edge` map for fast lookups.

### Debugging and Export
-   `printSizes()`: Prints the number of elements in each data vector.
-   `printVerts()`, `printEdges()`, `printSelection()`: Print detailed information for debugging.
-   `export_pos()`, `export_edges()`: Export geometry data to external arrays.
-   `addPointCross(const Vec3d& p, double d)`: Adds a 3D cross shape at a given position `p` with size `d` for debugging visualization.
-   `addArrow(const Vec3d& p1, const Vec3d& p2, double d)`: Adds an arrow from `p1` to `p2` with a head size `d` for debugging visualization.

## 5. Code Quality & Refactoring Suggestions

The `MeshBuilder2` class is highly functional but could benefit from some refactoring to improve clarity and reduce code duplication.

### Suggestion 1: Refactor `snap...Face` Methods

The methods `snapPrismFace` and `snapFrustrumFace` in `MeshBuilder2.cpp` share a significant amount of logic for finding base vertices and creating new vertices and edges.

**Problem:** Code duplication makes maintenance harder and increases the chance of bugs.

**Proposed Solution:** Create a private helper function that takes the common parameters (base position, rotation, vertices to find, new vertex offsets) and handles the vertex creation and basic edge connections. The public `snap...Face` methods would then call this helper and add only the logic specific to their shape (i.e., creating the final polygon chunks).

```cpp
// Potential helper function signature in MeshBuilder2.h (private)
void snapConnectorBase(const Vec3d& p0, const Mat3d& rot, double R_snap,
                       const std::vector<Vec3d>& find_offsets, std::vector<int>& found_ivs,
                       const std::vector<Vec3d>& new_v_ps,   std::vector<int>& new_ivs);
```

### Suggestion 2: Improve Naming Clarity in `ConstructionBlock.h`

The function `replaceId` in `ConstructionBlock.h` has a name that does not clearly communicate its purpose. It actually *replaces* a logical edge ID with a geometric chunk ID and *returns* the old ID.

**Problem:** A confusing function name makes the code harder to understand, especially in the critical `replace_chunk` workflow.

**Proposed Solution:** Rename the function to something more descriptive, such as `assignChunkToSlot` or `mapEdgeToChunk`.

```diff
--- a/home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/ConstructionBlock.h
+++ b/home/prokophapala/git/SimpleSimulationEngine/cpp/common/geometry/ConstructionBlock.h
@@ -253,9 +253,11 @@
         return true;
     }
 
-    int replaceId( int id, Vec2i where ){
+    // Replaces the logical edge ID in a slot with a new ID (typically a geometric chunk ID)
+    // and returns the old logical ID that was stored there.
+    int assignChunkToSlot( int new_id, Vec2i where ){
         int oid = faces[where.i].ids[where.j];
-        if( oid<0 ) faces[where.i].nid++;
+        faces[where.i].ids[where.j] = new_id;
         return oid;
     }
 
```
*Note: The logic in the original `replaceId` seems to have a side effect of incrementing `nid` if the slot was empty, which might be a bug. The suggested change clarifies the primary intent.*