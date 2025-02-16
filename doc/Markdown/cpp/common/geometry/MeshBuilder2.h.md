# MeshBuilder2.h

This file defines the `Mesh::Builder2` class, a versatile tool for constructing and manipulating mesh data. It consolidates functionalities previously found in separate classes like `MeshBuilder`, `Truss`, `OMesh`, and `GLMeshBuilder`, providing a unified interface for creating various geometric structures, including meshes and truss systems. The class supports different primitive types, vertex attributes (position, normal, UV coordinates), and selection mechanisms for editing and manipulation.

## Includes

- `<vector>`: Provides dynamic array capabilities.
- `<unordered_map>`: Provides hash table implementation for efficient key-value storage.
- `<unordered_set>`: Provides hash set implementation for storing unique elements.
- `fastmath.h`: Provides optimized mathematical functions.
- `Vec2.h`: Defines the `Vec2` class for 2D vector operations.
- `Vec3.h`: Defines the `Vec3` class for 3D vector operations.
- `Mat3.h`: Defines the `Mat3` class for 3x3 matrix operations.
- `quaternion.h`: Defines the `Quat4d` and `Quat4f` classes for quaternion operations, used for representing rotations.
- `datatypes.h`: Defines fundamental data types used throughout the project.
- `Slots.h`: Defines the `Slots` class for managing a fixed number of elements, used for efficient storage and retrieval.
- `MeshBuilder2.h`: (Self-inclusion) Declares the `Mesh::Builder2` class.
- `<algorithm>`: Provides a collection of functions implementing algorithms.
- `arrayAlgs.h`: Defines various array-based algorithms.
- `raytrace.h`: Defines functions for ray tracing, used for picking and collision detection.
- `testUtils.h`: Includes utilities for testing and debugging the code.

---
## Types (classes and structs)
---

### class `Mesh::Builder2`

The `Mesh::Builder2` class is the core class for constructing and manipulating mesh data. It provides a unified interface for creating various geometric structures, including meshes and truss systems. It supports different primitive types, vertex attributes (position, normal, UV coordinates), and selection mechanisms for editing and manipulation.

#### properties

##### Mesh Data Structures
- `blocks`: `std::vector<Quat4i>` - Stores blocks of data, where each block represents a complex object and contains the starting indices for vertices, edges, triangles, and chunks.
- `verts`: `std::vector<Vert>` - Stores the vertices of the mesh, including position, normal, and UV coordinates.
- `edges`: `std::vector<Quat4i>` - Stores the edges of the mesh, represented by the indices of the two connected vertices and optional type information or adjacent face indices.
- `tris`: `std::vector<Quat4i>` - Stores the triangles of the mesh, represented by the indices of the three vertices and optional type information.
- `chunks`: `std::vector<Quat4i>` - Stores chunks of primitives, which can represent polygons, line strips, or triangle strips, and contains the starting index, number of primitives, and type information.
- `strips`: `std::vector<int>` - Stores the indices of primitives within the chunks, providing the connectivity information for the mesh.

##### Vertex Attributes
- `Vert`: `using` - Alias for `VertT<double>`, representing a vertex with double-precision position, normal, and UV coordinates.

##### Configuration Flags
- `use_vert2edge`: `bool` - Flag indicating whether to use a hash map for efficient edge lookup by vertex pairs.
- `bPolygonToTris`: `bool` - Flag indicating whether to automatically triangulate polygons upon creation.
- `bExitError`: `bool` - Flag indicating whether to exit the program upon encountering an error.
- `bAdditiveSelect`: `bool` - Flag indicating whether selection operations should add to or replace the current selection.

##### Selection
- `selection_mode`: `int` - Integer representing the selection mode (0-none, 1-vert, 2-edge, 3-face).
- `selection`: `std::vector<int>` - Stores the indices of selected vertices, edges, or faces, depending on the selection mode.
- `selset`: `std::unordered_set<int>` - Stores the indices of selected elements in a hash set for efficient membership testing.

##### Geometric Parameters
- `max_size`: `double` - Maximum size used for adaptive tessellation or subdivision.
- `R_snapVert`: `double` - Tolerance for snapping vertices together during operations like finding or merging vertices.

##### Primitive Types
- `face_type`: `int` - Default type assigned to newly created faces.
- `edge_type`: `Vec2i` - Default type assigned to newly created edges, potentially representing UV edge types.

##### State Variables
- `ov`: `int` - Index of the previous vertex, used for creating connected primitives like strips.
- `oov`: `int` - Index of the vertex before the previous vertex, used for creating connected primitives like strips.

##### Limits
- `ngon_max`: `int` - Maximum number of vertices allowed in a polygon before a warning is issued.

##### Rendering
- `penColor`: `Vec3f` - Color used for drawing the mesh.

#### methods

##### Mesh Element Creation
- `block`: Creates a new block and returns its index.
- `vert`: Creates a new vertex with the given position, normal, and UV coordinates and returns its index.
- `edge`: Creates a new edge connecting two vertices with optional type information and returns its index.
- `tri`: Creates a new triangle from three vertices with optional type information and returns its index.
- `stick`: Creates a new edge (stick) between two specified 3D points, adding the points as new vertices.
- `stickTo`: Creates a new edge (stick) from an existing vertex to a new vertex at a specified 3D point.
- `edgst`: Creates an edge from the previous vertex (`ov`) to a new vertex, updating `ov` and returning the edge index.
- `trist`: Creates a triangle from the previous vertex (`ov`) and a new vertex, updating `oov` and `ov` and returning the triangle index.
- `chunk`: Creates a new chunk with the given parameters and returns its index.
- `quad`: Creates a quad (two triangles) from four vertices, assigning face and edge types.

##### Mesh Element Access
- `getOtherEdgeVert`: Returns the index of the vertex at the other end of an edge, given the edge index and one vertex index.
- `getChunkStrip`: Returns a pointer to the array of primitive indices within a specified chunk.
- `getChunkNormal`: Calculates and returns the normal vector of a face chunk.

##### Mesh Element Export
- `export_pos`: Exports vertex positions to an array of `Vec3d` objects.
- `export_pos`: Exports vertex positions to an array of `float4` objects.
- `export_edges`: Exports edge vertex indices to an array of `Vec2i` objects.
- `export_tris`: Exports triangle vertex indices to an array of `Quat4i` objects.

##### Mesh Building Blocks
- `rope`: Creates a series of connected edges (a rope) between two existing vertices.
- `ring`: Creates a ring of vertices around a center point, connecting them with edges.
- `tube`: Creates a tube-like structure by creating rings along a line segment.
- `plate`: Creates a flat surface (plate) between four corner points, subdividing it into smaller elements.
- `ngon`: Creates an n-sided polygon (ngon) centered between two points.
- `panel`: Creates a panel of truss elements between four corner points.

##### Mesh Modification
- `extrudeVertLoop`: Extrudes a loop of vertices along a direction vector, creating new vertices, edges, and faces.
- `extrudeFace`: Extrudes a face (chunk) along its normal vector, creating new vertices, edges, and faces.
- `plateBetweenVertStrips`: Creates a plate (surface) between two vertex strips.
- `plateBetweenEdges`: Creates a plate (surface) between two selected edges.
- `alling_polygons`: Aligns two polygons by rotating and translating one to match the other.
- `bridge_quads`: Creates a bridge (tunnel) between two quads, connecting corresponding vertices with edges.

##### Selection and Picking
- `selectVertsAlongLine`: Selects vertices that lie along a line segment within a specified radius.
- `selectVertsAlongPolyline`: Selects vertices that lie along a polyline within a specified radius.
- `select_in_cylinder`: Selects vertices within a cylinder defined by a center point, direction, radius, and length.
- `select_in_box`: Selects vertices within a box defined by a center point, forward vector, up vector, minimum lengths, and maximum lengths.
- `make_anchor_point`: Creates an anchor point by connecting a new vertex to all selected vertices within a cylinder.
- `pickVertex`: Picks the vertex closest to a ray.
- `pickEdge`: Picks the edge closest to a ray.
- `pickTriangle`: Picks the triangle intersected by a ray.
- `pickSelect`: Picks and selects an object (vertex, edge, or face) based on the current selection mode.
- `selectRectEdge`: Selects edges within a rectangular region.
- `selectRectVert`: Selects vertices within a rectangular region.
- `selectRect`: Selects objects (vertices or edges) within a rectangular region based on the current selection mode.
- `findClosestVert`: Finds the vertex closest to a given point.
- `findVert`: Finds a vertex near a given point within a specified radius.
- `toggleSelSet`: Toggles the selection state of a vertex, adding it to or removing it from the selection set.
- `selectionToFace`: Converts the current vertex selection to a face.
- `clearSelection`: Clears the current selection.
- `makeSelectrionUnique`: Removes duplicate indices from the selection list.

##### Edge Management
- `sortPotentialEdgeLoop`: Sorts a set of edges to form an edge loop, based on shared vertices.
- `sortEdgeLoop`: Sorts a set of edge indices to form an edge loop, based on shared vertices.
- `findEdgeByVerts_brute`: Finds an edge by its two vertices using a brute-force search.
- `findEdgeByVerts_map`: Finds an edge by its two vertices using a hash map.
- `findEdgeByVerts`: Finds an edge by its two vertices, using either a brute-force search or a hash map.
- `findOrAddEdges`: Finds an edge by its two vertices, and creates it if it doesn't exist.
- `buildVerts2Edge`: Builds a hash map for efficient edge lookup by vertex pairs.

##### Geometric Primitives
- `box`: Creates a box (cuboid) with specified dimensions and orientation.
- `snapBoxFace`: Snaps a box face to existing vertices within a tolerance.
- `frustrumFace`: Creates a frustum face with specified dimensions and orientation.
- `snapFrustrumFace`: Snaps a frustum face to existing vertices within a tolerance.
- `prismFace`: Creates a prism face with specified dimensions and orientation.
- `snapPrismFace`: Snaps a prism face to existing vertices within a tolerance.

##### Truss Structures
- `girder1`: Creates a girder (beam) between two points with specified width and subdivisions.
- `triangle_strip`: Creates a flat zig-zag triangle strip made of equal sized triangles.
- `plateOnGriders`: Creates a plate on girders.
- `girder1_caps`: Creates caps at the ends of a girder.
- `wheel`: Creates a wheel-shaped truss structure.
- `rope`: Creates a rope-like structure between two points.
- `panel`: Creates a panel of truss elements between four corner points.

##### Debugging and Information
- `printSelection`: Prints the current selection of vertices, edges, or faces.
- `printSelectedVerts`: Prints the positions of the selected vertices.
- `printSizes`: Prints the sizes of the main data structures (vertices, edges, triangles, chunks, strips).
- `printVerts`: Prints the positions of all vertices.
- `printEdges`: Prints the connectivity information for all edges.

##### Low Level
- `latsBlock`: Returns a `Quat4i` containing the current sizes of the vertex, edge, triangle, and chunk lists.