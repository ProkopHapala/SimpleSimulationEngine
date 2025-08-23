# Low-Resolution to High-Resolution Truss Generation

## 1. Overview

The primary objective of this system is to convert a simple, abstract graph representation of a spacecraft into a detailed, high-resolution, and physically plausible truss mesh. This allows designers to define complex structures using a simple set of nodes and edges, while the system automatically handles the geometric complexity of the joints and connections.

The process involves two main stages:
1.  **Low-Resolution Graph Definition**: A simple graph is created where vertices represent nodes (junctions) and edges represent the connections (girders or ropes).
2.  **High-Resolution Mesh Generation**: This graph is then processed to generate a detailed mesh using the `ConstructionBlock` and `MeshBuilder2` toolkit.

## 2. Data Representation (Low-Resolution Graph)

The low-resolution graph will be stored using a `Mesh::Builder2` instance, supplemented by an auxiliary data vector.

-   **`Mesh::Builder2 low_res_mesh`**:
    -   **`low_res_mesh.verts`**: Each vertex in this vector represents a node in the spacecraft's skeleton. Its `pos` attribute stores the 3D position of the node.
    -   **`low_res_mesh.edges`**: Each edge represents a connection between two nodes. The `Quat4i` structure `{v1, v2, type, type2}` is used to store the connection. We will utilize the `type` (`w` component) to distinguish between different kinds of connections:
        -   `type = 1`: A structural girder.
        -   `type = 2`: A flexible rope or cable.

-   **`std::vector<double> node_sizes`**:
    -   This vector runs in parallel to `low_res_mesh.verts`. `node_sizes[i]` stores the characteristic size (e.g., radius or half-extent) of the `ConstructionBlock` that will be generated at the location of `low_res_mesh.verts[i]`. This is necessary because `Mesh::Builder2` does not have a native way to store per-vertex size attributes.

This combined representation is lightweight and sufficient to describe the entire skeleton of the spacecraft.

## 3. The Conversion Process

The conversion is handled by a function, let's call it `convertLowToHigh`, which takes the low-resolution graph and outputs a new, high-resolution `Mesh::Builder2` instance. The process is as follows:

**Input:**
-   `const Mesh::Builder2& low_res_mesh`
-   `const std::vector<double>& node_sizes`

**Output:**
-   `Mesh::Builder2& high_res_mesh`

### Step-by-Step Algorithm:

1.  **Initialize `BlockBuilder`**:
    -   Create an empty `BlockBuilder skelet`. This will serve as the high-level blueprint for the structure.

2.  **Create Blocks (Nodes)**:
    -   Iterate through the vertices of the `low_res_mesh` from `i = 0` to `n-1`.
    -   For each vertex `i`, call `skelet.addBlock(low_res_mesh.verts[i].pos, node_sizes[i])`.
    -   This creates a `ConstructionBlock` for each node. The index of the block in `skelet.blocks` will directly correspond to the vertex index `i` from the low-resolution mesh, which is a crucial link.

3.  **Create Girder Connections**:
    -   Iterate through the edges in `low_res_mesh.edges`.
    -   For each edge `e` connecting vertices `e.x` and `e.y`:
        -   Check if the edge `type` (`e.w`) signifies a **girder**.
        -   If it is a girder, call `skelet.connectBlocks(e.x, e.y)`. This function automatically finds the appropriate faces on the two `ConstructionBlock`s and registers the logical connection.
    -   *Note: Ropes are intentionally skipped in this step, as they are not part of the `BlockBuilder`'s rigid connection system.*

4.  **Generate High-Resolution Blocks and Girders**:
    -   Create an instance of `Mesh::ConstructionBlockToMeshBuilder cbm`.
    -   Assign the output mesh: `cbm.mesh = &high_res_mesh`.
    -   Call `cbm.drawBlockBuilder(skelet)`. This is the core generation step. It will:
        -   Create the geometry for each `ConstructionBlock` (including any special face connectors).
        -   Create the detailed truss geometry for each girder by calling `mesh.bridge_quads` between the appropriate block faces.

5.  **Add Ropes to High-Resolution Mesh**:
    -   After the rigid structure is built, iterate through `low_res_mesh.edges` a second time.
    -   For each edge `e` where the `type` signifies a **rope**:
        -   Retrieve the positions of the two connected nodes: `p0 = low_res_mesh.verts[e.x].pos` and `p1 = low_res_mesh.verts[e.y].pos`.
        -   Call `high_res_mesh.rope(p0, p1, n_segments)`, where `n_segments` is a desired level of subdivision for the rope. This adds the rope geometry directly to the final mesh, connecting the centers of the generated blocks.

## 4. Conclusion

This two-pass approach cleanly separates the generation of the rigid frame (blocks and girders) from the addition of secondary elements (ropes). It effectively leverages the `BlockBuilder` system for its intended purpose—complex rigid connections—while handling simpler elements directly within `MeshBuilder2`. The result is a robust and extensible pipeline for creating detailed spacecraft meshes from a simple abstract definition.

---