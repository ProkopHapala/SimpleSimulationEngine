 # Spacecraft Construction with Blocks: Enhancing Structural Rigidity

 ## 1. Introduction

 This document outlines a new approach to spacecraft construction within the SimpleSimulationEngine. The primary goal is to improve the structural rigidity of spacecraft designs, addressing limitations observed with the previous single-node connection method.

 ## 2. Motivation: Limitations of Single-Node Connections

 The original spacecraft construction system relied on connecting structural components (girders, ropes, etc.) at single points (Nodes). While this approach offered simplicity in design, it resulted in structures that were prone to wobbling, swaying, and oscillations during simulation. This lack of rigidity made it difficult to stabilize spacecraft and achieve predictable behavior.

 ## 3. New Approach: Block-Based Construction

 To address these limitations, we introduce a new construction paradigm that utilizes *blocks* as connection hubs. These blocks, typically Platonic solids (e.g., cubes, octahedrons), are represented by small meshes and provide multiple connection points on their faces. This allows for more robust and distributed connections between structural components, significantly increasing the overall rigidity of the spacecraft.

 ## 4. Key Components and Functions

 *   **`ConstructionBlock`** ( `ConstructionBlock.h` ): Represents a single block in the structure. It defines the block's geometry and manages the connection points on its faces.
 *   **`BlockBuilder`** ( `BlockBuilder.h` ): Manages the overall structure of interconnected blocks.
 *   **`Mesh::Builder2`** ( `MeshBuilder2.h` ): The core class for building the final mesh.
 *   **`Mesh::ConstructionBlockToMeshBuilder`** ( `ConstructionBlock.h` ): Converts the `BlockBuilder` structure into a detailed `Mesh::Builder2` mesh.

 ## 5. Crucial Functions for Block-Based Construction

 *   **`Mesh::Builder2::extrudeFace`** ( `MeshBuilder2.h` ): Extrudes a face of a mesh, creating a new connected volume. This is used to create connector geometries on the block faces.
 *   **`Mesh::Builder2::bridge_quads`** ( `MeshBuilder2.h` ): Creates a bridge (a connecting tunnel) between two quad faces. This is used to connect girders between blocks.
 *   **`Mesh::Builder2::make_anchor_point`** ( `MeshBuilder2.h` ): Creates a new vertex and connects it to nearby vertices, providing a robust anchor point for other components.
 *   **`Mesh::Builder2::make_anchor_points`** ( `MeshBuilder2.h` ): Creates multiple anchor points.

 ## 6. Anchoring Mechanism: `make_anchor_point` and `make_anchor_points`

 These functions are essential for attaching ropes and potentially other structural elements to the block-based structure.

 *   **Snapping (Collapsing):** If a desired anchor point is located within a certain tolerance (`Rcolapse`) of an existing vertex, the new anchor point will *snap* to the existing vertex, avoiding the creation of duplicate vertices.
 *   **New Vertex Creation:** If no suitable vertex is found nearby, a new vertex is created at the desired anchor point. This new vertex is then connected to nearby vertices within a specified region (spherical or cylindrical) using new edges. This distributes the connection forces and provides a more stable anchor.

 ## 7. Mapping Vertex and Edge Ranges: Maintaining Component Relationships

 A critical aspect of the system is maintaining a reliable mapping between high-level spacecraft components (e.g., girders, ropes, modules) and the low-level mesh data (vertices and edges). This mapping is essential for:

 *   **Geometry Updates:** When the shape or configuration of a component changes, the corresponding mesh data needs to be updated efficiently.
 *   **Collision Detection:** The physics engine needs to quickly identify the mesh elements associated with a given component to accurately simulate collisions.
 *   **Interactive Editing:** The editor needs to highlight or manipulate the mesh elements corresponding to a selected component.

 In the previous system, this mapping was achieved by storing the starting and ending indices of the vertices and edges belonging to each component (`pointRange` and `stickRange` properties of `ShipComponent` class in `SpaceCraftComponents.h`). With the new block-based construction method, we need to ensure that these ranges are still correctly assigned, especially when using functions like `extrudeFace` and `bridge_quads`, which generate geometry in a different way than the old component-specific functions.

 ## 8. Problem Analysis: Challenges and Considerations

 *   **Finding the Nearest Side:** When attaching a component to a block, we need to determine the "nearest side" of the block to ensure a consistent and predictable attachment. The existing `nearSide` methods in `StructuralComponent` need to be adapted or extended to work with the new block-based geometry.
 *   **Optimal Slider Anchor Position:** For components like sliders (e.g., wheels attached to girders), finding the optimal anchor position is crucial for proper functionality. The existing analytical calculations for finding the intersection point between a circle (wheel) and a line segment (girder) may need to be adapted to work with the new block-based geometry.

 ## 9. Next Steps

 This document provides a starting point for understanding the new block-based construction system. The next steps involve:

 *   **Implementing the Mapping of Vertex/Edge Ranges:** Develop a robust and reliable method for assigning vertex and edge ranges to components created using `extrudeFace` and `bridge_quads`.
 *   **Adapting `nearSide` Methods:** Adapt or extend the existing `nearSide` methods to work with the new block-based geometry.
 *   **Developing Slider Anchor Calculations:** Develop analytical or numerical methods for finding the optimal anchor position for sliders attached to block-based structures.
 *   **Testing and Validation:** Thoroughly test the new system to ensure that it produces stable and predictable results.

---

## 10. Concrete Problem: Attaching a Wheel to a Block-Based Hull

This section provides a focused analysis of a specific, practical challenge: attaching a rotating wheel (a `Ring` component) to a hull constructed from `ConstructionBlock` nodes and bridged girders. This process will serve as a template for attaching other complex, movable components like sliders.

### 10.1. The Challenge

The old system attached wheels to simple, line-like girders using the `l_Ring2` Lua function. This function relied on two key operations:
1.  `circle_3point`: To define the wheel's geometry from three initial anchor points.
2.  `intersect_RingGirder`: To find the precise attachment points for the sliders by analytically calculating the intersection of the wheel's circle with the girder's line segment.

With the new block-based system, the "girders" are no longer simple lines but complex truss structures generated by `bridge_quads`. This presents two primary problems:
1.  **Lack of Component-to-Mesh Mapping:** The current block-building functions (`facingNodes`, `bridgeFacingPolygons`) do not create high-level `Girder` components or store their corresponding `pointRange` and `stickRange` in the `Mesh::Builder2` mesh. Without this mapping, we cannot programmatically identify which part of the mesh constitutes a specific girder.
2.  **Identifying Intersection Paths:** The `intersect_RingGirder` function is still highly relevant, as the new girders are fundamentally straight struts. However, a girder generated by `bridge_quads` is a 3D truss with four main longitudinal "corner" struts. We need a mechanism to identify which of these four struts should be treated as the "line segment" for the intersection calculation.

### 10.2. Step 1: Tracking Girder Geometry (`pointRange` and `stickRange`)

The first and most critical step is to ensure that when we build the hull, we also create corresponding high-level `Girder` components and track their geometry. The `bridgeFacingPolygons` function needs to be wrapped in a higher-level function that manages this mapping.

**Proposed Workflow:**

Instead of just calling `mesh.bridge_quads` directly, we need a new workflow that creates a `Girder` object for each bridged connection.



By implementing this, we re-establish the crucial link between the high-level `Girder` component and its low-level mesh representation.

### 10.4. Re-implementing `BuildCraft_truss` for Block-Based Construction

**Design Goal:** The primary goal is to re-implement the `BuildCraft_truss` function (located in `cpp/common/Orbital/SpaceCraft2Mesh2.h`) to support the new block-based construction paradigm. This re-implementation will draw inspiration from the mesh generation logic demonstrated in `cpp/apps/OrbitalWar/constructionBlockApp.cpp` within the `funcs["-oct_nodes"]` section.

**Key Principles:**
*   **Decoupling:** Maintain the clear separation between the `SpaceCraft` object (which defines the high-level structure of nodes and girders) and the `Mesh::Builder2` (which handles the low-level mesh generation).
*   **Iteration:** `BuildCraft_truss` will iterate through the `craft.nodes` (where `Node` objects now represent `ConstructionBlock`s or similar block-like entities) and `craft.girders` (representing connections between these blocks).
*   **Mesh Generation:** Within `BuildCraft_truss`, for each `Node` and `Girder` object, the corresponding mesh geometry will be generated using `Mesh::Builder2` functions (e.g., `truss.facingNodes`, `truss.bridgeFacingPolygons`).
*   **Property Assignment:** Crucially, the `pointRange` and `stickRange` properties of `Node` and `Girder` objects will be assigned during this mesh generation phase within `BuildCraft_truss`, reflecting the actual vertex and edge ranges in the generated mesh.

This approach ensures that the `SpaceCraft` object remains a high-level, abstract representation of the spacecraft structure, while `BuildCraft_truss` is responsible for translating that abstraction into a concrete mesh.

### 10.3. Step 2: Adapting Girder Functionality for Path Extraction

We do not need to redefine the Girder class. Instead, we need to provide a concrete implementation for its existing virtual methods, specifically `sideToPath`, that can work with the new mesh topology generated by `bridge_quads`.

A standard girder created by `bridge_quads` has four main longitudinal "corner" struts. These are the natural paths for sliders or for line-based intersection tests. The existing `sideToPath` concept from `StructuralComponent` is the perfect way to extract these.

**Proposed `Girder::sideToPath` Implementation:**

```cpp
// In the new Girder class, derived from StructuralComponent

virtual int sideToPath(int side, int*& inds, ...) const override {
    // 'side' would be an index from 0 to 3, for the four corner paths.
    // This method needs to traverse the girder's mesh (using its pointRange
    // and stickRange) to find the ordered sequence of vertex indices
    // that form the specified corner path.
    // This is a graph traversal problem on a subset of the mesh defined by the ranges.
    // The implementation will be complex but is essential.
    // It would start at one of the four vertices of the starting quad face
    // and walk along the longitudinal edges to the corresponding end face.
    // The traversal should be efficient, possibly using a breadth-first search (BFS) algorithm.
    // The BFS would start from the initial vertex and explore all neighboring vertices
    // until it reaches the end face, keeping track of the visited vertices to avoid loops.
    // ...
    // return number of vertices in the path
}
```

### 10.4. Slider Anchoring Logic

It is crucial to understand the fundamental behavior of sliders. A slider is typically attached at a fixed point of one structural component and slides along a path defined by the vertices of another structural component. This attachment mode is inherently asymmetric with respect to the two objects it connects.

In the specific case of wheel attachment, we aim to firmly attach 3-4 sliders to the girders, allowing them to slide along the wheel (i.e., the path is made of vertices of the wheel). This mechanism is consistent with how it functioned in the old ship-building system, emphasizing continuity and reuse of established principles.

### 10.5. Step 3: The New Wheel Attachment Workflow

With the above pieces in place, we can define the full workflow for attaching a wheel.

The attachment logic is asymmetric: the sliders are fixed to the girders and slide along the wheel. This is the same principle as in the old system.

1.  **Create the Ring Component:** A `Ring` component is created, and its mesh is generated using `mesh.wheel()`. Its `pointRange` and `stickRange` are stored. The `Ring`'s `sideToPath` method will provide the circular paths for the sliders.

2.  **Find Attachment Points on Girders:** For each of the 3-4 sliders that will connect the wheel to the hull:
    a.  Identify the target `Girder` on the hull.
    b.  Use the girder's new `sideToPath()` implementation to get one of its main longitudinal struts (e.g., side=0). This strut is a polyline (a sequence of vertices).
    c.  Approximate this strut as a single line segment from its start to end point.
    d.  Use the existing `intersect_RingGirder` function. It takes the `Ring` and the `Girder` (from which the line segment is derived) and calculates the precise intersection point in 3D space. This point lies on the girder.

3.  **Create Stable Anchors on Girders:**
    a.  At each of the precise intersection locations found above, call `mesh.make_anchor_point(pos, ...)`.
    b.  This function creates a robust physical attachment point on the girder's truss structure by adding a new vertex and connecting it to nearby existing vertices of the girder with new sticks. It returns the index of the new primary vertex (`ivert`). This is the fixed point for the slider.

4.  **Create and Configure Sliders:**
    a.  Create a `Slider` component for each anchor point.
    b.  The slider's position is now fixed in the mesh. Set its `ivert` to the vertex index returned by `make_anchor_point`.
    c.  The slider's `boundTo` property can be `nullptr`, as its position is absolute within the mesh, not relative to another component.
    d.  Initialize the slider's path of movement by calling `sideToPath` on the `Ring` component. This populates `slider.path.ps` with the vertex indices of the wheel's circumference.

This workflow correctly uses the `Ring` for the sliding path and the `Girder` for the fixed attachment points, reusing `intersect_RingGirder` and adapting the `Girder` class through its virtual methods.
5.  **Create Final Anchors and Sliders:**
    a.  At each of the 4 precise locations found above, call `mesh.make_anchor_point()` again. This creates the final, robust physical attachment points on the girder.
    b.  Create a `Slider` component for each anchor.
    c.  The `Slider`'s `ivert` is set to the vertex index returned by `make_anchor_point`.
    d.  The `Slider`'s `path` is initialized by calling `sideToPath` on the girder it's attached to. Its initial `path.cur` is set based on the intersection calculation.

### 10.6. Step 4: Rethinking Node and Slider Anchoring

A key challenge is tracking the geometry created by `mesh.make_anchor_point`. This function adds new vertices and edges to the mesh, and this geometry should be associated with a high-level `ShipComponent` to maintain a complete mapping.

The `Node` class, inheriting from `Object`, lacks `pointRange` and `stickRange`. While changing its inheritance is possible, a cleaner solution is to introduce a new component that explicitly represents the physical structure of an anchor.

**Proposed `Anchor` Component:**

We can create a new class `Anchor : public StructuralComponent`.
-   When `mesh.make_anchor_point()` is called, it not only creates the mesh geometry but also a corresponding `Anchor` object.
-   This `Anchor` object would store the `pointRange` and `stickRange` of the supporting sticks created by `make_anchor_point`.
-   It would also store the index of its primary vertex (`ivert`), which is the actual attachment point.
-   A `Slider`'s `ivert` would be set to this primary vertex. The `Anchor` component itself simply serves to own the associated mesh geometry.

This approach avoids modifying the `Node` class, keeps the concepts clean, and correctly maps all generated geometry to high-level components, fulfilling a core requirement of the new system.

This approach avoids modifying the `Node` class, keeps the concepts clean (a `Slider` binds to an `Anchor`, which is a physical part), and correctly maps all generated geometry to high-level components.