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