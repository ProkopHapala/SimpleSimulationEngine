 # Spacecraft Building System: Architectural Overview

 This document aims to provide a high-level understanding of the SimpleSimulationEngine's spacecraft building system. It outlines the core components, their relationships, and the workflow from high-level design to low-level physical simulation. Special attention is given to the anchoring mechanisms, particularly the "Slider" component, as it represents a complex interaction between the design and simulation layers.

 ## 1. Overall Architecture: The Big Picture

 The spacecraft system is structured in layers, moving from abstract design to concrete simulation data:

 *   **High-Level Design (Lua/`SpaceCraft.h`):** Defines the logical structure and components of the spacecraft.
 *   **Conversion Layer (`SpaceCraft2Mesh2.h`):** Translates the logical spacecraft into a geometric mesh suitable for simulation.
 *   **Low-Level Simulation (`Mesh::Builder2`, `TrussDynamics_d/f`):** Represents the spacecraft as a collection of interconnected points and sticks for physics calculations.

 ## 2. Core Components: Building Blocks of a Spacecraft

 The fundamental building blocks are defined in `SpaceCraftComponents.h` and managed by the central `SpaceCraft` class in `SpaceCraft.h`.

 ### 2.1. `SpaceCraft` Class (`SpaceCraft.h`)

 This is the central container for all spacecraft components. It holds `std::vector`s of pointers to various `ShipComponent` derived types.

 **Key Responsibilities:**
 *   Manages collections of structural (girders, ropes, rings) and functional (radiators, thrusters, tanks) components.
 *   Provides methods (`add_Node`, `add_Rope`, `add_Girder`, etc.) to programmatically add components.
 *   Includes utility functions for geometric queries (e.g., `rayPlate`, `rayLinkLine` for picking).
 *   Maintains a `SpaceCraftWorkshop` for material definitions.
 *   Crucially, it manages `Node` objects, which serve as connection points, and `Slider` objects, which represent movable connections.

 ### 2.2. Component Hierarchy (`SpaceCraftComponents.h`)

 The system uses a class hierarchy to categorize and manage different types of spacecraft parts.

 *   **`Object`**: Base class (provides `id`, `kind`).
 *   **`ShipComponent`**: Inherits from `Object`. Base for all spacecraft parts.
     *   **Properties**: `id`, `kind`, `shape`, `face_mat` (material ID), `mass`, `pointRange` (start/end vertex index in mesh), `stickRange` (start/end edge index in mesh).
 *   **`StructuralComponent`**: Inherits from `ShipComponent`. Base for components that form the physical structure and can have `Node`s attached.
     *   **Properties**: `nodes` (a `vec4<Node*>` for up to 4 connected nodes), `mvert` (verts per segment).
     *   **Methods**: `update_nodes`, `updateSlidersPaths`, `findNearestPoint`, `toBuckets` (for spatial partitioning).
 *   **`NodeLinker`**: Inherits from `StructuralComponent`. Base for linear structural elements.
     *   **Derived Classes**:
         *   **`Girder`**: Rigid structural beam.
         *   **`Rope`**: Flexible cable-like structure.
 *   **`Ring`**: Inherits from `StructuralComponent`. Circular structural element.
 *   **`Plate`**: Inherits from `ShipComponent`. Flat surface components.
      *   **Derived Classes**:
          *   **`Radiator`**: Heat dissipation panel.
          *   **`Shield`**: Protective panel.
          *   **`Collector`**: Resource collection panel.
 *   **`Modul`**: Inherits from `ShipComponent`. Modular components with defined volume/pose.
      *   **Derived Classes**:
        *   **`Tank`**: Storage for resources.
        *   **`Balloon`**: Inflatable structure.
        *   **`Rock`**: Asteroid/debris shielding.
        *   **`Thruster`**: Propulsion system.
 *   **`Weld`**: Inherits from `ShipComponent`. Represents connections between structural components.
 *   **`Node`**: Inherits from `Object`. A crucial connection point in the structure.
     *   **Properties**: `ivert` (vertex index in the mesh), `pos`, `boundTo` (pointer to `StructuralComponent` if bound), `calong` (position along `boundTo`), `along` (vertex index along `boundTo`).
     *   **Methods**: `updateBound` (updates position/vertex index based on binding).
 *   **`Slider`**: Inherits from `Node`. Represents a movable connection along a path. (Detailed below in Anchoring section).
 *   **`Accelerator`**: Inherits from `ShipComponent`. Base for projectile/particle propulsion.
     *   **Derived Class**:
         *   **`Gun`**: Weapon system.

 ### 2.3. `SpaceCraftWorkshop` (`SpaceCraftComponents.h`)

 Manages catalogs of materials, commodities, and component types (e.g., `Material`, `StickMaterial`, `PanelMaterial`, `ThrusterType`, `GunType`). This allows for consistent material properties across the spacecraft.

 ## 3. High-Level Construction: Lua Interface (`EditSpaceCraft.h`)

 The `EditSpaceCraft.h` file provides the C++ interface for Lua scripting, allowing spacecraft to be defined using Lua scripts (e.g., `data/ship_ICF_marksman_2.lua`). See also `spacecraft.def.lua` for details about Lua side of the interface listing all the available functions and their parameters.

 ### 3.1 Key Lua Functions (and their C++ counterparts)
 *   `Material`, `StickMaterial`: Define material properties.
 *   `Node(...)`: Creates a free `Node` at a given position.
 *   `BoundNode(...)`: Creates a `Node` bound to an existing `StructuralComponent` at a relative position (`calong`). `p0` is a hint for finding the nearest side.
 *   `Rope(...)`: Creates a `Rope` between two `Node`s.
 *   `Rope2(...)`: Creates a `Rope` between points on existing `Girder`s. This is more complex as it can create new `BoundNode`s.
 *   `Girder(...)`: Creates a `Girder` between two `Node`s.
 *   `Ring(...)`: Creates a free `Ring`.
 *   `Ring2(...)`: Creates a `Ring` anchored to existing `Girder`s. This is where `Slider`s are often created implicitly to manage the attachment points.
 *   `Slider(...)`: Creates a `Slider` bound to a `StructuralComponent`.
 *   `Weld(...)`: Creates a `Weld` between two `StructuralComponent`s.
 *   `Radiator`, `Shield`, `Tank`, `Thruster`, `Gun`, `Rock`, `Balloon`: Create respective components.

 ## 4. Conversion to Simulation Mesh (`SpaceCraft2Mesh2.h`)

 This is the crucial step where the high-level `SpaceCraft` object is translated into a low-level geometric mesh suitable for physics simulation. The `Mesh::Builder2` class (`MeshBuilder2.md`) is the target data structure.

 ### 4.1. `BuildCraft_truss` Function

 This function takes a `Mesh::Builder2` object, a `SpaceCraft` object, and optionally a `max_size` for mesh elements.

 **Workflow:**
 1.  **Nodes:** Iterates through `craft.nodes`. If a `Node` is not `boundTo` another component, a new vertex is created in `mesh.verts` and its `ivert` (vertex index) is stored in the `Node` object. If it is bound, its `ivert` will be set later by `Node::updateBound()`.
 2.  **Girders:** For each `Girder`, it calls `o->update_nodes()` (to ensure `ivert`s are up-to-date for bound nodes). Then, `mesh.girder1()` is called to generate the detailed truss geometry (vertices and edges) for the girder. The `pointRange` and `stickRange` of the `Girder` object are updated to reflect its corresponding indices in `mesh.verts` and `mesh.edges`. `mesh.girder1_caps()` adds connections to the end nodes.
 3.  **Ropes:** Similar to girders, `o->update_nodes()` is called, then `mesh.rope()` generates the rope's segments. `pointRange` and `stickRange` are updated.
 4.  **Rings:** `o->update_nodes()` is called, then `mesh.wheel()` generates the ring's geometry. `pointRange` and `stickRange` are updated.
 5.  **Radiators:** `mesh.plateOnGriders()` is used to create the plate geometry, connecting points from the specified girders. `pointRange` and `stickRange` are updated.
 6.  **Welds:** `mesh.bondsBetweenVertRanges()` creates edges between vertices in the specified `pointRange`s of the connected components, within a certain `Rmax`. `pointRange` and `stickRange` are updated.

 **Output:** The `mesh` object (a `Mesh::Builder2` instance) now contains the full geometric representation of the spacecraft as vertices and edges. The `pointRange` and `stickRange` members of the `StructuralComponent`s in the `SpaceCraft` object provide a mapping from the high-level components to their low-level mesh data.

 ### 4.2. `exportSim` Function (Overloaded for `TrussDynamics_d` and `TrussDynamics_f`)

 This function takes the `Mesh::Builder2` object and the `SpaceCraftWorkshop` to populate a `TrussDynamics_d` (or `_f`) simulation object.

 **Key Steps:**
 1.  **Recalculates Neighbors:** Determines the maximum number of neighbors for any vertex in the mesh to allocate `sim.neighs` arrays.
 2.  **Initializes Points:** Copies vertex positions from `mesh.verts` to `sim.points`. Initializes masses (`sim.points[i].w`) to zero.
 3.  **Fills Neighbors and Bonds:** Iterates through `mesh.edges`. For each edge:
    *   Populates `sim.neighs` (neighbor indices) for both endpoints.
    *   Populates `sim.neighBs` (neighbor bond indices) and `sim.neighB2s` (neighbor bond direction) if available.
    *   Calculates the rest length (`l0`) of the bond based on the initial distance between vertices and the `preStrain` of the `StickMaterial`.
    *   Distributes half the bond's mass (`l0 * linearDensity`) to each endpoint's `sim.points[i].w`.
    *   Calculates `param` (l0, Kpush/l0, Kpull/l0, damping) for the bond based on `StickMaterial` properties.
    *   Populates `sim.bonds`, `sim.bparams`, and `sim.maxStrain` arrays for direct bond access.
 4.  **Initializes Simulation State:** Calls `sim.cleanForce()` and `sim.cleanVel()`.

 **Output:** The `sim` object (a `TrussDynamics_d` or `_f` instance) is now ready for physics calculations, containing point masses, initial positions, and bond properties.

 ## 5. Anchoring and Sliders: Connecting High-Level Design to Low-Level Simulation

 This is a critical and complex part of the system, enabling dynamic connections between spacecraft components.

 ### 5.1. `Node` Binding

 A `Node` object can be "bound" to a `StructuralComponent` (e.g., a `Girder`, `Ring`, or `Rope`). This means its position is not fixed in space but is determined by its `boundTo` component and a parameter `calong` (a normalized position along the component, from 0.0 to 1.0).

 *   **`Node::updateBound(Vec3d p0)`:**
     *   Called during `BuildCraft_truss` and `EditSpaceCraft` to update a `Node`'s `pos` and `ivert` (vertex index in the mesh).
     *   It uses `boundTo->nearSide(p0)` to determine which "side" of the structural component the node is attached to (e.g., for a girder, which of its 4 corners in a cross-section).
     *   It then calls `boundTo->pointAlong(calong, along.y, &pos)` to calculate the exact 3D position (`pos`) and the corresponding vertex index (`along.x`) within the `StructuralComponent`'s mesh representation.
     *   Finally, `Node::update_vert()` sets the `Node`'s `ivert` to the global vertex index in `Mesh::Builder2`.

 ### 5.2. `Slider` Component

 A `Slider` is a specialized `Node` that can move along a predefined `Path` of vertices within the `TrussDynamics` simulation.

 *   **`Slider::path`**: A `Path` object (defined in `SpaceCraftComponents.h`) which is essentially an ordered list of vertex indices (`ps`) that the slider can traverse. It also stores `cur` (current position along the path, a float) and `closed` (if the path is a loop).
 *   **`StructuralComponent::updateSlidersPaths(...)`:**
     *   This method is called to initialize or update the `Slider::path` based on the mesh.
     *   It uses `o->sideToPath(side, path.ps)` to populate the `Slider::path.ps` array with the vertex indices along a specific "side" of the `StructuralComponent`'s mesh. For example, a `Girder` has 4 sides, and a `Ring` also has 4 sides (representing its perimeter).
     *   If `Quat4d* ps` (the `TrussDynamics` points array) is provided, it also updates the `Slider::path.cur` by finding the nearest point on the path to the slider's current `pos`, and then updates the slider's `pos` to snap it exactly onto the path.

 ### 5.3. Connecting Sliders to `TrussDynamics`

 The connection between `Slider`s and the `TrussDynamics` simulation is managed by `EdgeVertBond` structures and specific functions in `SpaceCraft2Mesh2.h`.

 *   **`EdgeVertBond` Structure (`TrussDynamics_d.h`)**:
     *   `Vec3i verts`: Stores three vertex indices: `verts.x` and `verts.y` are the two endpoints of the current segment on the slider's path, and `verts.z` is the index of the slider's own vertex in the `TrussDynamics` point array.
     *   `double c`: The interpolation parameter (0.0 to 1.0) indicating the slider's position along the segment defined by `verts.x` and `verts.y`.
     *   `double K`: Stiffness constant for the constraint.
     *   `Vec3d f`: Force applied by the constraint.

 *   **`sliders2edgeverts(SpaceCraft& craft, TrussDynamics_d& sim)` (`SpaceCraft2Mesh2.h`)**:
     *   Called once during simulation initialization (e.g., in `spaceCraftSimulator.h`).
     *   Iterates through all `craft.sliders`.
     *   For each `Slider`, it creates an `EdgeVertBond` in `sim.edgeVertBonds`.
     *   `ev.c = o->path.fromCur( ev.verts.x, ev.verts.y );`: This is crucial. It determines the current segment (`ev.verts.x`, `ev.verts.y`) on the slider's `path` based on its `path.cur` value, and calculates the interpolation parameter `ev.c`.
     *   `ev.verts.z = o->ivert;`: The slider's own vertex (`o->ivert`, which was set during `BuildCraft_truss`) is assigned to `verts.z`.
     *   `ev.K = ...`: Sets the stiffness of the constraint.

 *   **`applySliders2sim(...)` (`SpaceCraft2Mesh2.h`)**:
     *   Called every simulation step (e.g., in `spaceCraftDynamics.cpp`).
     *   Iterates through each `EdgeVertBond` in `sim.edgeVertBonds`.
     *   It retrieves the current positions (`sim.points`) and velocities (`sim.vel`) of the three involved vertices (`ev.verts.x`, `ev.verts.y`, `ev.verts.z`).
     *   It calculates the force (`ev.f`) required to constrain `ev.verts.z` to the line segment defined by `ev.verts.x` and `ev.verts.y` at the position `ev.c`. This force is then applied to the respective points in `sim.forces`.
     *   `o->move(dt, l, v, f)`: The `Slider`'s internal `move` method is called, which updates `o->path.cur` based on `control_speed` and the forces/velocities along the path.
     *   `ev.c = o->path.fromCur( ev.verts.x, ev.verts.y );`: After the slider's `path.cur` is updated, this line recalculates the new segment (`ev.verts.x`, `ev.verts.y`) and interpolation parameter `ev.c` for the next simulation step.

 ## 6. Conclusion

 The spacecraft building system integrates high-level design (Lua scripts defining `SpaceCraft` components) with low-level physics simulation (`TrussDynamics`) through a multi-stage conversion process. The `Mesh::Builder2` acts as an intermediate geometric representation. Complex interactions like `Slider`s are handled by mapping their logical paths to sequences of vertices in the mesh and then applying constraints in the physics engine via `EdgeVertBond`s, ensuring dynamic and realistic behavior. This modular approach allows for flexible design and efficient simulation of intricate spacecraft structures.