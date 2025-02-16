# SpaceCraft.h

This file defines the `SpaceCrafting::SpaceCraft` class, which represents a spacecraft composed of various interconnected components. It provides functionalities for managing and interacting with these components, including structural elements, functional modules, and collision detection. The class inherits from `CatalogItem`, suggesting its integration into a larger catalog or inventory system.

## Includes

- `<vector>`: Provides dynamic array capabilities.
- `Vec2.h`: Defines the `Vec2` class for 2D vector operations.
- `Vec3.h`: Defines the `Vec3` class for 3D vector operations.
- `Mat3.h`: Defines the `Mat3` class for 3x3 matrix operations.
- `quaternion.h`: Defines the `Quat4d` and `Quat4f` classes for quaternion operations, used for representing rotations.
- `geom3D.h`: Defines geometric primitives and functions for 3D geometry.
- `SpaceCraftComponents.h`: Defines the base class `ShipComponent` and derived classes representing various spacecraft components (e.g., `Node`, `Rope`, `Girder`, `Radiator`).

---

## Types (classes and structs)

### class `SpaceCraft`

The `SpaceCrafting::SpaceCraft` class represents a spacecraft as a collection of interconnected components. It manages the structure and functionality of the spacecraft, providing methods for adding, removing, and interacting with its constituent parts. The class inherits from `CatalogItem`, indicating that spacecraft instances can be managed within a catalog or inventory system.

**Inheritance**

- CatalogItem

#### properties

##### Configuration
- `SpaceCrafting::bool`: `bPrint` - Flag to enable or disable printing of debug information during object creation and modification.

##### External References
- `SpaceCrafting::workshop`: `SpaceCraftWorkshop*` - Pointer to the `SpaceCraftWorkshop` object, providing access to resources like materials and construction tools.

##### Structural Components
- `SpaceCrafting::nodes`: `std::vector<Node*>` - Vector of `Node` pointers, representing the connection points within the spacecraft structure.
- `SpaceCrafting::ropes`: `std::vector<Rope*>` - Vector of `Rope` pointers, representing cable-like structural elements connecting nodes.
- `SpaceCrafting::girders`: `std::vector<Girder*>` - Vector of `Girder` pointers, representing rigid structural beams connecting nodes.
- `SpaceCrafting::rings`: `std::vector<Ring*>` - Vector of `Ring` pointers, representing circular structural elements.
- `SpaceCrafting::guns`: `std::vector<Gun*>` - Vector of `Gun` pointers, representing weapon systems attached to the spacecraft.
- `SpaceCrafting::sliders`: `std::vector<Slider*>` - Vector of `Slider` pointers, representing components that can move along a defined path.

##### Functional Components
- `SpaceCrafting::radiators`: `std::vector<Radiator*>` - Vector of `Radiator` pointers, representing heat dissipation systems.
- `SpaceCrafting::shields`: `std::vector<Shield*>` - Vector of `Shield` pointers, representing protective shields.
- `SpaceCrafting::tanks`: `std::vector<Tank*>` - Vector of `Tank` pointers, representing storage containers for resources.
- `SpaceCrafting::pipes`: `std::vector<Pipe*>` - Vector of `Pipe` pointers, representing conduits for transporting resources.
- `SpaceCrafting::thrusters`: `std::vector<Thruster*>` - Vector of `Thruster` pointers, representing propulsion systems.
- `SpaceCrafting::balloons`: `std::vector<Balloon*>` - Vector of `Balloon` pointers, representing inflatable structures.
- `SpaceCrafting::rocks`: `std::vector<Rock*>` - Vector of `Rock` pointers, representing asteroid or debris shielding.

##### Construction and Organization
- `SpaceCrafting::build_order`: `std::vector<ShipComponent*>` - Vector of `ShipComponent` pointers, defining the order in which components should be built or assembled.
- `SpaceCrafting::welds`: `std::vector<Weld*>` - Vector of `Weld` pointers, representing connections between components.
- `SpaceCrafting::LODs`: `std::vector<int>` - Vector of integers representing levels of detail for OpenGL rendering.

##### Picking
- `SpaceCrafting::pickedTyp`: `int` - Stores the type of the last picked component, used for identifying the component during interaction.

#### methods

##### Component Access
- `getStructuralComponent`: Retrieves a structural component (Girder, Ring, or Rope) based on its ID and type.
- `getPicked`: Returns the picked ship component based on the `pickedTyp` variable.

##### Component Management
- `clear`: Deallocates and clears all component vectors, effectively resetting the spacecraft.
- `add_Node`: Creates a new `Node` object, adds it to the `nodes` vector, and returns its ID.
- `add_Rope`: Creates a new `Rope` object, adds it to the `ropes` vector, and returns its ID.
- `add_Girder`: Creates a new `Girder` object, adds it to the `girders` vector, and returns its ID.
- `add_Ring`: Creates a new `Ring` object, adds it to the `rings` vector, and returns its ID.
- `add_Radiator`: Creates a new `Radiator` object, adds it to the `radiators` vector, and returns its ID.
- `add_Shield`: Creates a new `Shield` object, adds it to the `shields` vector, and returns its ID.
- `add_Tank`: Creates a new `Tank` object, adds it to the `tanks` vector, and returns its ID.
- `add_Thruster`: Creates a new `Thruster` object, adds it to the `thrusters` vector, and returns its ID.
- `add_Gun`: Creates a new `Gun` object, adds it to the `guns` vector, and returns its ID.

##### Geometric Calculations
- `pointOnGirder`: Calculates a point along a girder based on a parameter `c`.
- `linker2line`: Converts a `NodeLinker` (containing two nodes) to a `Line3d` object.
- `plate2quad`: Converts a `Plate` object to a `Quad3d` object, defining the four corners of the plate.
- `rayPlate`: Calculates the intersection point of a ray with a plate, returning the distance along the ray.
- `rayLinkLine`: Calculates the intersection point of a ray with a line segment (representing a rope or girder), returning the distance along the ray.

##### Printing and Debugging
- `printAll_nodes`: Prints information about all nodes in the spacecraft.
- `printAll_ropes`: Prints information about all ropes in the spacecraft.
- `printAll_girders`: Prints information about all girders in the spacecraft.
- `printAll_rings`: Prints information about all rings in the spacecraft.
- `printAll_guns`: Prints information about all guns in the spacecraft.
- `printAll_sliders`: Prints information about all sliders in the spacecraft.
- `printAll_radiators`: Prints information about all radiators in the spacecraft.
- `printAll_shields`: Prints information about all shields in the spacecraft.
- `printAll_pipes`: Prints information about all pipes in the spacecraft.
- `printAll_balloons`: Prints information about all balloons in the spacecraft.
- `printAll_rocks`: Prints information about all rocks in the spacecraft.
- `checkIntegrity`: Performs checks to ensure the integrity of the spacecraft data, such as verifying node bindings.

##### Picking and Selection
- `pick`: Picks the closest component intersected by a ray, used for mouse picking or other selection mechanisms.
- `getVertAlong`: Returns a vertex index along a ship component (Girder, Ring, Rope) based on a parameter `c`.

##### Material Handling
- `updatePanelMaterials`: Updates the material properties of all panel materials in the workshop.