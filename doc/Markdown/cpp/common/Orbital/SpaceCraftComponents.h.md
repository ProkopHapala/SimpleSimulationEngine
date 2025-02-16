# SpaceCraftComponents.h

This file defines the core component classes used in the `SpaceCrafting` namespace for building and managing spacecraft. It includes definitions for structural elements, functional modules, materials, and other related data structures. These components are designed to be interconnected and managed by the `SpaceCraft` class.

## Includes

- `<stdio.h>`: Provides standard input/output functions.
- `<string.h>`: Provides string manipulation functions.
- `<string>`: Provides string class.
- `<vector>`: Provides dynamic array capabilities.
- `<unordered_map>`: Provides hash table implementation for efficient key-value storage.
- `datatypes.h`: Defines fundamental data types used throughout the project.
- `Vec2.h`: Defines the `Vec2` class for 2D vector operations.
- `Vec3.h`: Defines the `Vec3` class for 3D vector operations.
- `Mat3.h`: Defines the `Mat3` class for 3x3 matrix operations.
- `quaternion.h`: Defines the `Quat4d` and `Quat4f` classes for quaternion operations, used for representing rotations.
- `geom3D.h`: Defines geometric primitives and functions for 3D geometry.
- `macroUtils.h`: Defines macros and utility functions.
- `containers.h`: Defines custom container classes.
- `Buckets.h`: Defines the `Buckets` class for spatial partitioning, used for collision detection.

---
## Types (classes and structs)

--
### class `ComponetKind`

This enum class defines the different types of components that can be used to construct a spacecraft. It provides a way to identify and categorize components within the system.

**Inheritance**

- int

--
### class `CatalogItem`

This class serves as a base class for all catalog items, such as materials, commodities, and ship components. It provides basic properties like ID, kind, and name.

#### properties

- `CatalogItem::int`: `id` - Unique identifier for the catalog item.
- `CatalogItem::kind`: `int` - Integer representing the type or category of the catalog item.
- `CatalogItem::name`: `char` - Character array storing the name of the catalog item.

#### methods

- `print`: Prints basic information about the catalog item.

---
### class `Material`

This class represents a material with physical properties such as density, strength, elastic modulus, reflectivity, and melting point.

**Inheritance**

- CatalogItem

#### properties

- `Material::double`: `density` - Density of the material in kg/m³.
- `Material::Spull`: `double` - Tensile strength of the material in Pascals (Pa).
- `Material::Kpull`: `double` - Tensile elastic modulus of the material in Pascals (Pa).
- `Material::reflectivity`: `double` - Reflectivity of the material (dimensionless).
- `Material::Tmelt`: `double` - Melting point of the material in Kelvin (K).
- `Material::Spush`: `double` - Compressive strength of the material in Pascals (Pa).
- `Material::Kpush`: `double` - Compressive elastic modulus of the material in Pascals (Pa).

#### methods

- `print`: Prints detailed information about the material, including its physical properties.

---
### class `Commodity`

This class represents a commodity, such as fuel or water, with properties like density.

**Inheritance**

- CatalogItem

#### properties

- `Commodity::double`: `density` - Density of the commodity in kg/m³.

#### methods

- `print`: Prints information about the commodity, including its density.

---
### class `FuelType`

This class represents a specific type of fuel, inheriting from `Commodity` and adding an energy density property.

**Inheritance**

- Commodity

#### properties

- `FuelType::double`: `EnergyDesity` - Energy density of the fuel in Joules per kilogram (J/kg).

#### methods

- `print`: Prints information about the fuel type, including its energy density.

---
### class `StickMaterial`

This class represents the material properties of a structural stick or beam, including diameter, wall thickness, area, linear density, strength, elastic modulus, reflectivity, and melting point.

**Inheritance**

- CatalogItem

#### properties

- `StickMaterial::int`: `materialId` - Index of the base material in the material catalog.
- `StickMaterial::diameter`: `double` - Diameter of the stick in meters (m).
- `StickMaterial::wallThickness`: `double` - Thickness of the stick's wall in meters (m).
- `StickMaterial::area`: `double` - Cross-sectional area of the stick in square meters (m²).
- `StickMaterial::linearDensity`: `double` - Linear density of the stick in kilograms per meter (kg/m).
- `StickMaterial::reflectivity`: `double` - Reflectivity of the stick (dimensionless).
- `StickMaterial::Tmelt`: `double` - Melting point of the stick in Kelvin (K).
- `StickMaterial::damping`: `double` - Damping coefficient of the stick.
- `StickMaterial::Kpull`: `double` - Tensile elastic modulus of the stick in Newtons (N).
- `StickMaterial::Spull`: `double` - Tensile strength of the stick in Newtons (N).
- `StickMaterial::preStrain`: `double` - Pre-strain of the stick (dimensionless).
- `StickMaterial::Spush`: `double` - Compressive strength of the stick in Newtons (N).
- `StickMaterial::Kpush`: `double` - Compressive elastic modulus of the stick in Newtons (N).

#### methods

- `update`: Updates the stick material properties based on the base material properties.
- `print`: Prints detailed information about the stick material, including its dimensions and physical properties.

---
### class `PanelMaterial`

This class represents the material properties of a panel, including its layers, area density, and supporting stick material.

**Inheritance**

- CatalogItem

#### properties

- `PanelMaterial::std`: `layers` - Vector of `PanelLayer` objects, defining the layers of the panel.
- `PanelMaterial::areaDensity`: `double` - Area density of the panel in kilograms per square meter (kg/m²).
- `PanelMaterial::stickMaterialId`: `int` - Index of the stick material used to support the panel in the material catalog.

#### methods

- `evalAreaDensity`: Calculates the area density of the panel based on its layers.
- `print`: Prints information about the panel material, including its area density and number of layers.

---
### class `ThrusterType`

This class represents the type of a thruster, including its efficiency, exhaust velocity range, fuel type, and propellant.

**Inheritance**

- CatalogItem

#### properties

- `ThrusterType::double`: `efficiency` - Efficiency of the thruster (dimensionless).
- `ThrusterType::veMin`: `double` - Minimum exhaust velocity of the thruster in meters per second (m/s).
- `ThrusterType::veMax`: `double` - Maximum exhaust velocity of the thruster in meters per second (m/s).
- `ThrusterType::exhaustFuel`: `bool` - Flag indicating whether the burned fuel is added to the propellant mass.
- `ThrusterType::fuel`: `FuelType*` - Pointer to the `FuelType` object used by the thruster.
- `ThrusterType::Propelant`: `Commodity*` - Pointer to the `Commodity` object used as propellant by the thruster.

#### methods

- `print`: Prints information about the thruster type, including its efficiency, exhaust velocity range, fuel, and propellant.

---
### class `GunType`

This class represents the type of a gun, including its recoil.

**Inheritance**

- CatalogItem

#### properties

- `GunType::double`: `recoil` - Recoil force of the gun.

#### methods

- `print`: Prints information about the gun type, including its recoil.

---
### class `ShipComponent`

This class serves as a base class for all spacecraft components, providing basic properties like ID, kind, shape, material, mass, and point range.

**Inheritance**

- Object

#### properties

- `ShipComponent::int`: `id` - Unique identifier for the ship component.
- `ShipComponent::face_mat`: `int` - Index of the material used for the component's faces in the material catalog.
- `ShipComponent::mass`: `double` - Mass of the component in kilograms (kg).
- `ShipComponent::pointRange`: `Vec2i` - Range of point indices associated with this component.
- `ShipComponent::shape`: `int` - Integer representing the shape of the component.

#### methods

- `print`: Prints basic information about the ship component.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `StructuralComponent`

This class represents a structural component of the spacecraft, inheriting from `ShipComponent` and adding properties for managing nodes and their relationships.

**Inheritance**

- ShipComponent

#### properties

- `Node::vec4`: `nodes` - Vector of `Node*` pointers, representing the nodes connected to this structural component.
- `Node::mvert`: `int` - Number of vertices per segment.

#### methods

- `rotMat`: Calculates the rotation matrix for the component.
- `nearSide`: Determines the nearest side of the component to a given point.
- `pointAlong`: Calculates a point along the component based on a parameter.
- `sideToPath`: Returns a path (array of vertex indices) along a specified side of the component.
- `update_nodes`: Updates the vertex indices of the component's nodes.
- `updateSlidersPaths`: Updates the paths of any sliders attached to this component.
- `findNearestPoint`: Finds the nearest point on the component to a given point.
- `toBuckets`: Assigns the component's points to buckets for spatial partitioning.
- `print`: Prints basic information about the structural component.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Weld`

This class represents a weld connection between two structural components.

**Inheritance**

- ShipComponent

#### properties

- `Weld::double`: `Rmax` - Maximum distance for the weld connection.
- `Weld::comps`: `vec2<StructuralComponent*>` - Vector of two `StructuralComponent*` pointers, representing the two components connected by the weld.

#### methods

- `print`: Prints information about the weld, including the IDs of the connected components.

---
### class `Node`

This class represents a node, which is a connection point in the spacecraft structure.

**Inheritance**

- Object

#### properties

- `Node::int`: `ivert` - Index of the vertex in the mesh/truss that corresponds to this node.
- `Node::pos`: `Vec3d` - Position of the node in 3D space.
- `Node::boundTo`: `StructuralComponent*` - Pointer to the structural component to which this node is bound (e.g., a girder or ring).
- `Node::calong`: `double` - Position along the bound component, used for positioning the node.
- `Node::length`: `double` - Length of the bound component.
- `Node::vec2`: `along` - Index of the point along the bound component.

#### methods

- `updateBound`: Updates the position of the node based on its binding to a structural component.
- `update_vert`: Updates the vertex index of the node based on its binding.
- `print`: Prints information about the node, including its position and binding information.

---
### class `NodeLinker`

This class represents a component that links two nodes together, such as a rope or girder.

**Inheritance**

- StructuralComponent

#### properties

- `NodeLinker::double`: `length` - Length of the connection between the two nodes.

#### methods

- `rotMat`: Calculates the rotation matrix for the component.
- `nearSide`: Determines the nearest side of the component to a given point.
- `pointAlong`: Calculates a point along the component based on a parameter.
- `print`: Prints basic information about the node linker.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Girder`

This class represents a rigid structural beam connecting two nodes.

**Inheritance**

- NodeLinker

#### properties

- `Girder::int`: `nseg` - Number of segments along the girder.
- `Girder::mseg`: `int` - Number of segments around the girder.
- `Girder::wh`: `Vec2d` - Width and height of the girder's cross-section.
- `Girder::up`: `Vec3d` - Up vector defining the orientation of the girder.
- `Girder::st`: `Quat4i` - Stick types for different edges of the girder.

#### methods

- `rotMat`: Calculates the rotation matrix for the girder.
- `nearSide`: Determines the nearest side of the girder to a given point.
- `pointAlong`: Calculates a point along the girder based on a parameter.
- `sideToPath`: Returns a path (array of vertex indices) along a specified side of the girder.
- `print`: Prints detailed information about the girder, including its dimensions and properties.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Ring`

This class represents a circular structural element.

**Inheritance**

- StructuralComponent

#### properties

- `Ring::BodyPose`: `pose` - Position and orientation of the ring.
- `Ring::nseg`: `int` - Number of segments around the ring.
- `Ring::R`: `double` - Radius of the ring.
- `Ring::wh`: `Vec2d` - Width and height of the ring's cross-section.
- `Ring::st`: `Quat4i` - Stick types for different edges of the ring.

#### methods

- `rotMat`: Calculates the rotation matrix for the ring.
- `nearSide`: Determines the nearest side of the ring to a given point.
- `pointAlong`: Calculates a point along the ring based on a parameter.
- `sideToPath`: Returns a path (array of vertex indices) along a specified side of the ring.
- `print`: Prints detailed information about the ring, including its dimensions and properties.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Rope`

This class represents a flexible cable-like structure connecting two nodes.

**Inheritance**

- NodeLinker

#### properties

- `Rope::double`: `thick` - Thickness of the rope.
- `Rope::nseg`: `int` - Number of segments along the rope.

#### methods

- `rotMat`: Calculates the rotation matrix for the rope.
- `nearSide`: Determines the nearest side of the rope to a given point.
- `pointAlong`: Calculates a point along the rope based on a parameter.
- `sideToPath`: Returns a path (array of vertex indices) along a specified side of the rope.
- `print`: Prints basic information about the rope.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Modul`

This class represents a modular component with a defined volume and pose.

**Inheritance**

- ShipComponent

#### properties

- `Modul::BodyPose`: `pose` - Position and orientation of the module.
- `Modul::bbox`: `Box` - Bounding box of the module.
- `Modul::span`: `Vec3d` - Span of the module.
- `Modul::volume`: `double` - Volume of the module.

#### methods

- `pick`: Performs picking operations on the module.
- `print`: Prints basic information about the module.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Tank`

This class represents a storage tank for resources, inheriting from `Modul` and adding properties for commodity type and fill level.

**Inheritance**

- Modul

#### properties

- `Tank::int`: `commodityId` - Index of the commodity stored in the tank in the commodity catalog.
- `Tank::double`: `filled` - Fill level of the tank (0.0 to 1.0).

#### methods

- `print`: Prints information about the tank, including its capacity and fill level.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Balloon`

This class represents an inflatable structure, inheriting from `Modul`.

**Inheritance**

- Modul

#### methods

- `print`: Prints basic information about the balloon.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Rock`

This class represents an asteroid or debris shield, inheriting from `Modul`.

**Inheritance**

- Modul

#### methods

- `print`: Prints basic information about the rock.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Pipe`

This class represents a conduit for transporting resources between components.

**Inheritance**

- ShipComponent

#### properties

- `Pipe::double`: `maxFlow` - Maximum flow rate through the pipe.
- `Pipe::path`: `Path` - Path along which the pipe runs.
- `Pipe::a`: `ShipComponent*` - Pointer to the starting component of the pipe.
- `Pipe::b`: `ShipComponent*` - Pointer to the ending component of the pipe.

#### methods

- `print`: Prints information about the pipe, including its capacity and connected components.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Plate`

This class represents a flat surface component, such as a radiator or shield.

**Inheritance**

- ShipComponent

#### properties

- `Plate::double`: `area` - Area of the plate.
- `Plate::g1`: `int` - Index of the first supporting girder.
- `Plate::g1span`: `Vec2d` - Span along the first girder.
- `Plate::g2`: `int` - Index of the second supporting girder.
- `Plate::g2span`: `Vec2d` - Span along the second girder.
- `Plate::plate_mat`: `int` - Material of the plate.

#### methods

- `print`: Prints basic information about the plate.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Radiator`

This class represents a heat dissipation panel, inheriting from `Plate` and adding a temperature property.

**Inheritance**

- Plate

#### properties

- `Radiator::double`: `temperature` - Temperature of the radiator.

#### methods

- `print`: Prints information about the radiator, including its temperature.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Shield`

This class represents a protective panel, inheriting from `Plate`.

**Inheritance**

- Plate

#### methods

- `print`: Prints basic information about the shield.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Collector`

This class represents a resource collection panel, inheriting from `Plate`.

**Inheritance**

- Plate

#### methods

- `print`: Prints basic information about the collector.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Thruster`

This class represents a propulsion system, inheriting from `Modul` and adding properties for thrust, power, and consumption.

**Inheritance**

- Modul

#### properties

- `Thruster::int`: `type` - Type of the thruster.
- `Thruster::double`: `thrust` - Thrust force of the thruster.
- `Thruster::power`: `double` - Power consumption of the thruster.
- `Thruster::consumption`: `double` - Consumption rate of the thruster.

#### methods

- `print`: Prints information about the thruster, including its performance characteristics.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Rotor`

This class represents a rotating component, inheriting from `ShipComponent` and adding properties for radius, power, torque, and inertia.

**Inheritance**

- ShipComponent

#### properties

- `Rotor::double`: `Radius` - Radius of the rotor.
- `Rotor::power`: `double` - Power of the rotor.
- `Rotor::torque`: `double` - Torque of the rotor.
- `Rotor::Inertia`: `double` - Moment of inertia of the rotor.

#### methods

- `print`: Prints information about the rotor, including its performance characteristics.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Slider`

This class represents a component that can slide along a path, such as a rail or girder.

**Inheritance**

- Node

#### properties

- `Slider::Path`: `path` - Path along which the slider can move.
- `Slider::double`: `maxDist` - Maximum distance the slider can deflect perpendicularly from the path.
- `Slider::forceMax`: `double` - Maximum force that can be exerted by the slider.
- `Slider::powerMax`: `double` - Maximum power that can be consumed by the slider.
- `Slider::maxSpeed`: `double` - Maximum speed of the slider.
- `Slider::speed`: `double` - Current speed of the slider.
- `Slider::springK`: `double` - Spring constant for the slider's restoring force.
- `Slider::Kdv`: `double` - Damping coefficient for the slider's velocity.
- `Slider::vel`: `double` - Current velocity of the slider.
- `Slider::mass`: `double` - Mass of the slider.
- `Slider::icontrol`: `int` - Index of the degree of freedom that controls the slider.

#### methods

- `move`: Moves the slider along its path based on applied forces and constraints.
- `updatePath`: Updates the path of the slider based on its bound structural component.
- `print`: Prints information about the slider, including its properties and constraints.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Accelerator`

This class represents an accelerator component, such as a gun, that propels particles or projectiles.

**Inheritance**

- ShipComponent

#### properties

- `Accelerator::int`: `suppType` - Type of support structure (e.g., ring or girder).
- `Accelerator::suppId`: `int` - ID of the support structure.
- `Accelerator::suppSpan`: `Vec2d` - Span along the support structure.
- `Accelerator::path`: `Path` - Path along which the projectile travels.
- `Accelerator::lenght`: `double` - Length of the accelerator.
- `Accelerator::PowerPeak`: `double` - Peak power of the accelerator.
- `Accelerator::PulseEnergy`: `double` - Energy per pulse of the accelerator.
- `Accelerator::PulseDuration`: `double` - Duration of each pulse.
- `Accelerator::PulsePeriod`: `double` - Period between pulses.

#### methods

- `print`: Prints information about the accelerator, including its performance characteristics.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `Gun`

This class represents a weapon system, inheriting from `Accelerator` and adding properties for aperture and divergence.

**Inheritance**

- Accelerator

#### properties

- `Gun::int`: `gunId` - Index of the gun type in the catalog.
- `Gun::Aperture`: `double` - Aperture of the gun in square meters (m²).
- `Gun::divergence`: `double` - Divergence angle of the gun's projectile (dimensionless).

#### methods

- `print`: Prints information about the gun, including its performance characteristics and projectile properties.
- `component_kind`: Returns the `ComponetKind` enum value representing the type of the component.

---
### class `SpaceCraftWorkshop`

This class manages the catalogs of materials, commodities, and other resources used in spacecraft construction.

#### properties

- `SpaceCraftWorkshop::bool`: `bPrint` - Flag to enable or disable printing of debug information.
- `SpaceCraftWorkshop::materials`: `Dict<Material>` - Dictionary of `Material` objects, indexed by name.
- `SpaceCraftWorkshop::commodities`: `Dict<Commodity>` - Dictionary of `Commodity` objects, indexed by name.
- `SpaceCraftWorkshop::fuels`: `Dict<FuelType>` - Dictionary of `FuelType` objects, indexed by name.
- `SpaceCraftWorkshop::panelMaterials`: `Dict<PanelMaterial>` - Dictionary of `PanelMaterial` objects, indexed by name.
- `SpaceCraftWorkshop::stickMaterials`: `Dict<StickMaterial>` - Dictionary of `StickMaterial` objects, indexed by name.
- `SpaceCraftWorkshop::thrusterTypes`: `Dict<ThrusterType>` - Dictionary of `ThrusterType` objects, indexed by name.
- `SpaceCraftWorkshop::gunTypes`: `Dict<GunType>` - Dictionary of `GunType` objects, indexed by name.

#### methods

- `add_Material`: Creates a new `Material` object, adds it to the `materials` dictionary, and returns its ID.
- `add_StickMaterial`: Creates a new `StickMaterial` object, adds it to the `stickMaterials` dictionary, and returns its ID.