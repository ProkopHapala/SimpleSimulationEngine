# SimpleSimulationEngine 3D Demos

This document provides an overview of the various 3D demo programs (`test_*`) available in the `sketches_SDL/3D` directory. Each demo showcases a specific feature, algorithm, or library within the Simple Simulation Engine.

### `test_AABBTree`
This demo visualizes an Axis-Aligned Bounding Box (AABB) tree, a hierarchical spatial data structure used to accelerate collision detection and other spatial queries. It demonstrates the use of `AABBTree3D.h` and `kBoxes.h` to partition a scene with many objects. The visualization likely shows a collection of geometric objects and their enclosing bounding boxes at different levels of the tree, allowing the user to inspect the culling efficiency. The scene can typically be navigated using the mouse to rotate the camera and the scroll wheel to zoom.

### `test_BlockBuilder`
This program is a demonstration of an interactive construction system based on a 3D spatial hash grid. It serves as a prototype for a "panel-housing game" where users can place and remove blocks in a 3D world. The core logic relies on a 3D hashing system for efficient block management. User interaction is likely handled via mouse clicks to add or remove blocks and keyboard shortcuts to select different block types.

### `test_Camera`
This is a foundational demo for any first-person application. It sets up a standard 1st-person camera with a crosshair, suitable for a shooting game or a flight simulator. The user can navigate the 3D space using standard WASD keys for movement and the mouse or arrow keys for changing the camera orientation (pitch/yaw). This test is a minimal example of the `AppSDL2OGL_3D` application framework.

### `test_Collision`
This demo showcases physics simulation, specifically the collision between a 3D rigid body (from `Body.h`) and a deformable terrain. The terrain is represented by a spline-based height map using `Terrain2D.h`, allowing for smooth, continuous surfaces. The highlight is observing the dynamic interaction and realistic collision response between a solid object and a complex, non-planar surface.

### `test_CompressiveParticles`
This program demonstrates a particle-based simulation of a compressible fluid, likely using a method similar to Smoothed Particle Hydrodynamics (SPH). It visualizes a high-velocity impact of the fluid particles against a static wedge-shaped boundary, showcasing effects like shockwaves, compression, and fluid dynamics without a mesh. The simulation is self-running, with camera controls for observation.

### `test_EditorGizmo`
This is a utility demo showcasing the `EditorGizmo.h` library, a standard tool for 3D content creation. It displays a manipulation gizmo (with arrows for translation and rings for rotation) that can be attached to selected objects or points in space. The user can typically select points by clicking the mouse and then drag the gizmo handles to move or rotate the selection, providing a direct and intuitive way to edit 3D scenes.

### `test_Elasticity`
This demo runs a simulation of a deformable structure using a legacy linearized elasticity solver from `Truss.h` and `SoftBody.h`. It shows how a truss or other soft body deforms under external forces. While functional, the description notes that the more modern and capable `TrussDynamics_d.h` is preferred for new projects. Users can likely interact by applying forces with the mouse.

### `test_Electromagnetic`
This is a physics simulation of charged particles in a magnetic field. It models a thin, hot plasma where ions are trapped within a "magnetic bottle". The simulation calculates the aggregate electric and magnetic fields using a Poisson solver (`Fourier.h`) and `PotentialFlow.h`, then integrates the motion of the ions. The description suggests pressing the 'm' key to run or control the simulation.

### `test_GUI`
This program is a comprehensive showcase of the engine's immediate-mode Graphical User Interface library, `GUI.h`. It displays a variety of UI widgets, including 2D and 3D text rendering, data tables, text input fields, 2D plots, dropdown lists, and scissor boxes for clipping UI elements. Users can interact with all the elements to see their functionality, making this a key demo for learning how to build UIs in the engine.

### `test_Mesh`
This demo focuses on fundamental mesh processing operations using the `Mesh.h` library. It loads or generates a 3D mesh and allows the user to perform and visualize various algorithms, such as finding all edges, collapsing a selected edge, or generating a convex polygon from a set of planes. It's a technical demonstration of the engine's core geometry manipulation capabilities.

### `test_MousePicking`
This is a focused demo on object selection in a 3D scene. It uses the `raytrace.h` library to cast a ray from the camera through the mouse cursor to determine which object is being pointed at. The test specifically highlights `raySphere()` for picking spherical objects. When an object is "picked," it is typically highlighted. This is a fundamental mechanic for any interactive 3D application.

### `test_MultipoleAccel`
This demo illustrates an advanced algorithm for accelerating N-body simulations. It uses the Fast Multipole Method (FMM) via `kBoxes.h` and `HierarchicalKpivot` to approximate the influence of distant particles, reducing the computational complexity from O(N^2) to O(N). The visualization would likely show a system of particles and the hierarchical grid used to group them for the multipole expansion.

### `test_MusculeEditor`
This is a creative modeling tool for generating organic, muscle-like shapes. It uses Hermite splines (`spline_hermite.h`) to define cross-sections which are then "lofted" to form a smooth 3D surface. The shape can be interactively edited by manipulating the spline control points using the `EditorGizmo.h` interface, allowing for intuitive, real-time sculpting.

### `test_Patches`
A graphics-focused demo that showcases advanced terrain rendering. It uses `spline_triC1.h` to render a surface made of C1 or C2 continuous curved triangular patches. This avoids the faceted look of standard triangle meshes and creates a much smoother, more realistic landscape from a heightmap, using `TerrainSimplex.h` for the underlying grid structure.

### `test_Projection`
This is a technical debugging tool for developers. It visualizes the viewing frustum of a camera by rendering a cloud of points that lie within it. Its primary purpose is to debug projection matrices and culling logic, especially for custom renderers suchas those using OpenCL.

### `test_QuatRotSampling`
This is a math and geometry visualization that deals with mapping orientations on a sphere. It uses an icosahedron as a base and, with the help of Quaternions (`Quat4f.h`), visualizes a consistent set of local coordinate systems (u,v directions) on the surface. This is a key technique for applying seamless textures or vector fields (like wind patterns) to a sphere without singularities at the poles.

### `test_Radiosity`
This demo implements the classic Radiosity global illumination algorithm using `Radiosity.h`. It calculates the diffuse inter-reflection of light within a scene, in this case, a labyrinth of corridors. The result is a soft, realistic lighting effect with color bleeding, where light bounces off surfaces and illuminates others. The calculation is iterative, and the user can likely watch the scene's lighting converge to the final solution.

### `test_RayScattererMMC`
This program demonstrates volumetric light scattering using a Monte Carlo method, implemented in `RayScatter.h`. It simulates the transport of light or particles through a medium by tracing the random paths of many individual rays as they scatter. This is useful for rendering effects like fog, smoke, or murky water.

### `test_Raytracing`
A fundamental graphics demo that showcases the core concepts of ray tracing from `raytrace.h`. It casts rays from a viewpoint and calculates their intersection with a scene composed of triangles. A key feature is demonstrating occlusion, where rays are blocked by objects. The visualization typically colors rays to indicate whether they hit their target or were occluded, as seen in the `visualizeRayTracing` function in `spaceCraftEditor.cpp`.

### `test_RigidBody`
This is a classic physics simulation of rigid body dynamics using the `Body.h` library. The scene consists of one or more solid objects that are connected by string-like constraints and are subject to forces like gravity. Users can typically interact with the simulation by grabbing and throwing the objects with the mouse.

### `test_Solids`
A simple demo program for displaying and testing the generation of primitive shapes. It uses `Solids.h` to create various polyhedra, including the Platonic solids (cube, icosahedron, etc.), and renders them using a basic mesh class (`CMesh.h`). It's a good starting point for verifying that the 3D rendering pipeline is working correctly.

### `test_Scatterer`
This demo implements a particle or light scattering simulation using the direct iterative approach found in `Scatterer2.h`. Instead of solving a large matrix like in radiosity, this method models the transport of flux through a network of pre-defined channels connecting different surface elements. It's an alternative approach to global illumination that can be more efficient for certain scene types.

### `test_SceneGraph`
This program demonstrates the use of a scene graph (`SceneGraph.h`) to organize objects in a 3D world. It showcases how to create parent-child relationships between objects (e.g., a moon orbiting a planet). When the parent object is transformed (moved or rotated), all of its children are transformed along with it, simplifying the management of complex, articulated scenes.

### `test_SphereGaussSeidel`
This is a physics demo showcasing a method for resolving collisions and packing objects together. It uses an iterative Gauss-Seidel solver (`SphereGaussSeidel.h`) to enforce non-penetration constraints between a collection of spheres or boxes. The user can watch as the objects, which may start in an overlapping state, quickly settle into a valid, tightly-packed configuration.

### `test_SphereSampling`
This demo procedurally generates a planet or asteroid using `SphereSampling.h` and `Noise.h`. It applies one or more layers of a noise function (like Simplex noise) to a sphere to create a heightmap, resulting in features like continents, mountains, and craters. The result is rendered using `DrawSphereMap.h`, and the user can rotate the generated celestial body.

### `test_SphereTree`
This program simulates the growth of a fractal structure using the Diffusion-Limited Aggregation (DLA) algorithm. It starts with a seed particle, and new particles are added one by one, performing a random walk until they collide and stick to the growing cluster. The result is a natural, tree-like or coral-like structure. The simulation is accelerated using a `CubicRuler.h` spatial grid.

### `test_Stick`
This is a physics simulation that uses a Molecular Mechanics Force Field (`MMFF.h`), typically used for atomic-scale simulations, to model the dynamics of macroscopic objects. It demonstrates "stick" molecules colliding with a large sphere, showcasing an interesting application of molecular dynamics principles to a different scale.

### `test_TrussBuilder`
An interactive construction and simulation tool. Using `TrussBuilder.h`, the user can build a 2D or 3D truss structure by adding nodes and beams with the numeric keypad (`KP_1` to `KP_7`). After building, pressing the `SPACE` bar starts a physics simulation (`SoftBody.h` or `TrussDynamics_d.h`), showing how the created structure behaves under stress. The `k` and `l` keys can be used to save and load designs.

### `test_VortexLattice`
This is a computational fluid dynamics (CFD) demo for aerodynamics. It uses the vortex lattice method (`PotentialFlow.h`) to calculate the airflow and lift generated by a wing. The wing is represented as a "lift-line," and the demo visualizes the resulting velocity field around it. It also includes functions for numerical integration to verify the analytical results.