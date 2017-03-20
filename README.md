## Simple Simulation Engine

a minimalistic engine for: 
- physical simulations 
- numerical math 
- game development 
- computer graphics
- educational purposes

written in C++11 with possible python and lua binding using SDL2 and OpenGL for Graphical user interface.

## Philosophy

Main motivation is to assemble various common algorithma from physics, numerical mathematics, computer graphics and computer science and build minimalistic tool for rapid development of physical simulation programs and games. 

**I should be:**
- It should be **didactic** - user is expected to be able to understand how things works inside ( i.e. it should not be a blackbox )
- It is focused on development of **small programs** ( like sketches ) 
- Each particular part of the engine should be illustrated by simplistic example program
- Minimalized dependence on 3rd party libraries 
- Minimalized inter-dependence between different modules within the engine
- No compicated instalation procedures required
- Many modules are designed to work independently, and be easily plugged into other projects - often just as header files (`.h`)
- It is rather a set of tools and pre-programed code samples rather than enclosed package
- It does not follow rigorous "[encapsulation](https://en.wikipedia.org/wiki/Encapsulation_(computer_programming))" scheme of Object Oriented Programing and other rigorous software enginering methods typical for development of large projects. e.g. all properties of objects are `public`. 
- User is expected and ecouraged to modify the code (including core parts) to match requirements of his particular project. It should be rather starting point of development rather than final product.

## Dependecies 
- The code aims to have as little dependecies as possible.
- The computational / simulation core should have no dependencies at all.
- However, part of the code is dedicated to visualization and user interface using [simple direct media layer 2](https://www.libsdl.org/) ( SDL2 ) and OpenGL
- other dependencies may be added in examples of particular use cases ( i.e. in praticular games and simulation programs )

# What is inside

<font color="green"> This is work in progress. Often particular modules of the engine are functional, but they are not put together to build an unified system and work together. </font>

## What is working

- Math
  - 2d vector, 3d vector, 3x3 matrix and Quaternion math
  - fast approximations of some functions ( sin, cos, tan2, erf ... ) and some special function (like various tresholds and sigmoides)
  - Various splines in 1D, 2D 3D, e.g. cubic hermite spline, and other methods for composing of complicated 1D,2D,3D functions from simple primitives 
  - Basic computational geometry, ray tracing, Covex polygons, polygon Mesh.
  - Basic linear algebra ( matrix multiplication, transpose, GaussJordan solver, Biconjugate gradinet iterative solver, Fitting, Jacobi matrix diagonalization)
- basic algorithms and datastrucutres  
  - Various methos for acceleration of local interactions and nearest neighbor search e.g. HashMap, nD BoundigSphere hierarchy, and grids.
  - Rasterization of various grids (trinagular, regular)
- Physical simulations
  - RigidBody dynamics in 2D and 3D - with application in RigidBody molecular dynamics, flight dynamics and vehicle on tarrain
  - SoftBody dynamics with Truss simulation (e.g. for something like BridgeBuilder)
  - simulation of hydraulic errosion for terrain generation
- Apps
  - Molecular editor with rigid body molecular dynamics (fast local interactions) with python interface
  - 3D aircraft flight dynamics simulation (3D rigid body and 3D airfoil) 
  - 2D ships dynamics with sail simulation (2D rigid body and 2D airfoil) with shooting cannons fighting in (1) age-of-sail combat and (2) wwII battleship
  - 3D Tank in terrain 
  - 2D top view simulation of large army (~30000 soldiers) for total-war like game (interactions accelerated by BoxBuffer)
  
## What is preliminary
  - multipole expansion

## What is planned
  - fast multipole method or FFT for long range interactions
  - Global optimization algorithms for molecular simulation
  - easy to use generic rendering and physical engine for games (with scene graph) which goes together seamlessly 
