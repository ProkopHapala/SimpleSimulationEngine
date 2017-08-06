# Simple Simulation Engine
a minimalistic engine for: 
- Physical simulations 
- Numerical math 
- Game development 
- Computer graphics
- Educational purposes

## Structure
* Computational core is written in C++11. 
* Standalone apps, demos and games use SDL2 and OpenGL 1.2 and 4.0 for Graphical user interface.  
* Some apps use embeded Lua3.2 for scripting. 
* Python 2.7 + numpy + matplotlib is used to make convenient science/engineering packeges partialy wrappaing the C/C++ core libraries. 
* In addition there are some documents using iPython (Jupyter) and javascript+WebGL for rendering interactive 3D graphics and editors of glsl shaders on web.

## Philosophy
Main motivation is to assemble various common algorithma from physics, numerical mathematics, computer graphics and computer science and build minimalistic tools for rapid development of physical simulation programs and games. 

**I should be:**
- **Didactic** - user is expected to be able to understand how things works inside ( i.e. it should not be a blackbox )
- Focused on development of **small programs** ( like sketches/demos ) 
- Each particular part of the engine should be **illustrated** by simplistic example program
- **Minimal dependencies** on 3rd party libraries.
- **Minimalized inter-dependencies** between different modules within the engine. Many modules are designed to work independently, and be easily plugged into other projects 
  * often just as header files (`.h`)
- **No compicated instalation** procedures required
- It is rather a **set of tools** and pre-programed **code samples** rather than enclosed package
- It does not follow rigorous "[encapsulation](https://en.wikipedia.org/wiki/Encapsulation_(computer_programming))" scheme of Object Oriented Programing and other rigorous software enginering methods typical for development of large projects. 
  - e.g. all properties of objects are `public`.
- C++11 is used in pragmatic way, not *ideologically*. 
  - E.g. pass by reference, `templates` and `lambda` expressions are used for programming convenience and computational performance. 
  - At the same time old-schoop C is often prefered: 
    - `printf` is used instead of `std:iostrem` 
    - `char*` instead of `string` 
    - raw pointers `double*` for arrays instead of `std:vector` or [smart pointers](https://msdn.microsoft.com/en-us/library/hh279674.aspx)  
- User is expected and ecouraged to modify the code (including core parts) to match requirements of his particular project. It should be rather starting point of development rather than final product.

#### Dependecies 
- The code aims to have as little dependecies as possible.
- The computational / simulation core should have no dependencies at all.
- However, part of the code is dedicated to visualization and user interface using [simple direct media layer 2](https://www.libsdl.org/) ( SDL2 ) and OpenGL
- other dependencies may be added in examples of particular use cases ( i.e. in praticular games and simulation programs )

# Content
<font color="green"> This is work in progress. Often particular modules of the engine are functional, but they are not put together to build an unified system and work together. </font>

## C++ part
- [Math](cpp/common/math)
  - [2d vector](cpp/common/math/Vec2.h), [3d vector](cpp/common/math/Vec3.h), [3x3 matrix](cpp/common/math/Mat3.h), [4x4 matrix](cpp/common/math/Mat4.h), [Quaternions](cpp/common/math/Mat3.h)
  - [Fast math rutines](cpp/common/math/fastmath.h) with approximations of some functions ( e.g. [goniometry](cpp/common/math/gonioApprox.h), Error function ... )
  - Functions with derivatives ( like [tresholds and sigmoides](cpp/common/math/fastmath.h), and Lorenzieans [etc.](cpp/common/math/functions.h) )
  - [Composing](cpp/common/math/warpFunction2D.h) of complicated 1D,2D,3D functions from simple primitives with analytic gradients 
  - Splines in 1D, 2D 3D, e.g. cubic [hermite spline](cpp/common/math/spline_hermite.h) <font color="green" size=2> (other splines needs to be systematized) </font>
  - [Simplex noise] (cpp/common/math/Noise.h)
  - Basic [linear algebra](cpp/common/math/Lingebra.h) ( matrix multiplication, transpose, GaussJordan solver, Biconjugate gradinet iterative solver, Fitting, Jacobi matrix diagonalization), and [N-dimensional vectors](cpp/common/math/VecN.h) 
  - Minimalistic [Fourier transform](cpp/common/math), with FFT for $2^N$ elemetns
- Math solvers
  - [Brent line search](cpp/common/optimization/lineSearch.h)
  - [Stochastic optimization](cpp/common/optimization/optimizer_random.h)  
  - iterative non-linear [gradient based optimization](cpp/common/math/DynamicOpt.cpp) using [FIRE](http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf) (faster and more robust CG, comparable and sometimes better than BFGS) 
  - Variable step [ODE integrator](cpp/common/dynamics/ODEintegrator.h) by [Runge–Kutta–Fehlberg](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method)
- Computational Geometry  
  - Basic [2d](cpp/common/math/geom2D.h) and [3d](cpp/common/math/geom3D.h) shapes [computational geometry](), [ray tracing](cpp/common/math/raytrace.h), [Covex polygons](cpp/common/math/Convex2d.h), [polygon Mesh](cpp/common/math/Mesh.h) and [regular polyhedra](cpp/common/math/Solids.h)
- [Spatial datastructures](cpp/common/maps) - <font color="green" size=2> This needs cleanup, there is too many duplicite methods/versions. Especially, Rulers should be splited from data. </font>
  - Various methos for acceleration of local interactions and nearest neighbor search e.g. [HashMap](cpp/common/dataStructures/HashMat.h), nD BoundigSphere hierarchy, and grids.
  - Various grids ([trinagular](cpp/common/dataStructures/SimplexGrid.h), [rectangular](cpp/common/dataStructures/Ruler2DFast.h), [cubic](cpp/common/maps/CubicRuler.h), and [triclinc](cpp/common/dataStructures/Grid.h) )
  - Some maps suport raymarching/rasterization e.g. [trinagular](cpp/common/dataStructures/SimplexGrid.h) and [cubic](cpp/common/maps/CubicRuler.h)
- Physical simulations
  - Rigid-Body dynamics in [2D](cpp/common/dynamics/Body2D.h) and [3D](cpp/common/dynamics/Body.h) - with application in [RigidBody molecular dynamics](cpp/common/dynamics/MolecularWorld.h), [Flight dynamics](cpp/common/dynamics/AeroCraft.h) and [vehicle on tarrain](/cpp/apps/Tanks), and [Sailing](cpp/apps/SailWar)
  - [SoftBody](cpp/common/dynamics/SoftBody.h) dynamics with [Truss simulation](cpp/apps/BlockHouseTactics) (e.g. for something like BridgeBuilder)
  - [Molecular mechanics](cpp/common/dynamics/MMFF.h) used [there](cpp/apps/MolecularEditor2)
  - simulation of [hydraulic errosion](cpp/common/maps/TerrainHydraulics.h) for [terrain generation](cpp/sketches_SDL/2D/test_TerrainHydraulics.cpp)
- GUI and interfaces
  - basic SDL2 application prefabricates in [2D](cpp/common_SDL/SDL2OGL/AppSDL2OGL.h) and [3D](cpp/common_SDL/SDL2OGL/AppSDL2OGL_3D.cpp)
  - Basic lightweight [GUI](cpp/common_SDL/SDL2OGL/GUI.h) components (buttons, sliders, input box ... ) for SDL2
  - basic [plotting](cpp/common_SDL/SDL2OGL/Plot2D.h) utilities for 1D,2D and 3D data (inspired by matplolib)
- Graphics
  - OpenGL 1.2/2
    - drawing many basic primitives and more complex shapes in [2D](cpp/common_SDL/SDL2OGL/Draw2D.h) and [3D](cpp/common_SDL/SDL2OGL/Draw3D.h)
  - OpenGL 3+
    - Examples of Rendering [meshes](cpp/sketches_SDL/test_MeshOGL3.cpp), [rendering to texture](test_RenderToTexture.cpp), [instances](cpp/sketches_SDL/test_Instances.cpp), [bilboards](cpp/sketches_SDL/cpp/sketches_SDL/test_Sprites.cpp)
    - [ScreenSpace ambient occlusion](cpp/sketches_SDL/OGL3/test_SSAO.cpp)
    - Ray-Traced analytic primitives e.g. Spherical [atoms](cpp/sketches_SDL/test_Atoms.cpp)
    - Ray-Marching of analytical functions e.g. [molecular orbitals](cpp/sketches_SDL/test_OrbitalRayMarch.cpp).
    - <font color="green" size=2>TODO: (in-shader Constructive solid geometry CGS)</font>
    - [Terrain rendering](cpp/sketches_SDL/test_LandScape.cpp) from heightmap texture with logaritmic depth buffer for large scenes
    - Angularly sensitive bilboards (impostors) for rendering [Vegetation](cpp/sketches_SDL/test_Vegetation.cpp) and clouds

- Apps and Games
  - Molecular editors with [Rigid Body Molecular dynamics](cpp/apps/MolecularEditor) (coarse grained=>faster, no need for intramoleculer forcefiel) with [python interface](python/pyMolecular), and normal [Soft Body Molecular mechanics](cpp/apps/MolecularEditor2)
  - 3D subsonic [aircraft combat simulator](cpp/apps/AeroCombat)
  - 2D [Sail-ships simulator](cpp/apps/SailWar) with [multiplayer](cpp/apps/SailWar_Multi) and [battleship simulator]() 
  - [Land-Combat tactical simulator](cpp/apps/LandTactics) simulator for WWII and modern era 
  - 3D [Tank simulator](cpp/apps/Tanks) 
  - 2D [historical battle simulator](cpp/apps/FormationTactics) (TotalWar like) with large army (~30000 soldiers)
  - Realistic [Spaceship desing and simulation](cpp/apps/OrbitalWar) game 
  
##### C++ Planned / Preliminary
  - multipole expansion
  - fast multipole method or FFT for long range interactions
  - Global optimization algorithms for molecular simulation
  - easy to use generic rendering and physical engine for games (with scene graph) which goes together seamlessly 

## Python
 * [pyMolecular](python/pyMolecular) - Rigid Body Molecular dynamics ([test](test_Molecular.py))
 * [pyRay](python/pyRay) - real-time Constructive Solid Geometry by ray-marching on GPU using pyOpenCL
 * [pySimE](python/pySimE) - Physics and chemistry suite with many useful functions
    * Pure python
      * [enumeration of chemical reactions](test_chem_fuels2table.py) and calculation of [entalpy bilance](test_chem_entalpy.py)
      * Comparison of performance for various [spacecraft propulsion](test_KosmoSuite_shipAccel.py)
    * Using C++ libraries
      * [Trajectory of electron](test_KosmoSuite_elmag.py) in electromagnetic field generated by arbitrary shape of coil 
      * phenomenological calculation of [nuclar bomb yield](test_KosmoSuite_FissionPulse.py) by time propagation of rate equations in 0D.
      * N-body simulations of [planet orbits](test_KosmoSuite_nbody.py) using [RKF45](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method)
      * Simulation of [Rocket launch](test_KosmoSuite_SpaceLaunchODE.py) from plant surface (considering also aerodynamic drag)
  * [pyShock](python/pyShock) time propagation shockwaves in 1D (planar, cylinrical or spherical symmetry) - e.g. for implosion of nuclear bombs (WIP)

## Web
projects using HTML + Javascritp + WebGL + THREE.js and/or Python Jupyter.nb

* [SpaceCombat](projects/SpaceCombat/readme.md) - ideas about combat in solar system supported by physical calculations in python (Jupyter.nb) ([see exampl](https://nbviewer.jupyter.org/github/ProkopHapala/SimpleSimulationEngine/blob/master/projects/SpaceCombat/ch1_basics.ipynb) and [3D graphics](http://htmlpreview.github.io/?https://github.com/ProkopHapala/SimpleSimulationEngine/blob/master/projects/SpaceCombat/HTML/StickSpaceCraft.html) (WebGL+Three.js) - it is strongly inspired by [Project Rho: Atomic rockets ](http://www.projectrho.com/public_html/rocket/)

* Web interfaces for building fragment shaders inspired by [ShaderToy](https://www.shadertoy.com/)
  * [RayTraced analytical primitives](js/GLSL_solid_modeling/ListOfPrimitives.html) for other shaders 
  * [PlanetDesigner](js/PlanetDesigner) - web interface for building glsl fragment shaders for procedural generation and rendering of planets seen from space
  * [ShaderDebug](js/PlanetDesigner) - debuging [domain wraping functions](http://iquilezles.org/www/articles/warp/warp.htm) used for procudural generation of other shaders 