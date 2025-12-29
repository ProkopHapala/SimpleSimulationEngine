# SimpleSimulationEngine 2D Demos

This document provides an overview of the various 2D demo programs (`test_*`) available in the `sketches_SDL/2D` directory. Each demo showcases a specific feature, algorithm, or library within the Simple Simulation Engine.

### `test_AppSDL2OGL`
This is a minimal boilerplate application demonstrating the basic setup of an `AppSDL2OGL` 2D window. It allows for simple camera navigation (pan, zoom) using WASD and arrow keys, and draws a crosshair, serving as a foundational starting point for other 2D demos.

### `test_HashMap2D`
This demo visualizes the functionality of `HashMap2D.h`, a 2D spatial hashing data structure. It populates a grid with numerous 2D circular bodies and displays the grid cells, highlighting bodies within the cell currently under the mouse cursor, demonstrating efficient spatial partitioning for object lookup.

### `test_SphereTree2D`
This program demonstrates `SphereTree2D.h`, a 2D hierarchical bounding sphere tree. It populates the scene with many circles and visualizes the tree structure, allowing the user to explore how objects are grouped spatially and to identify the tree node containing the mouse cursor.

### `test_TileTree2D`
This demo showcases `TileTree2D.h`, a 2D spatial data structure optimized for managing and querying rectangular tiles. It displays a grid of tiles, highlighting the tile currently under the mouse cursor, illustrating its use for efficient spatial indexing of rectangular regions.

### `test_Clustering2D`
This program demonstrates a 2D clustering algorithm, likely K-means or a similar iterative method, using `Clustering2D.h`. It generates random points and iteratively refines their cluster assignments, allowing the user to observe the clustering process step-by-step by pressing the `SPACE` key.

### `test_HashMap2D_uniformity`
This demo evaluates the uniformity of point distribution within a `HashMap2D` grid. It populates the hash map with points and visualizes the density or count of points per cell, providing insight into the effectiveness of the hashing function for even spatial distribution.

### `test_HashMap2D_3`
This program is a performance test for `HashMap2D.h`, focusing on the efficiency of adding a large number of points to the 2D hash map. It visualizes the populated grid and measures the computational time required for the insertion process.

### `test_NBodyColHashMap`
This demo simulates 2D N-body collisions, leveraging `HashMap2D.h` for optimized broad-phase collision detection. It visualizes the dynamic interaction of numerous circular bodies, where the hash map efficiently identifies potential collision pairs, and the simulation resolves their physical responses. The simulation can be toggled with the `SPACE` key, and bodies can be dragged with the mouse.

### `test_NBodyWorld`
This program simulates a 2D N-body system using `NBodyWorld2D.h`, where particles interact under forces like gravity. It visualizes the trajectories and dynamic behavior of the bodies, allowing the user to pause/resume the simulation with the `SPACE` key and interactively manipulate individual particles with the mouse.

### `test_PolyLine`
This demo showcases the `PolyLine.h` class, which represents a sequence of connected line segments. It visualizes the polyline and its control points, allowing for interactive manipulation of these points (add, remove, drag with mouse) to demonstrate how the polyline can be dynamically reshaped.

### `test_Voronoi`
This program generates and visualizes 2D Voronoi diagrams using `Voronoi.h`. It displays the cells associated with a set of seed points (sites), allowing the user to add new sites interactively with mouse clicks and observe how the diagram dynamically updates.

### `test_Voronoi2`
This demo provides an enhanced visualization of 2D Voronoi diagrams, building upon the `Voronoi.h` library. It likely offers more detailed rendering options, such as coloring individual cells, displaying centroids, or overlaying the dual Delaunay triangulation, to better illustrate the properties of Voronoi tessellations.

### `test_BranchFract`
This program generates and visualizes a 2D fractal branching pattern, such as a fractal tree. It uses recursive drawing functions, allowing the user to interactively modify parameters like branch angle, length, and recursion depth using keyboard controls to explore different fractal forms.

### `test_Mech2D`
This demo simulates a 2D mechanical system using `Mech2D.h`, which likely involves particles and forces like springs or gravity. It visualizes the dynamic behavior of the system, allowing the user to pause/resume the simulation with the `SPACE` key and interactively apply forces by dragging particles with the mouse.

### `test_MechEuler2D`
This program demonstrates a 2D mechanical simulation using `Mech2D.h`, specifically highlighting the Euler integration method for time-stepping. It visualizes the system's dynamics, allowing comparison of the Euler method's stability and accuracy with other integration schemes.

### `test_AutoMesh2D`
This demo showcases `AutoMesh2D.h`, a library for automatic 2D mesh generation. It allows the user to interactively define input points with mouse clicks, then visualizes the resulting mesh, demonstrating the algorithm's ability to create a tessellation of the 2D space. The mesh can be regenerated by pressing the `SPACE` key.

### `test_SuperSonic2D`
This program simulates 2D supersonic fluid flow, likely using a numerical method to solve the governing equations. It visualizes the complex flow patterns, including shock waves and pressure distributions around objects, demonstrating principles of aerodynamics.

### `test_CommodityNetwork`
This demo simulates a 2D commodity network using `CommodityNetwork.h`. It visualizes the flow of resources or goods between interconnected nodes (e.g., cities, factories) and along edges (e.g., roads, pipelines), demonstrating economic or logistical network dynamics.

### `test_SimplexGrid`
This program demonstrates `SimplexGrid.h`, a 2D data structure based on a triangular grid. It visualizes the grid and allows for the display of scalar data (e.g., height, color) interpolated across its triangular cells, useful for terrain representation or scientific visualization.

### `test_Fluid2D`
This demo simulates 2D fluid dynamics using `Fluid2D.h`. It visualizes the behavior of a fluid on a grid, showing properties like velocity, density, or pressure, and allows for interactive manipulation (e.g., introducing disturbances with the mouse) to observe fluid flow phenomena. The simulation can be toggled with the `SPACE` key.

### `test_TerrainHydraulics`
This program simulates hydraulic erosion on a 2D terrain, using `TerrainHydraulics.h` to model water flow and sediment transport. It visualizes the dynamic evolution of the landscape as water carves channels and deposits material, demonstrating procedural terrain generation. The simulation can be toggled with the `SPACE` key.

### `test_AnalyticalMushroomVortex`
This demo visualizes an analytical solution for a mushroom vortex, a specific fluid dynamics pattern. It displays the velocity field or streamlines of the vortex, allowing the user to explore the mathematical properties of this fluid flow.

### `test_TerrainCubic`
This program demonstrates `TerrainCubic.h`, a method for representing and rendering 2D terrain using cubic interpolation or a cubic grid structure. It visualizes the terrain, potentially leveraging `TiledView.h` for efficient rendering of large landscapes by dividing them into manageable tiles.

### `test_TerrainRBF`
This demo showcases `TerrainRBF.h`, a technique for generating or interpolating 2D terrain using Radial Basis Functions. It visualizes the smoothly interpolated terrain surface, allowing the user to interactively define or adjust the RBF control points that shape the landscape.

### `test_CityGen`
This program demonstrates procedural 2D city generation. It creates a city layout, including road networks, building plots, and potentially simple building representations, allowing the user to explore different urban patterns generated algorithmically.

### `test_Tris2Rect`
This demo explores the conversion between triangular and rectangular data representations in 2D, using `Tris2Rect.h`. It visualizes how data from a triangular mesh can be mapped onto a regular rectangular grid, or vice-versa, which is useful for various data processing tasks.

### `test_GlobOpt2D`
This program demonstrates 2D global optimization for molecular systems using `DynamicOpt.h` and `MoleculeWorld2D.h`. It visualizes the iterative search for minimum energy configurations of molecules, allowing the user to observe the system settling into stable arrangements. The optimization process can be toggled with the `SPACE` key, and the system can be perturbed with the mouse.

### `test_Plotting2D`
This demo showcases the `Plot2D.h` library, a utility for rendering 2D graphs and data visualizations. It displays multiple data series with customizable styles, axes, and labels, demonstrating the library's capabilities for plotting scientific or analytical data.

### `test_Integration1D`
This program demonstrates various numerical integration methods for 1D functions. It visualizes the function and the approximation process (e.g., Riemann sums, trapezoidal rule), allowing the user to compare the accuracy and behavior of different integration techniques.

### `test_ConvexApprox1D`
This demo illustrates the concept of convex approximation for 1D functions. It visualizes a given function and its piecewise linear convex hull or approximation, allowing the user to observe how the approximation quality changes with various parameters.

### `test_PixelGlyphs`
This program demonstrates pixel-based glyph rendering for displaying text. It loads a font texture and renders characters on screen, showcasing basic text rendering capabilities and potentially exploring effects like anti-aliasing or color variations.