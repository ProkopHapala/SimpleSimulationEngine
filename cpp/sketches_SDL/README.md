# cpp/sketches_SDL — Visual Demos & Test Sketches

Interactive SDL2/OpenGL applications demonstrating engine features. Each `test_*.cpp` is a standalone visual demo — not an automated test. Build via CMake (`WITH_SDL=ON`), run via `tests_bash/sketches_SDL/` shell scripts.

## 2D — `2D/`

- **test_AppSDL2OGL** — minimal boilerplate 2D SDL application
- **test_HashMap2D** — 2D spatial hashing visualization
- **test_HashMap2D_3** — performance test for adding points to HashMap2D
- **test_HashMap2D_uniformity** — evaluates uniformity of point distribution in HashMap2D
- **test_SphereTree2D** — 2D hierarchical bounding sphere tree
- **test_TileTree2D** — 2D spatial data structure for rectangular tiles
- **test_Clustering2D** — 2D clustering algorithm demonstration
- **test_NBodyColHashMap** — 2D N-body collision simulation using HashMap2D
- **test_NBodyWorld** — 2D N-body gravitational system simulation
- **test_PolyLine** — PolyLine class for connected line segments
- **test_Voronoi** — generates and visualizes 2D Voronoi diagrams
- **test_Voronoi2** — enhanced visualization of 2D Voronoi diagrams
- **test_BranchFract** — generates and visualizes a 2D fractal branching pattern
- **test_Mech2D** — simulates a 2D mechanical system
- **test_MechEuler2D** — 2D mechanical simulation highlighting Euler integration
- **test_AutoMesh2D** — demonstrates automatic 2D mesh generation
- **test_SuperSonic2D** — simulates 2D supersonic fluid flow (uses Lingebra)
- **test_CommodityNetwork** — simulates a 2D commodity transport network
- **test_SimplexGrid** — demonstrates a 2D data structure based on triangular grid
- **test_Fluid2D** — simulates 2D fluid dynamics (uses Fluid2D library)
- **test_TerrainHydraulics** — simulates hydraulic erosion on a 2D terrain (uses TerrainGrid2D, Noise)
- **test_AnalyticalMushroomVortex** — visualizes an analytical solution for a mushroom vortex
- **test_TerrainCubic** — 2D terrain rendering using cubic interpolation (uses TerrainCubic, TiledView)
- **test_TerrainRBF** — generates 2D terrain using Radial Basis Functions (uses TerrainRBF, TiledView)
- **test_CityGen** — procedural 2D city generation (uses TerrainRBF, TiledView)
- **test_Tris2Rect** — explores conversion between triangular and rectangular data (uses TiledView)
- **test_GlobOpt2D** — 2D global optimization for molecular systems (uses MoleculeWorld2D, DynamicOpt, Body2D)
- **test_Plotting2D** — showcases the Plot2D library for 2D graphs
- **test_Integration1D** — numerical integration methods for 1D functions
- **test_ConvexApprox1D** — convex approximation for 1D functions
- **test_PixelGlyphs** — pixel-based glyph rendering for displaying text

## 3D — `3D/`

- **test_QuatRotSampling** — visualize orientation of surface directions (u,v) on icosahedron using Quat4f.h and RotationMesh.h
- **test_Projection** — debugging OpenCL projection frustum by point-clouds
- **test_Solids** — demo of Solids.h (Platonic solids) and CMesh.h (simple mesh class)
- **test_SceneGraph** — showcase SceneGraph.h, especially Scene::Group
- **test_SphereSampling** — generate heightmap on planet/asteroid using SphereSampling.h, Noise.h, DrawSphereMap.h (supports `-testRand`, `-octmap` args)
- **test_Patches** — render C1,C2-continuous triangle patches with heightmap using spline_triC1.h and TerrainSimplex.h
- **test_GUI** — GUI.h demo (2D/3D text, Table, GUITextInput, Plot2D, DropDownList, ScisorBox)
- **test_EditorGizmo** — manipulate selected points in space using EditorGizmo.h
- **test_MusculeEditor** — edit organic muscle-like shapes using loft of spline_hermite.h and EditorGizmo.h
- **test_Camera** — 1st-person camera and crosshair as base for shooting/action combat simulator
- **test_MousePicking** — mouse picking of 3D object using raytrace.h raySphere()
- **test_SphereTree** — 3D diffusion-limited aggregation (DLA) using CubicRuler.h and unordered_multimap grid
- **test_Raytracing** — raytrace.h rayTriangles() to calculate ray-triangle intersection and occlusion
- **test_Radiosity** — Radiosity.h to calculate radiosity of light labyrinth of rectangular corridors
- **test_Scatterer** — scattering of particles on thin surfaces using Scatterer2.h, flux transport through channel network
- **test_RayScattererMMC** — RayScatter.h volumetric scattering of light/particles with Monte Carlo method
- **test_Elasticity** — legacy linearized elasticity solver using Truss.h and SoftBody.h (prefer TrussDynamics_d.h)
- **test_Electromagnetic** — simulation of thin hot plasma with ions in magnetic bottle using PotentialFlow.h and Fourier.h poisson solver
- **test_MultipoleAccel** — multipole acceleration using kBoxes.h HierarchicalKpivot
- **test_Mesh** — various mesh operations using Mesh.h (convex polygon by fromPlanes(), findEdges(), collapseEdge())
- **test_RigidBody** — rigid body dynamics hanging on strings using Body.h RigidBody class
- **test_VortexLattice** — PotentialFlow.h velocity field of vortex lattice for flow around wing (lift-line), numerical integration verification
- **test_Collision** — collision of RigidBody with spline HeightMap using Body.h and Terrain2D.h
- **test_CompressiveParticles** — high velocity impact of compressible fluid on wedge boundary using CompressiveParticles.h
- **test_AABBTree** — axis-aligned bounding boxes using AABBTree3D.h and kBoxes.h
- **test_BlockBuilder** — 3D-hashing box-building system as base for panel-housing game
- **test_TrussBuilder** — create truss by keyboard input [KP_1-7] using TrussBuilder.h, export to SoftBody.h and run simulation [SPACE], [k,l] save/load
- **test_SphereGaussSeidel** — SphereGaussSeidel.h to quickly pack boxes close to each other
- **test_Stick** — dynamics of straight sticks hitting a sphere using MMFF.h

## Molecular — `Molecular/`

- **test_sp3space** — relaxation in sp3 orbital space (uses DynamicOpt)
- **test_spRotations** — testing rotations of sp3 orbitals like Slater-Koster tables, uses Approx::AutoApprox (uses DynamicOpt, Lingebra)
- **test_CLCFSF** — Compact Linear Combination of Finite Spherical Functions (multi-center electron force field) (uses DynamicOpt, Lingebra)
- **test_CLCFGO** — Compact Linear Combination of Floating Gaussian Orbitals (multi-center electron force field) (uses DynamicOpt, Lingebra)
- **test_ConfDynamics** — sampling empty free space in collision potential (uses DynamicOpt)
- **test_Multipoles** — multipole expansion demonstration
- **test_SoftMolecularDynamics** — soft molecular dynamics simulation (uses DynamicOpt)
- **test_eFF** — Electron Force Field simulation (uses DynamicOpt)
- **test_eFF_old** — older version of eFF (uses DynamicOpt)
- **test_eFFMC** — compares eFF (eff) and CLCFGO (ff) force fields (uses DynamicOpt)
- **test_MMFFmini** — minimal MMFF molecular force field (uses DynamicOpt) — *known issue: exits with "nff.pairMask is not sorted"*
- **test_FARFF** — Flexible Atom sp-hybridization Reactive Force Field from FlexibleAtomReactiveFF.h (uses DynamicOpt)
- **test_RARFF** — Rigid Atom Reactive Force Field from RARFF.h (uses DynamicOpt)
- **test_RARFF2** — Rigid Atom Reactive Force Field 2 from RARFF2.h (uses DynamicOpt)
- **test_RARFFarr** — Rigid Atom Reactive Force Field using arrays from RARFFarr.h (uses DynamicOpt)
- **test_RARFF_SR** — large scale simulation of grid-based Bucket accelerated RARFF_SR.h (uses DynamicOpt)
- **test_PBD_LJ_cluster** — Position-Based Dynamics with short-range approximate Lennard-Jones potential and cluster constraints (uses DynamicOpt)
- **test_RspFF** — Flexible Atom sp-hybridization forcefield from FspFFclean.h (uses DynamicOpt)
- **test_FTRFF** — Fixed Type Reactive Force-field from FTRFF.h (uses DynamicOpt) — *not rendering*
- **test_RRFF** — Reactive Reactive Force-field from RRFF.h (uses DynamicOpt) — *not rendering*
- **test_EOFF** — Electron Orbital Force Field from EOFF.h (uses DynamicOpt) — *not rendering*
- **test_ESFF** — Electron Spin Force Field from ESFF.h (uses DynamicOpt) — *not rendering*
- **test_BondAdaptedMesh** — optimal sampling grid around molecule

## OGL3 — `OGL3/` (OpenGL 3+, requires GLEW)

- **test_DiffractShader** — diffraction shader demo
- **test_SphereShader** — raycasted sphere impostor demo
- **test_ShaderDepth** — depth / shadow map testing (uses Noise)
- **test_LandScape** — landscape rendering (uses Noise)
- **test_DrawOGL3** — general DrawOGL3 tests
- **test_Instances** — instanced meshes demo (uses Noise)
- **test_Atoms** — atom impostor spheres rendering
- **test_Vegetation** — vegetation rendering (uses Noise)
- **test_Horizont** — sky/horizon rendering (uses Noise)
- **test_VAOs** — VAO (Vertex Array Object) usage demo
- **test_PatchesOGL3** — patch rendering with OpenGL 3
- **test_StencilTextures** — stencil texture demo (uses Noise)
- **test_VolumetricTexture** — volumetric textures (uses Noise)
- **test_RenderToTexture** — render-to-texture pipeline (uses Noise)
- **test_Sprites** — billboard sprites demo (uses Noise)
- **test_AntiAliasing** — anti-aliasing demo
- **test_MeshOGL3** — mesh rendering with OpenGL 3
- **test_GeometryShader** — geometry shader demo
- **test_Tubes** — tube rendering
- **test_OrbitalRayMarch** — orbital ray-marching demo
- **test_Texture** — basic textured rendering
- **test_SSAO** — Screen Space Ambient Occlusion demo
- **test_ScreenOGL3** — ScreenSDL2OGL3 base class demo
- **meshviewer** — generic mesh viewer using MeshRenderOGL3 (takes `.obj` file as argument)

## Shooter — `Shooter/` (collision detection broad-phase algorithms)

- **test_ShotHit** — shooting and hit detection demo
- **test_SweepAndPrune** — sweep-and-prune broad-phase collision detection
- **test_BoxAndSweep** — box-and-sweep collision detection variant
- **test_GridHash** — grid-based hashing for broad-phase collision detection

## Math — `math/`

- **test_TresholdFunc** — threshold function visualization
- **test_FastMath** — fast math approximation comparison (sin, cos, exp, etc.)
- **test_SchroedingerLeapFrog** — 1D Schrödinger equation integration using leap-frog method
- **test_Convex2d** — 2D convex hull and convex operations (uses Convex2d library)
- **test_opt2d** — 2D function optimization demo
- **test_SpaceFillingCurves** — space-filling curves (Hilbert, Morton, etc.) visualization

## Music — `music/` (requires SDL2_mixer)

- **test_playMP3** — MP3 playback demo
- **test_visualize** — music visualization with frequency spectrum

## Network — `network/` (requires SDL2_net)

- **test_UDPNode_client** — UDP networking client demo
- **test_UDPNode_server** — UDP networking server demo

## TEMP — `TEMP/` (experimental, standalone Makefiles)

- **Schroedinger1D** — 1D Schrödinger equation solver (standalone, not in CMake build)
- **function_optimization** — 2D function optimization experiment (standalone, not in CMake build)
