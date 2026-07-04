# tests_bash — Build & Run Scripts for Visual Demos

Shell scripts that compile (via `make` into `cpp/Build/`) and run each demo with ASan preloaded. Each script symlinks `common_resources` and `data/` as needed. Edit the script to uncomment the desired `name=...` line, then run `./script.sh`.

## apps/ — Full Application Demos

- **AeroCombat.sh** — `AeroCombat_main`: aircraft simulator with realistic aerodynamics, mouse+keyboard control, shooting targets (icosahedra on sticks above terrain). Also has `AeroCraft_editor` / `AeroCraft_editor2` for aircraft model editing.
- **BlockHouseTactics.sh** — `BlockHouseTactics_main`: build blockhouse using boxes stored in 3D hashmap with faces
- **CAD.sh** — `cad2d`: simple 2D CAD tool and analytical geometry demonstration
- **CastleBuilder.sh** — `CastleBuilder_main`: castle/building construction demo
- **FormationTactics.sh** — `FormationTactics_main`: sophisticated TotalWar-like simulation of historical battles
- **LandCraft.sh** — `LandCraft_main`: economy and logistics simulation game (inspired by Transport Tycoon, Cities: Skylines, Rise of Industry, Factorio). Goal: efficiently use terrain and natural resources to build civilization.
- **LandTactics.sh** — `LandTactics_main`: tactical simulation of modern battle (WWII) with units using terrain cover, mobility in different terrain, line-of-sight in 2D top view
- **MinimalTactics.sh** — `MinimalTactics_main`: quick test of minimal tactical simulation with simplified combat models
- **MolecularEditor.sh** — `MolecularEditor_main`: interactive molecular editor and MD simulation
- **MultiFight3D.sh** — `MultiFight3D_main`: 3D multiplayer combat simulation
- **NavalBattle.sh** — `NavalBattle_main`: simulation of battleship with multiple cannon turrets
- **NonIneritial.sh** — `NonInertial_main`: combat simulation (like Liero/Worms real-time) in non-inertial frame of reference
- **SailWar.sh** — `SailWar_main`: simulation of sailing ship — naval battle from age of sail
- **ShapePainter.sh** — `ShapePainter_main`: drawing program where brush makes useful shapes (quick hybrid between pixel and vector graphics)
- **SwordPlay.sh** — `SwordPlay_main`: soft-body sword play — martial combat system inspired by Mortal Kombat and Blade Symphony
- **Tanks.sh** — `Tanks_main`: tank simulation game with turret control, vehicle movement on cubic-spline terrain, and armor penetration calculation

## sketches_SDL/ — Engine Feature Demos

### test_2D.sh — 2D Sketches

Compile and run 2D SDL demos. Uncomment one `name=` line in the script. Demos:

- `test_AppSDL2OGL` — minimal boilerplate 2D SDL application
- `test_HashMap2D` — 2D spatial hashing visualization
- `test_HashMap2D_3` — performance test for adding points to HashMap2D
- `test_HashMap2D_uniformity` — evaluates uniformity of point distribution in HashMap2D
- `test_SphereTree2D` — 2D hierarchical bounding sphere tree
- `test_TileTree2D` — 2D spatial data structure for rectangular tiles
- `test_Clustering2D` — 2D clustering algorithm demonstration
- `test_NBodyColHashMap` — 2D N-body collision simulation using HashMap2D
- `test_NBodyWorld` — 2D N-body gravitational system simulation
- `test_PolyLine` — PolyLine class for connected line segments
- `test_Voronoi` / `test_Voronoi2` — 2D Voronoi diagram generation and visualization
- `test_BranchFract` — 2D fractal branching pattern
- `test_Mech2D` — 2D mechanical system simulation
- `test_MechEuler2D` — 2D mechanical simulation with Euler integration
- `test_AutoMesh2D` — automatic 2D mesh generation
- `test_SuperSonic2D` — 2D supersonic fluid flow
- `test_CommodityNetwork` — 2D commodity transport network simulation
- `test_SimplexGrid` — triangular grid data structure demo
- `test_Fluid2D` — 2D fluid dynamics simulation
- `test_TerrainHydraulics` — hydraulic erosion on 2D terrain
- `test_AnalyticalMushroomVortex` — analytical mushroom vortex visualization
- `test_TerrainCubic` — 2D terrain with cubic interpolation
- `test_TerrainRBF` — 2D terrain using Radial Basis Functions
- `test_CityGen` — procedural 2D city generation
- `test_Tris2Rect` — triangular to rectangular data conversion
- `test_GlobOpt2D` — 2D global optimization for molecular systems
- `test_Plotting2D` — Plot2D library for 2D graphs
- `test_Integration1D` — numerical integration methods for 1D functions
- `test_ConvexApprox1D` — convex approximation for 1D functions
- `test_PixelGlyphs` — pixel-based glyph rendering for text display

### test_3D.sh — 3D Sketches

Compile and run 3D SDL demos. Uncomment one `name=` line. Supports `args` variable. Demos:

- `test_QuatRotSampling` — orientation of surface directions on icosahedron (Quat4f.h, RotationMesh.h)
- `test_Projection` — debugging OpenCL projection frustum by point-clouds
- `test_Solids` — Platonic solids (Solids.h) and simple mesh (CMesh.h)
- `test_SceneGraph` — SceneGraph.h, especially Scene::Group
- `test_SphereSampling` — heightmap generation on planet/asteroid (SphereSampling.h, Noise.h). Args: `-testRand`, `-octmap`
- `test_Patches` — C1,C2-continuous triangle patches with heightmap
- `test_GUI` — GUI.h demo (text, Table, GUITextInput, Plot2D, DropDownList, ScisorBox)
- `test_EditorGizmo` — manipulate selected points in space (EditorGizmo.h)
- `test_MusculeEditor` — organic muscle-like shape editing (spline_hermite loft)
- `test_Camera` — 1st-person camera and crosshair for combat simulator
- `test_MousePicking` — mouse picking of 3D objects (raytrace.h raySphere())
- `test_SphereTree` — 3D diffusion-limited aggregation (DLA)
- `test_Raytracing` — ray-triangle intersection and occlusion
- `test_Radiosity` — radiosity of light labyrinth of rectangular corridors
- `test_Scatterer` — particle scattering on thin surfaces, flux transport through channels
- `test_RayScattererMMC` — volumetric scattering with Monte Carlo method
- `test_Elasticity` — legacy linearized elasticity (Truss.h, SoftBody.h) — prefer TrussDynamics_d.h
- `test_Electromagnetic` — thin hot plasma with ions in magnetic bottle (PotentialFlow.h, Fourier.h)
- `test_MultipoleAccel` — multipole acceleration (kBoxes.h HierarchicalKpivot)
- `test_Mesh` — mesh operations (fromPlanes(), findEdges(), collapseEdge())
- `test_RigidBody` — rigid body dynamics hanging on strings (Body.h)
- `test_VortexLattice` — vortex lattice velocity field for flow around wing (lift-line)
- `test_Collision` — RigidBody collision with spline HeightMap
- `test_CompressiveParticles` — high velocity compressible fluid impact on wedge
- `test_AABBTree` — axis-aligned bounding boxes (AABBTree3D.h, kBoxes.h)
- `test_BlockBuilder` — 3D-hashing box-building system for panel-housing game
- `test_TrussBuilder` — interactive truss builder [KP_1-7], export to SoftBody simulation [SPACE], save/load [k,l]
- `test_SphereGaussSeidel` — pack boxes close together respecting bounding boxes
- `test_Stick` — dynamics of straight sticks hitting a sphere (MMFF.h)

### test_OGL3.sh — OpenGL 3+ Demos

Compile and run OGL3 demos (requires GLEW). Uncomment one `name=` line. Supports `args` variable. Demos:

- `test_DiffractShader` — diffraction shader
- `test_SphereShader` — raycasted sphere impostor
- `test_ShaderDepth` — depth / shadow map testing
- `test_LandScape` — landscape rendering
- `test_DrawOGL3` — general DrawOGL3 tests
- `test_Instances` — instanced meshes
- `test_Atoms` — atom impostor spheres
- `test_Vegetation` — vegetation rendering
- `test_Horizont` — sky/horizon rendering
- `test_VAOs` — VAO usage
- `test_PatchesOGL3` — patch rendering
- `test_StencilTextures` — stencil textures
- `test_VolumetricTexture` — volumetric textures
- `test_RenderToTexture` — render-to-texture pipeline
- `test_Sprites` — billboard sprites
- `test_AntiAliasing` — anti-aliasing
- `test_MeshOGL3` — mesh rendering
- `test_GeometryShader` — geometry shader
- `test_Tubes` — tube rendering
- `test_OrbitalRayMarch` — orbital ray-marching
- `test_Texture` — basic textured rendering
- `test_SSAO` — Screen Space Ambient Occlusion
- `test_ScreenOGL3` — ScreenSDL2OGL3 base class
- `meshviewer` — generic mesh viewer (args: `common_resources/my_mesh.obj`)

### Molecular/run_test.sh — Molecular Force Field Demos

Compile and run molecular dynamics / force field demos. Uncomment one `name_def=` line or pass name as argument. Demos:

- `test_eFF` — Electron Force Field simulation
- `test_eFFMC` — compares eFF and CLCFGO force fields
- `test_eFF_old` — older version of eFF
- `test_CLCFGO` — Compact Linear Combination of Floating Gaussian Orbitals (multi-center eFF)
- `test_CLCFSF` — Compact Linear Combination of Finite Spherical Functions (multi-center eFF)
- `test_ConfDynamics` — sampling empty free space in collision potential
- `test_SoftMolecularDynamics` — soft molecular dynamics (known issue: large balls, large forces)
- `test_FARFF` — Flexible Atom sp-hybridization Reactive Force Field
- `test_RARFF` — Rigid Atom Reactive Force Field
- `test_RARFF2` — Rigid Atom Reactive Force Field 2
- `test_RARFFarr` — Rigid Atom Reactive Force Field (array-based)
- `test_RARFF_SR` — large scale grid-based Bucket accelerated RARFF
- `test_PBD_LJ_cluster` — Position-Based Dynamics with Lennard-Jones potential and cluster constraints
- `test_RspFF` — Flexible Atom sp-hybridization forcefield (FspFFclean.h)
- `test_sp3space` — relaxation in sp3 orbital space
- `test_spRotations` — rotations of sp3 orbitals (Slater-Koster tables, AutoApprox)
- `test_BondAdaptedMesh` — optimal sampling grid around molecule
- `test_MMFFmini` — minimal MMFF force field (*known issue: exits with pairMask error*)
- `test_Multipoles` — multipole expansion (*not working*)
- `test_EOFF` / `test_ESFF` / `test_FTRFF` / `test_RRFF` — various force field demos (*not rendering*)

## Orbital/ — Spacecraft & Orbital Mechanics Demos

- **SolarSystemMap.sh** — `SolarSystemMap`: solar system map viewer
- **spaceTactics.sh** — `spaceTactics`: space combat game with N-body physics, trajectory splines, weapons (railguns, lasers, Whipple shields), time scrubbing
- **spaceCraftEditor.sh** — `spaceCraftEditor`: interactive spacecraft editor. Loads Lua ship scripts (e.g. `ship_ICF_marksman_2.lua`). Args: `-dt`, `-method`, `-perframe`, `-s`
- **SpaceCraftEditorNew.sh** — `spaceCraftEditorNew`: newer version of spacecraft editor
- **spaceCraftDynamics.sh** — `spaceCraftDynamics`: spacecraft dynamics simulation (truss physics, weapons, damage)
- **spaceCraftDynamicsOCL.sh** — `spaceCraftDynamicsOCL`: GPU-accelerated (OpenCL) spacecraft dynamics
- **spaceCraftMeshExport.sh** — `spaceCraftMeshExport`: export spacecraft mesh to `.obj`/`.truss` files (headless, no GUI)
- **SpaceTacticsLib.sh** — `SpaceTacticsLib`: shared library test for space tactics combat models
- **optContinuousThrust.sh** — `test_OptContinuousThrust`: trajectory optimization with continuous thrust
- **constructionBlock.sh** — `constructionBlockApp`: block-based spacecraft construction demo
- **trussSimBatch.sh** — `trussSimBatch`: batch truss simulation (headless)
- **test_2D.sh** / **test_3D.sh** — copies of the 2D/3D sketch scripts (same demos as sketches_SDL/)
- Also contains output logs, `.obj` mesh exports, and `test_spacetactics.py` Python test

## LandCraft/ — LandCraft Library & Terrain Demos

- **LandCraft.sh** — `LandCraft_main`: full LandCraft economy/logistics simulation game
- **LandCraftLib.sh** — `LandCraftLib`: shared library test for LandCraft C API (terrain, hydraulics, roads, vehicles, economy, pathfinding)
- **test_landcraft.py** — Python test script using ctypes bindings to LandCraftLib, saves ground/water arrays as PGM images and plots hydraulic data with matplotlib
- Contains `ground.bin` / `water.bin` binary terrain data (2M each)

## Particle_In_Cell/ — PIC Method Notes

- **README.md** — documentation only (no runnable demo). Describes Particle-In-Cell algorithms: projection, grid-field solver, interpolation. Links to ShaderToy PIC implementations.
