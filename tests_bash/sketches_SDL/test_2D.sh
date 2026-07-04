#!/bin/bash

# This script compiles and runs a specified 2D demo from the SimpleSimulationEngine.
# It mirrors the structure of test_3D.sh.

# ====== Paths

# Uncomment one of the following lines to select the demo to run:
#name=test_AppSDL2OGL              # minimal boilerplate 2D SDL application
#name=test_HashMap2D               # 2D spatial hashing visualization
#name=test_SphereTree2D            # 2D hierarchical bounding sphere tree
#name=test_TileTree2D              # 2D spatial data structure for rectangular tiles
#name=test_Clustering2D            # 2D clustering algorithm demonstration
#name=test_HashMap2D_uniformity    # Evaluates uniformity of point distribution in HashMap2D
#name=test_HashMap2D_3             # Performance test for adding points to HashMap2D
#name=test_NBodyColHashMap          # 2D N-body collision simulation using HashMap2D
name=test_NBodyWorld              # 2D N-body system simulation
#name=test_PolyLine                # Showcases the PolyLine class for connected line segments
#name=test_Voronoi                 # Generates and visualizes 2D Voronoi diagrams
#name=test_Voronoi2                # Enhanced visualization of 2D Voronoi diagrams
#name=test_BranchFract             # Generates and visualizes a 2D fractal branching pattern
#name=test_Mech2D                  # Simulates a 2D mechanical system
#name=test_MechEuler2D             # 2D mechanical simulation highlighting Euler integration
#name=test_AutoMesh2D              # Demonstrates automatic 2D mesh generation
#name=test_SuperSonic2D            # Simulates 2D supersonic fluid flow
#name=test_CommodityNetwork        # Simulates a 2D commodity network
#name=test_SimplexGrid             # Demonstrates a 2D data structure based on a triangular grid
#name=test_Fluid2D                 # Simulates 2D fluid dynamics
#name=test_TerrainHydraulics       # Simulates hydraulic erosion on a 2D terrain
#name=test_AnalyticalMushroomVortex # Visualizes an analytical solution for a mushroom vortex
#name=test_TerrainCubic            # Demonstrates 2D terrain rendering using cubic interpolation
#name=test_TerrainRBF              # Generates 2D terrain using Radial Basis Functions
#name=test_CityGen                 # Demonstrates procedural 2D city generation
#name=test_Tris2Rect               # Explores conversion between triangular and rectangular data
#name=test_GlobOpt2D               # Demonstrates 2D global optimization for molecular systems
#name=test_Plotting2D              # Showcases the Plot2D library for 2D graphs
#name=test_Integration1D           # Demonstrates numerical integration methods for 1D functions
#name=test_ConvexApprox1D          # Illustrates convex approximation for 1D functions
#name=test_PixelGlyphs             # Demonstrates pixel-based glyph rendering for displaying text

dir=../../cpp/Build/sketches_SDL/2D
ln -s ../../cpp/common_resources

# ====== Multiprocessing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# ====== Compilation
wd=`pwd`
cd $dir
pwd
rm -f $name # Use -f to avoid error if file doesn't exist
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
rm $name.x
ln -s $dir/$name ./$name.x

# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD
#export ASAN_OPTIONS=detect_leaks=1

# ====== RUN

./$name.x