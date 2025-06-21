#!/bin/bash

# ====== Paths

#name=test_AABBTree              # axis-aligned bounding boxes using AABBTree3D.h and kBoxes.h
#name=test_BlockBuilder          # using 3D-hashing box-building system as base for panel-housing game
#name=test_Camera                # 1st-person camera and crosshair as a base for shooting or action air/space combat simulator
#name=test_Collision             # collision of RigidBody with spline HeightMap using Body.h and Terrain2D.h
#name=test_CompressiveParticles  # high velocity impact of compressible fluid on wedge boundary using CompressiveParticles.h
#name=test_EditorGizmo           # using EditorGizmo.h to manipulate selected points in space 
#name=test_Elasticity            # legacy linearized elasticity solver using Truss.h and  SoftBody.h (NOTE: prefer TrussDynamics_d.h instead)
#name=test_Electromagnetic       # [m] to run simulation of thin hot plasma with ions moving in magnetic bottle (feeling agregate electric and magnetic fields) using PotentialFlow.h and poisson solver using Fourier.h 
#name=test_GUI                   # using GUI.h ( 2D/3D text, Table, GUITextInput  Plot2D, DropDownList, ScisorBox ) 
#name=test_Mesh                  # various operations on mesh using Mesh.h, e.g. convex polygon by OMesh::fromPlanes(), findEdges(), colapseEdge()
#name=test_MousePicking          # mouse picking of 3D object using raytrace.h raySphere(), NOTE: very similar to test_Camera, should merge (?)
#name=test_MultipoleAccel        # Multipole acceleration using kBoxes.h HierarchicalKpivot 
#name=test_MusculeEditor         # edit Organic Muscule-like shapes using loft of spline_hermite.h and EditorGizmo.h 
#name=test_Patches               # rendering C1,C2-continuous triangle patches with heightmap using spline_triC1.h and TerrainSimplex.h 
#name=test_Projection            # debugging OpenCL projection frustum by point-clouds
#name=test_QuatRotSampling       # visualize orietation of surface diractions (u,v) on surface of icosahedron using Quat4f.h and RotationMesh.h, usefull for continuous maps on surface of sphere
name=test_Radiosity             # use Radiosity.h to calculate radiosity of light labirynth of rectangular corridors sampled by surface points on uniform grid
#name=test_RayScattererMMC       # use RayScatter.h to calculate volumetric scattering of light or particle with Monte Carlo method
#name=test_Raytracing            # use raytrace.h rayTriangles() to calculate intersection of ray with triangles and occlusion of other polygon
#name=test_RigidBody             # rigid body dynamics hanging on strings using Body.h RigidBody class
#name=test_Solids                # demo and testing of Solids.h (with Platonic solids) and CMesh.h (with simple mesh class) 
#name=test_Scatterer             # scattering of particles on thin surfaces using Scatterer2.h to calculate flux transport through a network of channels connecting scattering elements
#name=test_SceneGraph            # showcase SceneGraph.h especially  Scene::Group
#name=test_SphereGaussSeidel     # showcase SphereGaussSeidel.h to quickly pack boxes close to each other (respecting their bounding boxes)
#name=test_SphereSampling        # demonstrate generation of heightmap on planet or asteroided using SphereSampling.h Noise.h and DrawSphereMap.h
#name=test_SphereTree            # 3D Diffusion-limited aggregation (DLA) using CubicRuler.h and std::unordered_multimap<int_fast64_t,int>  grid;
#name=test_Stick                 # dynamics of straight sticks hitting a sphere using MMFF.h
#name=test_TrussBuilder          # Create Truss by keyboard input [KP_1-7] using TrussBuilder.h and export it to SoftBody.h by and run simulation [SPACE], [k,l] to save/load from file
#name=test_VortexLattice         # PotentialFlow.h to calculate velocity field of vortex lattice for flow around wing (represented by lift-line), contains function for numerical integration of velocity field (to verify analytical formulas)




dir=../../cpp/Build/sketches_SDL/3D
ln -s ../../cpp/apps/OrbitalWar/data
ln -s ../../cpp/common_resources

# ====== Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# ====== Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
rm $name.x
ln -s $dir/$name ./$name.x

# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

# ====== RUN

./$name.x



