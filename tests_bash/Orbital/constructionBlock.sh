#!/bin/bash

# ====== Paths

name=constructionBlockApp
dir=../../cpp/Build/apps/OrbitalWar
ln -s ../../cpp/apps/OrbitalWar/data
ln -s ../../cpp/common_resources

# ====== Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# ====== Compilation
start_time=$(date +%s.%N)
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
rm $name.x
ln -s $dir/$name ./$name.x
end_time=$(date +%s.%N)
elapsed_time=$(echo "scale=2; $end_time - $start_time" | bc)
echo "Compilation time: $elapsed_time [s]"
#exit 0


# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

# ====== RUN

#./$name.x -skelet 1
#./$name.x -skelet 0
#./$name.x -parabola
#./$name.x -blocks
#./$name.x -extrude_octahedron
#./$name.x -oct_nodes | tee OUT-constructionBlockApp
#./$name.x -cube_nodes | tee OUT-constructionBlockApp
#./$name.x -bevel | tee OUT
#./$name.x -panel | tee OUT


# TorusSheet(truss, {4,6}, {-0.125,0.0}, {0.875,1.0}, {5.0,20.0}, 0b0011, 0.0 );
#./$name.x -TorusSheet 4,6 -0.125,0.0 0.875,1.0 5.0,20.0 11 | tee OUT-constructionBlockApp
# TubeSheet(truss, {4,10}, {0.0,0.0}, {1.0,1.0}, {10.0,10.0}, 10.0, 0b1011, 0.5 );
#./$name.x -TubeSheet 4,10 0.0,0.0 1.0,1.0 10.0,10.0,10.0 1101 | tee OUT-constructionBlockApp

#QuadSheet(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, 0b1011, Quat4i{0,0,0,0}, 3,15 );
#./$name.x -QuadSheet 10,10 0.0,0.0,0.0 86.602540378,50.0,0.0 0.0,100.0,0.0 86.602540378,150.0,0.0 1101 3,15 | tee OUT-constructionBlockApp

# SlabTube(truss, {2,16}, {0.0,0.0}, {0.2,1.5*M_PI}, {10.0,10.0}, 10.0, {0.333333,0.333333,2.0}, 0b101010111, Quat4i{0,0,0,0} );
#./$name.x -SlabTube 2,16 0.0,0.0 0.2,6.0 10.0,10.0,10.0 101010111 | tee OUT-constructionBlockApp
#./$name.x -SlabTube 2,16 0.0,0.0 0.2,6.0 10.0,10.0,10.0  0.333333,0.333333,2.0 111010101 | tee OUT-constructionBlockApp

#QuadSlab(truss, {10,10}, {0.0,0.0}, {1.0,1.0}, {0.0,0.0,0.0}, {86.602540378,50.0,0.0}, {0.0,100.0,0.0}, {86.602540378,150.0,0.0}, {0.333333,0.333333,7.0}, 0b101010111, Quat4i{0,0,0,0} );
#./$name.x -QuadSlab 10,10 0.0,0.0,0.0 86.602540378,50.0,0.0 0.0,100.0,0.0 86.602540378,150.0,0.0 0.333333,0.333333,7.0 101010111 | tee OUT-constructionBlockApp
./$name.x -QuadSlab 10,10 0.0,0.0,0.0 86.602540378,50.0,0.0 0.0,100.0,0.0 86.602540378,150.0,0.0 0.333333,0.333333,7.0  111010101 | tee OUT-constructionBlockApp

# ====== Parabolic Mesh Generators (ported from JS MeshesUV.js)
# Run headless: generate mesh, export .obj + .svg, exit without GUI loop

# ParabolaSheet: single-layer parabolic truss (magnetic nozzle)
# ParabolaSheet(n={nx,ny}, UVmin, UVmax, R, L, dirMask, twist)
./$name.x -headless -ParabolaSheet 6,10 0.2,0.0 1.0,1.0 10.0 10.0 1011 | tee OUT-parabola_sheet

# ParabolaSlab: double-layer parabolic with bracing (nozzle + damper structure)
# ParabolaSlab(n={nx,ny}, UVmin, UVmax, R, L, up={ox,oy,oz}, dirMask)
./$name.x -headless -ParabolaSlab 6,8 0.0,0.0 1.0,1.0 10.0 10.0 0.0,0.0,1.2 10001001111 | tee OUT-parabola_slab

# ParametricParabola: annular patch with independent inner/outer tessellation
# ParametricParabola(nTop, nBottom, nRows, R1, R2, L)
./$name.x -headless -ParametricParabola 13 4 4 0.5 3.0 3.0 | tee OUT-parabola_patch




