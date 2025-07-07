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

#./$name.x -skelet 1
#./$name.x -skelet 0
#./$name.x -parabola
./$name.x -panel

#./$name.x -blocks
#./$name.x -extrude_octahedron
#./$name.x -oct_nodes
#./$name.x -bevel | tee OUT


