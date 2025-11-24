#!/bin/bash

# ====== Paths
name=spaceCraftMeshExport
dir=../../cpp/Build/apps/OrbitalWar

ln -s ../../cpp/apps/OrbitalWar/data          2>/dev/null
ln -s ../../cpp/common_resources              2>/dev/null

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
rm -f $name
make -j$ncpu $name
cd $wd
ln -sf $dir/$name ./$name.x

# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

# ====== RUN
# Simple default: OBJ only
#./$name.x -s data/ship_ICF_marksman_2.lua -o ship_ICF_marksman_2.obj -v 1 | tee OUT-spaceCraftMeshExport

# OBJ + TRUSS
./$name.x -s data/ship_ICF_marksman_2.lua \
    -o ship_ICF_marksman_2.obj \
    -t ship_ICF_marksman_2.truss \
    -v 1 | tee OUT-spaceCraftMeshExport
