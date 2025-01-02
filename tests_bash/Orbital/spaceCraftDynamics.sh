#!/bin/bash

# ====== Paths
name=spaceCraftDynamics
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
ln -s $dir/$name ./$name.x

# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

# ====== RUN

#./$name.x -float -fix 2
./$name.x -double -fix 2 -shape 2,200 -method 4 -perframe 10 -nsolve 1

#./$name.x -s data/ship_ICF_interceptor_1.lua
#./$name.x -s data/ship_ICF_marksman_1.lua
#./$name.x -s data/ship_ICF_marksman_2.lua
#./$name.x -s data/ship_NFPP_pendulum_1.lua
#./$name.x -s data/ship_NTR_marksman_1.lua
#./$name.x -s data/colony_1.lua




