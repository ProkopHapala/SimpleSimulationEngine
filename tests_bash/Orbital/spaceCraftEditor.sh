#!/bin/bash

# ====== Paths
name=spaceCraftEditor
dir=../../cpp/Build/apps/OrbitalWar
ln -s ../../cpp/apps/OrbitalWar/data
ln -s ../../cpp/common_resources

# ====== Multiprocesing
ncpu=`nproc`
#ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
#echo "compile using ncpu="$ncpu
#OMP_NUM_THREADS=$ncpu
#export OMP_NUM_THREADS

# ====== Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
ln -s $dir/$name .

# ====== ASan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD

# ====== RUN

#./$name -s data/ship_ICF_interceptor_1.lua
./$name -s data/ship_ICF_marksman_1.lua
#./$name -s data/ship_NFPP_pendulum_1.lua
#./$name -s data/ship_NTR_marksman_1.lua
#./$name -s data/colony_1.lua




