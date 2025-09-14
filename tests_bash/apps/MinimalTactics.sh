#!/bin/bash

name=MinimalTactics_main  # quick test of minimal tactical simulation, testing simplified combat-models

dir=../../cpp/Build/apps/MinimalTactics   
ln -s ../../cpp/common_resources
rm -rf data
ln -s ../../cpp/apps/MinimalTactics/data

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