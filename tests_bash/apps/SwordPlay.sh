#!/bin/bash

# Experimenting with soft-body sword play - Martial Combat system inspired by Mortal Combat and Blade Symphony

name=SwordPlay_main                 
dir=../../cpp/Build/apps/SwordPlay
ln -s ../../cpp/common_resources
rm -rf data
ln -s ../../cpp/apps/SwordPlay/data

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