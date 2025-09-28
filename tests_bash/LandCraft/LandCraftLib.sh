#!/bin/bash

# Economy and Logistics simulation game inspired by Transport Tycoon, City Skylenes, Rise of Indrustry, Factorio etc. 
# The goal is to efficiently use terrain and natural resorces to build civilization

name=LandCraftLib
dir=../../cpp/Build/libs/CombatModels   # 
ln -s ../../cpp/common_resources
rm -rf data
ln -s ../../cpp/apps/LandCraft/data

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
# LD_PRELOAD=$(g++ -print-file-name=libasan.so)
# echo   $LD_PRELOAD
# export LD_PRELOAD
# Silence link-order warning and leaks spam; keep hard errors
export ASAN_OPTIONS=verify_asan_link_order=0:detect_leaks=0:halt_on_error=1:abort_on_error=1:print_stats=0
export LSAN_OPTIONS=detect_leaks=0:report_objects=0:verbosity=0
#export ASAN_OPTIONS=detect_leaks=1

# ====== RUN

python3 -u test_landcraft.py 2>&1 | tee test_landcraft.log