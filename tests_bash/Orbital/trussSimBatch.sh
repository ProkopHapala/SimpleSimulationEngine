#!/bin/bash

name=trussSimBatch
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



# default input/output (user can edit as needed)
IN_TRUSS=${1:-ship_ICF_marksman_2.truss}
OUT_TRAJ=${2:-ship_ICF_marksman_2_traj.xyz}

./$name.x -i "$IN_TRUSS" -o "$OUT_TRAJ" -n 1000 -dt 5e-4 -psave 10 -v 1 | tee OUT-trussSimBatch
