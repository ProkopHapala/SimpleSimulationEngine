#!/bin/bash

name=SpaceTacticsLib
dir=../../cpp/Build/libs/OrbitalWar

# Link shared resources if missing
ln -s ../../cpp/common_resources 2>/dev/null

# ====== Multiprocessing
ncpu=`nproc`
ncpu=$(($ncpu - 1))
if [ $ncpu -lt 1 ]; then
    ncpu=1
fi
echo "compile using ncpu="$ncpu
export OMP_NUM_THREADS=$ncpu

# ====== Compilation
wd=`pwd`
cd $dir || exit 1
pwd
rm -f $name
make -j$ncpu $name || exit 1
cd $wd
rm -f $name.x
ln -s $dir/lib$name.so ./$name.x 2>/dev/null

# ====== ASan config (match LandCraft pattern)
export ASAN_OPTIONS=verify_asan_link_order=0:detect_leaks=0:halt_on_error=1:abort_on_error=1:print_stats=0
export LSAN_OPTIONS=detect_leaks=0:report_objects=0:verbosity=0

# ====== RUN
python3 -u test_spacetactics.py "$@" 2>&1 | tee test_spacetactics.log
