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

# --- Examples (uncomment one below to run) ---
# NOTE: All examples call tests_bash/LandCraft/test_landcraft.py
# Terrain size defaults to 512x512 unless overridden via --nx --ny

# Example A: Default (load cache or generate), contour basin, apply, plot overlays
#python3 -u test_landcraft.py --basin contour --basin-apply 1 --plot 2

# Example B: Dijkstra (priority-flood) with apply and interactive plot
#python3 -u test_landcraft.py --basin dijkstra --basin-apply 1 --plot 2

# Example C: Bellman-Ford boundary with max 200 iterations, apply
#python3 -u test_landcraft.py --basin bellman --basin-iters 200 --basin-apply 1 --plot 2

# Example D: Contour with explicit level cap (no PQ); cap=300.0
#python3 -u test_landcraft.py --basin contour --basin-level-cap 300.0 --basin-apply 1 --plot 2

# Example E: Add explicit drain seed at (256,256) with initial water level 200
#python3 -u test_landcraft.py --basin dijkstra --drain-x 256 --drain-y 256 --drain-level 200 --basin-apply 1 --plot 2

# Example F: Different neighborhood (4-neigh) and smaller grid
#python3 -u test_landcraft.py --nx 256 --ny 256 --neighbors 8 --basin dijkstra --basin-apply 1 --plot 2

# Example G: Generate terrain via C++ pipeline (noise+erosion), then dijkstra
#python3 -u test_landcraft.py --terrain cpp --basin dijkstra --basin-apply 1 --plot 2

# Example H: Bisect noise with droplet erosion, then dijkstra
#python3 -u test_landcraft.py --terrain cpp_bisec --erosion-iters 200 --erosion-drops 20000 --erosion-steps 50 --basin dijkstra --basin-apply 1 --plot 2

# Example I: Relax/outflow baseline (no basin), set outflow and iterate
#python3 -u test_landcraft.py --basin none --outflow-x 0 --outflow-y 0 --relax-iters 200 --plot 2

# Example J: Use cache load/save folders
#python3 -u test_landcraft.py --load-prefix ./ --save-prefix ./ --basin dijkstra --basin-apply 1 --plot 2

python3 -u test_landcraft.py 2>&1 | tee test_landcraft.log