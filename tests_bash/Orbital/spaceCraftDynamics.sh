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
# LinSolveMethod{ CG=0, CGsparse=1, Cholesky=2, CholeskySparse=3, Jacobi=4, GaussSeidel=5 };
#./$name.x -float -fix 2
#./$name.x -double -fix 2 -shape 2,200 -method 2 -perframe 10 

#./$name.x -dt 0.1 -double -shape 2,20 -method 2 -perframe 10
#./$name.x -dt 0.01 -double -fix 0 -shape 2,10 -method 2 -perframe 10
#./$name.x -dt 0.0025 -double -fix 2 -shape 2,10 -method 2 -perframe 400
#./$name.x -dt 0.01 -double -fix 2 -shape 2,10 -method 2 -perframe 100  -G 0.0,-9.81,0.0

#./$name.x -dt 0.01 -double -fix 2 -shape 2,10 -method 4  -perframe 10  -G 0.0,-9.81,0.0  -nsolve 30
#./$name.x -dt 0.01 -double -fix 2 -shape 2,10 -method 11 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 30
#./$name.x -dt 0.01 -double -fix 2 -shape 2,10 -method 12 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 30
#./$name.x -dt 0.01 -double -fix 2 -shape 2,50 -method 12 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 10
#./$name.x -dt 0.01 -double -fix 2 -shape 2,50 -method 2  -perframe 10  -G 0.0,-9.81,0.0  

#./$name.x -dt 0.05 -double -fix 2 -shape 2,10 -method 2 -perframe 100   -G 0.0,-9.81,0.0


#./$name.x -dt 0.01 -double -fix 0 -shape 4,2 -method 2 -perframe 10
#./$name.x -dt 0.001 -double -shape 4,2 -method 2 -perframe 1
#./$name.x -dt 0.01 -double -fix 0 -shape 4,20 -method 2 -perframe 100


#./$name.x -double -fix 2 -shape 2,200 -method 6 -perframe 1 -nsolve 30 -bmix 3,0.95
#./$name.x -double -fix 2 -shape 2,200 -method 8 -perframe 10 -nsolve 20 -bmix 3,0.75
#./$name.x -double -fix 2 -shape 2,200 -method 8 -perframe 10 -nsolve 20 -bmix 3,0.75


#./$name.x -s data/ship_ICF_interceptor_1.lua
#./$name.x -s data/ship_ICF_marksman_1.lua
#./$name.x -fix 2 -s data/ship_ICF_marksman_2.lua

#./$name.x -dt 0.01 -omega 0.0,0.0,0.05 -s data/ship_ICF_marksman_2.lua  -method 12 -perframe 10  -nsolve 10
./$name.x -dt 0.01 -omega 0.0,0.0,0.05 -s data/ship_ICF_marksman_2.lua  -method 2 -perframe 10  


#./$name.x -fix 2 -s data/ship_ICF_marksman_2.lua   -method  6   -nsolve 50 -bmix 3,0.95
#./$name.x -s data/ship_NFPP_pendulum_1.lua
#./$name.x -s data/ship_NTR_marksman_1.lua
#./$name.x -s data/colony_1.lua




