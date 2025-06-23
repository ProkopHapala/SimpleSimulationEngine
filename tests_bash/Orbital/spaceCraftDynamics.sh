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
#./$name.x -method 2 -double -fix 2 -shape 2,200  -perframe 10 

#./$name.x -method 2 -dt 0.1    -double -shape 2,20  -perframe 10
#./$name.x -method 2 -dt 0.01   -double -fix 0 -shape 2,10 -perframe 10
#./$name.x -method 2 -dt 0.0025 -double -fix 2 -shape 2,10 -perframe 400
#./$name.x -method 2 -dt 0.01   -double -fix 2 -shape 2,10 -perframe 100  -G 0.0,-9.81,0.0

#./$name.x -method 4  -dt 0.01 -double -fix 2 -shape 2,10  -perframe 10  -G 0.0,-9.81,0.0  -nsolve 30
#./$name.x -method 11 -dt 0.01 -double -fix 2 -shape 2,10  -perframe 10  -G 0.0,-9.81,0.0  -nsolve 30
#./$name.x -method 12 -dt 0.01 -double -fix 2 -shape 2,10  -perframe 10  -G 0.0,-9.81,0.0  -nsolve 30

#./$name.x -method 2  -dt 0.01 -double -fix 2 -shape 2,50  -perframe 10  -G 0.0,-9.81,0.0 
#./$name.x -method 11 -dt 0.01 -double -fix 2 -shape 2,50  -perframe 10  -G 0.0,-9.81,0.0  -nsolve 10
#./$name.x -method 12 -dt 0.01 -double -fix 2 -shape 2,50  -perframe 10  -G 0.0,-9.81,0.0  -nsolve 10 -bmix 0.9


#./$name.x -method 2 -dt 0.01 -double -fix 2 -shape 2,20 -perframe 10  -G 0.0,-9.81,0.0 
#./$name.x -method 3 -dt 0.01 -double -fix 2 -shape 2,20 -perframe 10  -G 0.0,-9.81,0.0 


#./$name.x -method 11 -dt 0.01 -double -fix 2 -shape 2,20 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 5
#./$name.x -method 12 -dt 0.01 -double -fix 2 -shape 2,20 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 5 -bmix 0.7
#./$name.x -method 12 -dt 0.01 -double -fix 2 -shape 2,20 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 5 -bmix 0.9
#./$name.x -method 12 -dt 0.01 -double -fix 2 -shape 2,20 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 5 -bmix 0.95





#./$name.x -method 13 -dt 0.01 -double -fix 2 -shape 2,50 -perframe 10  -G 0.0,-9.81,0.0  -nsolve 10
#./$name.x -method 2  -dt 0.01 -double -fix 2 -shape 2,50 -perframe 10  -G 0.0,-9.81,0.0  


#./$name.x -method 2  -dt 0.05 -double -fix 2 -shape 2,2 -perframe 100 -G 0.0,-9.81,0.0
#./$name.x -method 3  -dt 0.05 -double -fix 2 -shape 2,2 -perframe 100 -G 0.0,-9.81,0.0


#./$name.x -method 2  -dt 0.01  -double -fix 0  -shape 4,2  -perframe 10
#./$name.x -method 2  -dt 0.001 -double         -shape 4,2  -perframe 1
#./$name.x -method 2  -dt 0.01  -double -fix 0  -shape 4,20 -perframe 100


#./$name.x -method 6 -double -fix 2 -shape 2,200 -perframe 1  -nsolve 30 -bmix 3,0.95
#./$name.x -method 8 -double -fix 2 -shape 2,200 -perframe 10 -nsolve 20 -bmix 3,0.75
#./$name.x -method 8 -double -fix 2 -shape 2,200 -perframe 10 -nsolve 20 -bmix 3,0.75


#./$name.x -s data/ship_ICF_interceptor_1.lua
#./$name.x -s data/ship_ICF_marksman_1.lua
#./$name.x -fix 2 -s data/ship_ICF_marksman_2.lua

#./$name.x -debug_orig   -method 12  -dt 0.01 -omega 0.0,0.0,0.05 -s data/ship_ICF_marksman_2.lua  -perframe 10  -nsolve 10 -bmix 3,0.8  | tee OUT_orig
#./$name.x               -method 12  -dt 0.01 -omega 0.0,0.0,0.05 -s data/ship_ICF_marksman_2.lua  -perframe 10  -nsolve 10 -bmix 3,0.8  | tee OUT_orig


./$name.x   -omega 0.0,0.0,50.5 -dt 0.01   -oct_nodes



#./$name.x -method 2   -dt 0.01 -omega 0.0,0.0,0.05 -s data/ship_ICF_marksman_2.lua  -perframe 10


#./$name.x  -method  6  -fix 2 -s data/ship_ICF_marksman_2.lua    -nsolve 50 -bmix 3,0.95
#./$name.x -s data/ship_NFPP_pendulum_1.lua
#./$name.x -s data/ship_NTR_marksman_1.lua
#./$name.x -s data/colony_1.lua




