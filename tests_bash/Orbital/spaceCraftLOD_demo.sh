#!/bin/bash
# spaceCraftLOD_demo.sh — end-to-end LOD review: sketch_prism OBJ → low-res + high-res .obj/.svg
# Outputs: review_lod_demo/prism_demo_{low,high}_res.{obj,svg} and low_res_tags.obj

name=spaceCraftMeshExport
dir=../../cpp/Build/apps/OrbitalWar
outdir=review_lod_demo
mkdir -p "$outdir"

ln -sf ../../cpp/apps/OrbitalWar/data data 2>/dev/null
ln -sf ../../cpp/common_resources common_resources 2>/dev/null

ncpu=`nproc`; ncpu=$(($ncpu - 1))
echo "compile using ncpu=$ncpu"
OMP_NUM_THREADS=$ncpu; export OMP_NUM_THREADS

wd=`pwd`
cd $dir && make -j$ncpu $name && cd $wd
ln -sf $dir/$name ./$name.x

LD_PRELOAD=$(g++ -print-file-name=libasan.so)
export LD_PRELOAD

ship=data/ship_lod_demo.lua
prefix=$outdir/prism_demo

echo "=== LOW-RES sketch LOD ==="
./$name.x -s $ship -lod sketch \
    -o ${prefix}_low_res.obj \
    -svg ${prefix}_low_res.svg \
    -g ${prefix}_low_res_tags.obj \
    -v 1 | tee OUT-spaceCraftLOD_demo-sketch

echo "=== HIGH-RES blocks LOD ==="
./$name.x -s $ship -lod blocks \
    -o ${prefix}_high_res.obj \
    -svg ${prefix}_high_res.svg \
    -v 1 | tee OUT-spaceCraftLOD-demo-blocks

echo ""
echo "Review artifacts in $outdir/:"
ls -la $outdir/
