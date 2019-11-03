#!/bin/bash

# ================ SETUP

#target="fast_atan2"
#target="invSphereMap"

target=$1

CPPC="g++"

build_path="Build"


ERRFLAGS="-Wno-missing-braces -Werror=return-type"

#CFLAGS="-std=c++17 -Og"
#CFLAGS="-std=c++17 -Og -Wall"
#CFLAGS="-std=c++17 -Og -Wall -Wno-missing-braces -Werror=return-type" 
#CFLAGS="-std=c++17 -O2    -march=native -mtune=native"
#CFLAGS="-std=c++17 -O3    -march=native -mtune=native"
CFLAGS="-std=c++17 -Ofast -march=native -mtune=native"

IFLAGS="-I../../common/math -I../../common_SDL/SDL2OGL -I../../common_SDL/ -I../../common/utils -I../../common/dataStructures -I../../common_SDL/SDL2OGL -I../../libs_SDL/GLView -I/usr/include" 
LFLAGS="-L../../Build/libs_SDL/GLView -lGLView -lGL -lSDL2"

# ================ Main

mkdir $build_path

"
dirbak=`pwd`
src_path=../../libs_SDL/GLView
bin_path=../../Build/libs_SDL/GLView
cd $bin_path
make GLView
#echo "!!! ==== GLView compiled ... go back .. "
cd $dirbak

ln -f -s $bin_path/libGLView.so  $build_path/libGLView.so
ln -f -s $bin_path/libGLView.so  libGLView.so
"


rm $build_path/$target.x
#g++ -o testGLV.x test.cpp $IFLAGS $LFLAGS 
$CPPC -o $build_path/$target.x $target.cpp $CFLAGS $ERRFLAGS $IFLAGS $LFLAGS 

#gcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2
#tcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2

./$build_path/$target.x

