#!/bin/bash

target="fast_atan2"

CPPC="g++"

#CFLAGS="-std=c++17 -Og"
#CFLAGS="-std=c++17 -Og -Wall"
CFLAGS="-std=c++17 -Og -Wall -Wno-missing-braces -Werror=return-type" 
#CFLAGS="-std=c++17 -O2    -march=native -mtune=native"
#CFLAGS="-std=c++17 -O3    -march=native -mtune=native"
#CFLAGS="-std=c++17 -Ofast -march=native -mtune=native"

IFLAGS="-I../../common/math -I../../common_SDL/SDL2OGL -I../../common_SDL/ -I../../common/utils -I../../common/dataStructures -I../../common_SDL/SDL2OGL -I../../libs_SDL/GLView -I/usr/include" 
LFLAGS="-L../../Build/libs_SDL/GLView -lGLView -lGL -lSDL2"

# ================ SETUP

dirbak=`pwd`
src_path=../../libs_SDL/GLView
bin_path=../../Build/libs_SDL/GLView
cd $bin_path
make GLView
#echo "!!! ==== GLView compiled ... go back .. "
cd $dirbak

ln -f -s $bin_path/libGLView.so  ./libGLView.so

rm $target.x
#g++ -o testGLV.x test.cpp $IFLAGS $LFLAGS 
$CPPC -o $target.x $target.cpp $IFLAGS $LFLAGS 

#gcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2
#tcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2

./target.x

