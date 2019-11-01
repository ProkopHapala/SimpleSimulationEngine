#!/bin/bash


dirbak=`pwd`
src_path=../../libs_SDL/GLView
bin_path=../../Build/libs_SDL/GLView

cd $bin_path
make GLView
#mv GLV_lib.so libGLV.so
echo "!!! ==== GLView compiled ... go back .. "
cd $dirbak

IFLAGS="-I../../common/math -I../../common_SDL/SDL2OGL -I../../common/utils -I../../common/dataStructures -I../../common_SDL/SDL2OGL -I../../libs_SDL/GLView -I/usr/include" 
LFLAGS="-L../../Build/libs_SDL/GLView -lGLView -lGL -lSDL2"
#LFLAGS="-lGLV -lGL -lSDL2"

echo "!!!!!!! "

ln -f -s $bin_path/libGLView.so  ./libGLView.so

echo "!!!!!!! "

g++ -o testGLV.x test.cpp $IFLAGS $LFLAGS 
#gcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2
#tcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2

./testGLV.x

