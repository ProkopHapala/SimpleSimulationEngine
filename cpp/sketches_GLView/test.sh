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

PREDIR="../"
PREI="-I"$PREDIR

echo "PREI " $PREI

#IFLAGS="-I../../common/math -I../../common_SDL/SDL2OGL -I../../common_SDL/ -I../../common/utils -I../../common/dataStructures -I../../common_SDL/SDL2OGL -I../../libs_SDL/GLView -I/usr/include" 
IFLAGS=$PREI"common/math "$PREI"common/CombatModels "$PREI"common_SDL/SDL2OGL "$PREI"common_SDL/ "$PREI"common/utils "$PREI"common/dataStructures "$PREI"common_SDL/SDL2OGL "$PREI"libs_SDL/GLView -I/usr/include" 

echo "IFLAGS " $IFLAGS

LFLAGS="-L"$PREDIR"Build/libs_SDL/GLView -lGLView -lGL -lSDL2"

src_path=$PREDIR/libs_SDL/GLView
bin_path=$PREDIR/Build/libs_SDL/GLView

# ================ Main

mkdir $build_path
dirbak=`pwd`
cd $bin_path
make GLView
#echo "!!! ==== GLView compiled ... go back .. "
cd $dirbak

#ln -f -s $bin_path/libGLView.so  libGLView.so
#ln -f -s libGLView.so  $build_path/libGLView.so
#cp libGLView.so $build_path
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$bin_path
echo $LD_LIBRARY_PATH



rm $build_path/$target.x
#g++ -o testGLV.x test.cpp $IFLAGS $LFLAGS 
$CPPC -o $build_path/$target.x $target.cpp $CFLAGS $ERRFLAGS $IFLAGS $LFLAGS 


#gcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2
#tcc -o testGLV.x test.c $LFLAGS -lGLV -lGL -lSDL2

./$build_path/$target.x

