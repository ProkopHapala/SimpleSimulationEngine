#!/bin/bash

name=OpenCL_OpenGL
rm $name.x
g++ -o $name.x $name.cpp -I/usr/include/CL -I/usr/include/GL -I/usr/include/SDL2 -lGL -lGLEW -lSDL2 -lOpenCL
./$name.x