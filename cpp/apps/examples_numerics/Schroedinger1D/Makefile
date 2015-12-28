OBJS = program.cpp
OBJ_NAME = program.x
LFLAGS=-I/usr/local/include/SDL2 -lSDL2

CFLAGS= -g -Og
#CFLAGS= -w -O3
#CFLAGS= -w -Ofast

all : $(OBJS)
	g++ $(OBJS) $(LFLAGS) $(CFLAGS) -o $(OBJ_NAME)
