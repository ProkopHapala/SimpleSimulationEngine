
#FLAGS = -std=c++11 -w
FLAGS = -std=c++11 -Og -g -w
#FLAGS = -std=c++11 -Og -g -Wall
#FLAGS = -std=c++11 -O3 -ftree-vectorize -unroll-loops
#FLAGS = -std=c++11 -O3 -ftree-vectorize -unroll-loops -ffast-math
#FLAGS = -std=c++11 -O3 -ftree-vectorize -unroll-loops --param max-unroll-times=4 -ffast-math -march=native -mtune=native -mssse3

INCLUDE_PATHS=-I../../../common/

all: clean program

program: 
	g++ $(FLAGS) -o program.x main.cpp $(INCLUDE_PATHS)  -lm -lSDL2 -lGL -lGLU

clean:
	rm -rf *.o program.x  

