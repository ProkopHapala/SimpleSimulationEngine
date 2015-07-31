## Simple Simulation Engine

a minimalistic engine for: 
- physical simulations 
- numerical math 
- game development 
- computer graphics
- educational purposes

written in C++11 with possible python and lua binding using SDL2 and OpenGL for Graphical user interface.

## Philosophy

Main motivation is to assemble various common algorithma from physics, numerical mathematics, computer graphics and computer science and build minimalistic tool for rapid development of physical simulation programs and games. 

**I should be:**
- It should be **didactic** - user is expected to be able to understand how things works inside ( i.e. it should not be a blackbox )
- It is focused on development of **small programs** ( like sketches ) 
- Each particular part of the engine should be illustrated by simplistic example program
- Minimalized dependence on 3rd party libraries 
- Minimalized dependence between different modules within the engine
- No instalation required - just include `.h` files directly into your program by `#include`
- It is rather a set of tools and pre-programed code samples rather than enclosed package
- It is expected that for particular project could be used just a few `.h` files independently ( i.e. not the engine as a whole )
- It does not follow rigorous "[encapsulation](https://en.wikipedia.org/wiki/Encapsulation_(computer_programming))" scheme of Object Oriented Programing and other rigorous software enginering method typical for development of large projects. e.g. all properties of objects are `public`. 
- User is expected and ecouraged to modify the code (including core parts) to match requirements of his particular project. It should be rather starting point of development rather than final product.

## Dependecies 
- The code aims to have as little dependecies as possible.
- The computational / simulation core should have no dependencies at all.
- However, part of the code is dedicated to visualization and user interface using [simple direct media layer 2](https://www.libsdl.org/) ( SDL2 ) and OpenGL
- other dependencies may be added in examples of particular use cases ( i.e. in praticular games and simulation programs )

<font color="green"> This is work in progress. Particular modules of the engine are functional, but they are not put together to build an unified system and work together. There are currectly no finished functional examples </font>
