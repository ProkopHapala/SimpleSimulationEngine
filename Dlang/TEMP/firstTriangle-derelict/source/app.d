import std.stdio;

//import derelict.sdl2.sdl;
//import derelict.opengl3.gl3;
//import derelict.util.exception;

import glUtils;

import shaders;

version(linux){pragma(lib, "dl");}

WindowSDL window1;
//GLuint programID;

Shader shader1;
GLuint vertexArrayID;
GLuint vbID_1;

bool frame(SDL_Window* window){

    // event loop
    SDL_Event e;
    while(SDL_PollEvent(&e) != 0){
        if(e.type == SDL_QUIT){
            return false;
        }else if(e.type == SDL_KEYDOWN){
            switch(e.key.keysym.sym){
                case SDLK_ESCAPE: return false;
                default:break;
            }
        }
    }

    // draw
    // Clear the back buffer with a red background (parameters are R, G, B, A)
    glClearColor( 1.0, 0.0, 0.0, 1.0 );
    glClear(GL_COLOR_BUFFER_BIT);

    //glUseProgram(programID);
    shader1.use();
    drawVertexArray( 6, vbID_1 );

    SDL_GL_SwapWindow(window);
    return true;
};

int main(string[] args){

    // init SDL/GL window and context
    initDerelictSDL2GL3();
    scope(exit){ closeDerelictSDL2GL3(); }
    window1 = new WindowSDL( "firstTriangle", 640, 480 );

    // init GL objects
    vertexArrayID = makeVertexArrays();
    //vbID_1      = makeVertexBuffer( [ Vec3f(-1, -1, 0), Vec3f(1, -1, 0), Vec3f(0, 1, 0) ] );
    vbID_1        = makeVertexBuffer([ 
        Vec3f(-1, -1, 0), Vec3f(1, -1, 0), Vec3f(-1, 1, 0), 
        Vec3f(1, -1, 0), Vec3f(-1, 1, 0), Vec3f( 1,  1, 0),
    ]);
    //programID     = compileShaders(vertexShaderSrc, fragmentShaderSrc);

    shader1 = new Shader( vertexShaderSrc, fragmentShaderSrc, null );

    // per-frame  update / draw loop
    for(;;){ if( !frame(window1.window) ){ window1.destroy(); return 0; } }

}
