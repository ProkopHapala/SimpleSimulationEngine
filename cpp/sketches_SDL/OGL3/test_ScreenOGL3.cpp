#include <stdlib.h>
#include <stdio.h>


#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include <SDL2/SDL.h>

#include <fastmath.h>
#include <Vec2.h>
#include <Vec3.h>
#include <Mat3.h>
#include <quaternion.h>
#include <raytrace.h>
//#include <Body.h>

#include "Shader.h"
#include "GLObject.h"
#include "SceneNode.h"
#include "ScreenSDL2OGL3.h"

// ========== functions

class TestAppScreenOGL3{
    public:

    int frameCount = 0;
    bool STOP = false;

    SceneNode3D    * thisNode;
    ScreenSDL2OGL3 * thisScreen;

    //static const int nBodies = 16;
    //PointBody bodies[nBodies];      // this would require to decouple "class Body" from SDL2OGL

    void inputHanding ();
    void init();
    void draw();
    void loop( int nframes );
    void quit();

    TestAppScreenOGL3();

};

TestAppScreenOGL3::TestAppScreenOGL3(){

    // ------------- object
/*
	object1 = new GLObject( );
	object1->nVert   = 4;
	object1->vertDim = 2;
	object1->vertexes = &vertexes[0][0];
	object1->init();
*/

    // ------------- shader

/*
	shader1=new Shader();
	//shader1->init( "shaders/plain_vert.c", "shaders/sphere_frag.c" );
	//shader1->init( "shaders/afine2D_vert.c", "shaders/sphere_frag.c" );
	shader1->init( "shaders/plain_vert.c", "shaders/texture_frag.c" );

    resolution[0] = (float)WIDTH;
    resolution[1] = (float)HEIGHT;

	glUseProgram(shader1->shaderprogram);

    GLuint uloc;
    uloc = glGetUniformLocation( shader1->shaderprogram, "resolution" );	glUniform2fv(uloc, 1, resolution  );
*/

}

// FUNCTION ======	inputHanding
void TestAppScreenOGL3::inputHanding(){
	SDL_Event event;
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE ) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE  ) { STOP=!STOP; }
		}
		if( event.type == SDL_QUIT){ quit();  };
	}
	int mouseX, mouseY;
	SDL_GetMouseState( &mouseX, &mouseY );
}

void TestAppScreenOGL3::quit(){
	//glDeleteVertexArrays(1, &vao);
    //if( context != NULL ) SDL_GL_DeleteContext( context );
    //if( window  != NULL ) SDL_DestroyWindow   ( window  );
    SDL_Quit();
	exit(0);
};

void TestAppScreenOGL3::loop( int nframes ){
    for ( int iframe=1; iframe<nframes; iframe++)    {
 		if( ( !STOP )&&( thisScreen != NULL ) )
            thisScreen->draw();
		inputHanding();
		frameCount++;
        SDL_Delay(10);
    }
}


// ================== main


TestAppScreenOGL3 * app;

int main(int argc, char *argv[]){

    app = new TestAppScreenOGL3( );
    //app->init();
	//app->setup();
    app->loop( 1000000 );
    app->quit();
    return 0;
}

