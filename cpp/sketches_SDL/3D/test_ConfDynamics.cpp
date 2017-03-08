
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"
#include "RBodyConfDyn.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

// ======================  TestApp

class TestConfDynamics : public AppSDL2OGL_3D {
	public:

	RBodyConfDyn confWorld;
	int ivob;

	// ---- function declarations

	virtual void draw();

	TestConfDynamics( int& id, int WIDTH_, int HEIGHT_ );

};

TestConfDynamics::TestConfDynamics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

	ivob=glGenLists(1);
	glNewList( ivob, GL_COMPILE );
        glDisable ( GL_LIGHTING );
        Draw3D::drawAxis(0.2);
	glEndList();

	confWorld.init();

}

void TestConfDynamics::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	Draw3D::drawShape( *confWorld.pos, *confWorld.rot, ivob );

	glDisable ( GL_LIGHTING );
	//glCallList( point_cloud );
	//Draw3D::drawAxis ( 3.0f );

};

// ===================== MAIN

TestConfDynamics * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestConfDynamics( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















