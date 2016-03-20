
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "Vec2.h"
#include "geom2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

// ==================== Declaration

class TestApp : public AppSDL2OGL {
	public:
	int frameCount = 0;
	bool loopEnd,STOP;


	// ---- function declarations 

	//virtual void loop( int nframes );
	virtual void draw   ();
	//virtual void drawHUD();
	TestApp( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp::TestApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

}

void TestApp::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    Vec2d A,B;
    A.set( -1.0, -5.0 );
    B.set( +3.0, +5.0 );
	Line2d cutline; cutline.set( A, B );

    glColor3f( 0.9f, 0.2f, 0.9f ); Draw2D::drawPointCross_d( A, 0.3 );
    glColor3f( 0.2f, 0.9f, 0.9f ); Draw2D::drawPointCross_d( B, 0.3 );    
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine_d( A, B );

};

// ===================== MAIN

TestApp * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestApp( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















