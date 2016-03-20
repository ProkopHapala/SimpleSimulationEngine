
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include "testUtils.h"

#include "tresholdFunctions.h"


// ======================  TestApp

class TestAppTerrainCubic : public AppSDL2OGL {
	public:

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();

	TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTerrainCubic::TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    //exit(0);
}

void TestAppTerrainCubic::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    //glColor3f(0.9f,0.9f,0.9f); Draw2D::drawGrid( -10, -10, 10, 10, 0.1, 0.1 );
    glColor3f(0.9f,0.9f,0.9f); Draw2D::drawGrid( -5, -5, 5, 5, 0.5, 0.5 );
    glColor3f(0.8f,0.8f,0.8f); Draw2D::drawGrid( -5, -5, 5, 5, 1.0, 1.0 );

    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p0i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p0  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p0d );

    glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p1i );
    glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p1  );
    glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p1d );
    glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p1i );
    glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p1  );

    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p2i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p2  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p2d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p2i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p2  );

    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p3i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p3  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p3d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p3i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p3  );


    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p4i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p4  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p4d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p4i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p4  );


    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p5i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p5  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::p5d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p5i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::p5  );


    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::r1i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::r1  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::r1d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::r1i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::r1  );

    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::r2i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::r2  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::r2d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::r2i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::r2  );

    //glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::root2i );
    //glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::root2  );
    //glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &Treshold::root2d );
    //glColor3f(0.2f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::root2i );
    //glColor3f(0.9f,0.9f,0.2f); Draw2D::drawFuncDeriv( -3, 3, 0.01, 1000, &Treshold::root2  );


    double xmin =  -0.5;
    double xmax =   1.5;
    int ncall   = 10000000;

    printf("\n");
    SPEED_TEST_FUNC( "void  ",              , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "void  ",              , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "void  ",              , xmin, xmax, ncall );

    printf("\n");
    SPEED_TEST_FUNC( "p1i   ", Treshold::p1i, xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "p1    ", Treshold::p1 , xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "p1d   ", Treshold::p1d, xmin, xmax, ncall  );

    printf("\n");
    SPEED_TEST_FUNC( "p2i   ", Treshold::p2i, xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "p2    ", Treshold::p2 , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "p2d   ", Treshold::p2d, xmin, xmax, ncall );

    printf("\n");
    SPEED_TEST_FUNC( "p3i   ", Treshold::p3i, xmin, xmax, ncall );
    SPEED_TEST_FUNC( "p3    ", Treshold::p3 , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "p3d   ", Treshold::p3d, xmin, xmax, ncall );

    printf("\n");
    SPEED_TEST_FUNC( "p5i   ", Treshold::p5i, xmin, xmax, ncall );
    SPEED_TEST_FUNC( "p5    ", Treshold::p5 , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "p5d   ", Treshold::p5d, xmin, xmax, ncall );

    printf("\n");
    SPEED_TEST_FUNC( "r1i   ", Treshold::r1i, xmin, xmax, ncall );
    SPEED_TEST_FUNC( "r1    ", Treshold::r1 , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "r1d   ", Treshold::r1d, xmin, xmax, ncall );

    printf("\n");
    SPEED_TEST_FUNC( "r2i   ", Treshold::r2i, xmin, xmax, ncall );
    SPEED_TEST_FUNC( "r2    ", Treshold::r2 , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "r2d   ", Treshold::r2d, xmin, xmax, ncall );

    printf("\n");
    SPEED_TEST_FUNC( "root2i", Treshold::root2i, xmin, xmax, ncall );
    SPEED_TEST_FUNC( "root2 ", Treshold::root2 , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "root2d", Treshold::root2d, xmin, xmax, ncall );

    STOP = true;

};

void TestAppTerrainCubic::drawHUD(){}


// ===================== MAIN

TestAppTerrainCubic * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppTerrainCubic( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















