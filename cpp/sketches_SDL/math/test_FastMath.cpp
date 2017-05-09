
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

inline double floor_exact( double x ){
    int ix   = floor( x );
    return x - ix;
}

inline double floor_fast( double x ){
    //int ix = static_cast<int>( x );
    //if( x<0 ) ix--;
    int ix = fastFloor( x );
    return x - ix;
}

inline double floor_veryfast( double x ){
    int ix = static_cast<int>( x + 10000 ) - 10000;
    return x - ix;
}

inline double floor_veryfast_c( double x ){
    int ix = (int)( x + 10000 ) - 10000;
    return x - ix;
}

inline double floor_justcast( double x ){
    int ix = static_cast<int>( x );
    return x - ix;
}

inline double floor_justcast_c( double x ){
    int ix = (int)( x );
    return x - ix;
}

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


    //glColor3f(0.9f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &floor_justcast    );
    //glColor3f(0.5f,0.2f,0.5f); Draw2D::drawFunc( -3, 3, 1000, &floor_justcast_c  );


    double xmin =  -0.5;
    double xmax =   1.5;
    int ncall   = 100000000;
    //int ncall   = 100;
    int narr = 10000;
    double * arr = new double[ narr ];
    genRandomArray( narr, arr, -1000.0, +1000.0 );

    printf("\n");

    glColor3f(0.2f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &floor_exact       );
    glColor3f(0.2f,0.2f,0.9f); Draw2D::drawFunc( -3, 3, 1000, &floor_fast        );
    glColor3f(0.9f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &floor_veryfast    );
    glColor3f(0.5f,0.2f,0.2f); Draw2D::drawFunc( -3, 3, 1000, &floor_veryfast_c  );
    SPEED_TEST_FUNC( "void       ",                  , xmin, xmax, ncall );
    SPEED_TEST_FUNC( "exact      ", floor_exact      , xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "fast       ", floor_fast       , xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "veryfast   ", floor_veryfast   , xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "veryfast_c ", floor_veryfast_c , xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "justcast   ", floor_justcast   , xmin, xmax, ncall  );
    SPEED_TEST_FUNC( "justcast_c ", floor_justcast_c , xmin, xmax, ncall  );
    int m = ncall/narr;
    printf("%i \n", m);
    SPEED_TEST_FUNC_ARRAY( "void       ",                  , arr,  narr,  m );
    SPEED_TEST_FUNC_ARRAY( "exact      ", floor_exact      , arr,  narr,  m );
    SPEED_TEST_FUNC_ARRAY( "fast       ", floor_fast       , arr,  narr,  m );
    SPEED_TEST_FUNC_ARRAY( "veryfast   ", floor_veryfast   , arr,  narr,  m );
    SPEED_TEST_FUNC_ARRAY( "veryfast_c ", floor_veryfast_c , arr,  narr,  m );
    SPEED_TEST_FUNC_ARRAY( "justcast   ", floor_justcast   , arr,  narr,  m );
    SPEED_TEST_FUNC_ARRAY( "justcast_c ", floor_justcast_c , arr,  narr,  m );

    //SPEED_TEST_FUNC_ARRAY( "acos       ",                  , arr,  narr,  m );

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
















