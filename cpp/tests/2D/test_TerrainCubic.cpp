
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "drawMath2D.h"
#include "AppSDL2OGL.h"

#include "TerrainCubic.h"
#include "TiledView.h"

// ======================  TestApp

class TestAppTerrainCubic : public AppSDL2OGL, public TiledView {
	public:
	int npoints;
	Vec2d * points;
	TerrainCubic terrain;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();

	virtual int tileToList( float x0, float y0, float x1, float y1 );

	TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTerrainCubic::TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    terrain.init( 100, 100, 1.0 );
    terrain.allocate( );
    terrain.generateRandom( 0.0, 1.0 );

    Map2D::init( 10, 10, 30 ); // initialize parent
    tiles    = new int[ nxy ];

}

void TestAppTerrainCubic::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    terrain.renderRect( 3, 3, 8, 8, 30 );

    //exit(0);

};


 int TestAppTerrainCubic::tileToList( float x0, float y0, float x1, float y1 ){
    terrain.renderRect( x0, y0, x1, y1,  30 );
 }

void TestAppTerrainCubic::drawHUD(){
}

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
















