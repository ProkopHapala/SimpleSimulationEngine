
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "TileTree2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

// ======================  TestApp

/*
class GridCell{
	public:
	int value;
	inline bool isEmpty(){ return value == 0; };
	GridCell(){ value = 0; }
};
*/

#define CELL_TYPE int

class TestAppTileTree2D : public AppSDL2OGL {
	public:
	//TileTree2D<GridCell,2,3,5> map;
	TileTree2D<CELL_TYPE,2,3,5> map;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
	TestAppTileTree2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTileTree2D::TestAppTileTree2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
	printf( " %i %i   %i %i %i \n", map.nsub, map.sub_mask, map.ntotx, map.ntoty, map.ntotxy );

	auto ptile = map.tiles[0];

	printf( " >> %i %i       \n", ptile->n, ptile->n2 );
	printf( " >> %i %i %i %i \n", ptile->index2D(1,0), ptile->index2D(0,1), ptile->index2D(ptile->n-1,0), ptile->index2D(0,ptile->n-1) );

	( *map.getValidPointer( 2, 3,  0 ) ) = 10; 
	( *map.getValidPointer( 4, 5,  0 ) ) = 11;
	( *map.getValidPointer( 6, 7,  0 ) ) = 12;
	
}

void TestAppTileTree2D::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	for( int iy=0; iy<map.ntoty;  iy++ ){
		for( int ix=0; ix<map.ntotx;  ix++ ){
			//GridCell cell;
			//map.getBare( ix, iy, cell );
			//if( cell.value > 0 ) printf( " %i %i = %i \n", ix, iy, cell  );
			CELL_TYPE* pcell = map.getPointer( ix, iy );
			if( pcell != NULL ){
				CELL_TYPE cell = *pcell;
				//if( cell > 0 ) printf( " %i %i = %i \n", ix, iy, cell );
				//printf( " %i %i = %i \n", ix, iy, cell );
			}
		}
	}

	exit(0);
	//Draw2D::drawPoints( npoints, points );

};

void TestAppTileTree2D::drawHUD(){
	//Draw2D::drawRectangle( 10,10, 100,200, true );
	//drawHashMapFilling( map );
}

// ===================== MAIN

TestAppTileTree2D * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppTileTree2D( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















