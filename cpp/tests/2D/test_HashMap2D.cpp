
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
#include "HashMap2D.h"
#include "drawMath2D.h"

// ==================== HashMap Debug utils

void drawHasMapTile( const HashMap2D<Vec2d>& map, ULONG index ){
	double x,y;
	map.unfoldBucket( index, x, y );
	Draw2D::drawRectangle( (float)x, (float)y, (float)(x+map.step), (float)(y+map.step), false );
	//printf( " drawHasMapTile %f %f %f %f \n", (float)x, (float)y, (float)(x+map.step), (float)(y+map.step) );
}


void drawHashMapFilling( const HashMap<Vec2d>& map ){
	int n = map.capacity;
	for( int i=0; i<n; i++ ){
		int ni = map.fields[i].n;
		Draw2D::drawLine( {i, 10}, {i,10+ni*3} );
	}
	//map.unfoldBucket( index, x, y );
}


void printHashMap( const HashMap2D<Vec2d>& map ){
	int n = map.capacity;
	//printf( "======================= \n");
	for( int i=0; i<n; i++ ){
		ULONG bucket = map.fields[i].bucket;
		UHALF ibx,iby;
		map.unfoldBucket( bucket, ibx,iby);
		if( map.fields[i].object == NULL ){
			//printf( "field %03i  :  %03i (%i,%i) NULL \n", i, map.fields[i].n, ibx, iby );
		}else{
			Vec2d* p = map.fields[i].object;
			//printf( "field %03i  :  %03i (%i,%i) (%3.3f,%3.3f) \n", i, map.fields[i].n, ibx, iby,  p->x, p->y  );
			//drawHasMapTile( map, bucket );
			double x,y;
			map.unfoldBucket( bucket, x, y );
			Draw2D::drawRectangle( (float)x, (float)y, (float)(x+map.step), (float)(y+map.step), false );
			//printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) \n", x, y,   p->x, p->y );
		}
	}
}

// ======================  TestApp

#include "AppSDL2OGL.h"
class TestApp : public AppSDL2OGL {
	public:	
	int npoints;
	Vec2d * points;
	HashMap2D<Vec2d> map; 

	// ---- function declarations 

	virtual void draw   ();
	virtual void drawHUD();
	TestApp( int& id, int WIDTH_, int HEIGHT_ );

};

TestApp::TestApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
	map.init( 1.0f, 6 );
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );
	npoints = 30;
	points  = new Vec2d[npoints];
	//randomize( 289 );
	for( int i=0; i<npoints; i++ ){ 
		
		double x  = randf(-5.0,5.0);
		double y  = randf(-5.0,5.0); 

/*
		double x_,y_;
		UHALF ix,iy,ix_,iy_;
		ix = map.getIx( x );
		iy = map.getIy( y );
		ULONG bucket = map.getBucket( ix, iy );
		map.unfoldBucket( bucket, ix_, iy_ );
		map.unfoldBucket( bucket, x_, y_ );
		printf( " (%3.3f,%3.3f), (%i,%i), %i, (%i,%i), (%3.3f,%3.3f) \n", x,y,  ix,iy, bucket, ix_,iy_,  x_,y_ );
*/		
		points[i].set( x, y ); 
		int index = map.insertNoTest( &points[i], points[i].x, points[i].y  );
		//int index = map.insertIfNew( &(points[i]), points[i].x, points[i].y  );
		//printf( " inserting %i-th point to index %i \n", i, index  );
		//printHashMap( map );
	};
	//exit(0);
}

void TestApp::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	Draw2D::drawPoints( npoints, points );

	Vec2d* out[16];
	printHashMap( map );
	UINT   nbucket  = map.getBucketObjects( mouse_begin_x, mouse_begin_y, &(out[0]) );
	printf( " mouse (%f,%f) nbucket %i filled %i \n", mouse_begin_x, mouse_begin_y, nbucket, map.filled );
	//exit(0);
	for( int i=0; i<nbucket; i++ ){
		Draw2D::drawPointCross_d( *(out[i]), 0.3 ); 
	}
	//STOP = true;
};

void TestApp::drawHUD(){
	//Draw2D::drawRectangle( 10,10, 100,200, true );
	drawHashMapFilling( map );
}

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
















