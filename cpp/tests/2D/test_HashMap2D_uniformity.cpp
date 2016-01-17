
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

class Histogram{
	public:
	double vmin,vmax;
	double step,invStep,range;
	int nbins;
	int * bins;

	void init( int nbins_, double vmin_, double vmax_ ){
		vmin  = vmin_; vmax = vmax_; nbins = nbins_;
		range = vmax - vmin;
		step  = range/nbins;
		invStep = 1/step;
		bins  = new int[nbins+1];
		for(int i=0;i<=nbins;i++){ bins[i]=0; }
	};

	int insert( double val ){
		int ibin = (val - vmin)*invStep;
		if( (ibin>=0)&&(ibin<nbins) ){
			//printf( " %3.3ff %i  \n", val, ibin );
			bins[ibin]++;
			return ibin;
		};
		return -1;
	};

};



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

void drawHistogram( int nbins, int * bins, float hstep ){
	for( int i=0; i<nbins; i++ ){
		Draw2D::drawLine( {i, 30}, {i,30+bins[i]*hstep} );
		//printf( " drawHistogram %i :  %i \n", i, bins[i]    );
	}
}

void printHistogram( const Histogram& hist ){
	printf( "printHistogram %i\n", hist.nbins );
	for( int i=0; i<hist.nbins; i++ ){
		printf( " drawHistogram %i %f:  %i \n", i, hist.vmin+i*hist.step, hist.bins[i]    );
	}
}

void printHashMap( const HashMap2D<Vec2d>& map ){
	int n = map.capacity;
	for( int i=0; i<n; i++ ){
		ULONG bucket = map.fields[i].bucket;
		UHALF ibx,iby;
		map.unfoldBucket( bucket, ibx,iby);
		if( map.fields[i].object == NULL ){
			printf( "field %03i  :  %03i (%i,%i) NULL \n", i, map.fields[i].n, ibx, iby );
		}else{
			Vec2d* p = map.fields[i].object;
			printf( "field %03i  :  %03i (%i,%i) (%3.3f,%3.3f) \n", i, map.fields[i].n, ibx, iby,  p->x, p->y  );
			//double x,y;
			//map.unfoldBucket( bucket, x, y );
			//Draw2D::drawRectangle( (float)x, (float)y, (float)(x+map.step), (float)(y+map.step), false );
			//printf( " (%3.3f,%3.3f) (%3.3f,%3.3f) \n", x, y,   p->x, p->y );
		}
	}
}

void drawHashMap( const HashMap2D<Vec2d>& map ){
	int n = map.capacity;
	for( int i=0; i<n; i++ ){
		ULONG bucket = map.fields[i].bucket;
		UHALF ibx,iby;
		map.unfoldBucket( bucket, ibx,iby);
		if( map.fields[i].object != NULL ){
			Vec2d* p = map.fields[i].object;
			double x,y;
			map.unfoldBucket( bucket, x, y );
			Draw2D::drawRectangle( (float)x, (float)y, (float)(x+map.step), (float)(y+map.step), false );
		}
	}
}

void testMapIndexing( const HashMap2D<Vec2d>& map, double x, double y ){
	double x_,y_;
	UHALF ix,iy,ix_,iy_;
	ix = map.getIx( x );
	iy = map.getIy( y );
	ULONG bucket = map.getBucket( ix, iy );
	map.unfoldBucket( bucket, ix_, iy_ );
	map.unfoldBucket( bucket, x_, y_ );
	printf( " (%3.3f,%3.3f), (%i,%i), %i, (%i,%i), (%3.3f,%3.3f) \n", x,y,  ix,iy, bucket, ix_,iy_,  x_,y_ );
}



// ======================  TestApp

#include "AppSDL2OGL.h"
class TestApp : public AppSDL2OGL {
	public:
	int npoints;
	Vec2d * points;
	HashMap2D<Vec2d> map;
	Histogram        hist;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
	TestApp( int& id, int WIDTH_, int HEIGHT_ );

};

TestApp::TestApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

	map.init( 0.5f, 16 );
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );
	int nside = 120;
	npoints = 4*nside*nside;
	points  = new Vec2d[npoints];
	int i = 0;
	for( int iy=-nside+1; iy<nside; iy++ ){
		for( int ix=-nside+1; ix<nside; ix++ ){
			i++;
			points[ i ].set( ( ix + randf() ) * map.step, ( iy + randf() ) * map.step );
			map.insertNoTest( &(points[i]), points[i].x, points[i].y  );
			//printf( " insering (%i,%i) %i (%3.3f,%3.3f) \n", ix, iy, i, points[i].x, points[i].y );
		};
	};
	printf( "map: %i %i %i %i \n", map.power, map.mask, map.capacity, map.filled );

	hist.init( 20, 0, 20 );
	for( int i=0; i<map.capacity; i++ ){
		hist.insert( map.fields[i].n + 0.5 );
	}
	printf( "now comes printHistogram( hist );\n" );
	printHistogram( hist );
}

void TestApp::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//Draw2D::drawPoints( npoints, points );

	//printHashMap( map );
	//drawHashMap( map );

	Vec2d* out[1<<16];
	UINT   nbucket  = map.getBucketObjects( mouse_begin_x, mouse_begin_y, &(out[0]) );
	//printf( " mouse (%f,%f) nbucket %i filled %i \n", mouse_begin_x, mouse_begin_y, nbucket, map.filled );
	for( int i=0; i<nbucket; i++ ){
		Draw2D::drawPointCross_d( *(out[i]), 0.3 );
	}

	// find points inside screen
	float fWIDTH = 0.98*10*ASPECT_RATIO;
	float fHEIHT = 0.98*10;
	float x0 = -camX0 - fWIDTH; float y0 = -camY0 - fHEIHT;
	float x1 = -camX0 + fWIDTH; float y1 = -camY0 + fHEIHT;
	printf( " x0,y0 (%3.3f,%3.3f) x1,y1 (%3.3f,%3.3f) %f\n", x0,y0, x1,y1, zoom );
	Draw2D::drawRectangle( x0, y0, x1, y1, false );
	UINT nfound = map.getObjectsInRect( x0, x1, y0, y1, out );
	glBegin(GL_POINTS);
	for( int i=0; i<nfound; i++ ){
		glVertex3f( (float) out[i]->x, (float)out[i]->y, 0.0f );
	}
	glEnd();

	//STOP = true;

	drawHistogram( hist.nbins+1, hist.bins, 3 );
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
















