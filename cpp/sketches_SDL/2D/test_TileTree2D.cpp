/// @file @brief This demo showcases `TileTree2D.h`, a 2D spatial data structure optimized for managing and querying rectangular tiles. It displays a grid of tiles, highlighting the tile currently under the mouse cursor, illustrating its use for efficient spatial indexing of rectangular regions.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "TileTree2D.h"
#include "ArrayMap2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

// ======================  TestApp
#define CELL_TYPE int

class TestAppTileTree2D : public AppSDL2OGL {
	public:

	//TileTree2D<GridCell,2,3,5> map;
	TileTree2D_d<CELL_TYPE,3,200,300> map1;
    ArrayMap2D  <CELL_TYPE,200*8,300*8> map2;

    bool mouse_left = false;
    bool mouse_right = false;

    long t1,t12;

	// ---- function declarations

	void         drawMap( );

	void   set_speed ( int ntry, double xmin, double ymin, double xmax, double ymax );
	double get_speed ( int ntry, double xmin, double ymin, double xmax, double ymax );

    void   set_speed_ref( int ntry, double xmin, double ymin, double xmax, double ymax );
	double get_speed_ref( int ntry, double xmin, double ymin, double xmax, double ymax );

	virtual void draw   ();
	void eventHandling( const SDL_Event& event );
	TestAppTileTree2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTileTree2D::TestAppTileTree2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

	printf( " %i %i   %i %i %i \n", map1.nsub, map1.sub_mask, map1.ntotx, map1.ntoty, map1.ntotxy );

    double step = 0.5;
    map1.setSpacing( step, step );

    int n_try = 100000000;

	t1 = getCPUticks();
    set_speed( n_try, -50.0, -50.0, 50.0, 50.0 );
	t12 = getCPUticks() - t1;
	printf( " set_speed %3.3f ticks/iter | %3.3f Mticks iterations %i \n", t12/((double)n_try), t12*1.0e-6, n_try );

	t1 = getCPUticks();
    get_speed( n_try, -50.0, -50.0, 50.0, 50.0 );
	t12 = getCPUticks() - t1;
	printf( " get_speed %3.3f ticks/iter | %3.3f Mticks iterations %i \n", t12/((double)n_try), t12*1.0e-6, n_try );

    map2.setSpacing( step, step );

	t1 = getCPUticks();
    set_speed_ref( n_try, -50.0, -50.0, 50.0, 50.0 );
	t12 = getCPUticks() - t1;
	printf( " set_speed_ref %3.3f ticks/iter | %3.3f Mticks iterations %i \n", t12/((double)n_try), t12*1.0e-6, n_try );

	t1 = getCPUticks();
    get_speed_ref( n_try, -50.0, -50.0, 50.0, 50.0 );
	t12 = getCPUticks() - t1;
	printf( " get_speed_ref %3.3f ticks/iter | %3.3f Mticks iterations %i \n", t12/((double)n_try), t12*1.0e-6, n_try );


}

void TestAppTileTree2D::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable( GL_DEPTH );

    if      ( mouse_left  ){ ( *map1.getValidPointer_d( mouse_begin_x, mouse_begin_y,  0 ) ) = 255; }
    else if ( mouse_right ){ ( *map1.getValidPointer_d( mouse_begin_x, mouse_begin_y,  0 ) ) = 0;   }

    drawMap();

};

void TestAppTileTree2D::drawMap(){
    for( int iy=0; iy<map1.ny; iy++ ){
		for( int ix=0; ix<map1.nx; ix++ ){
		    auto ptile = map1.tiles[ map1.isup2D( ix, iy ) ];
            if( ptile != NULL ){
                double xi = map1.getX( (ix<<map1.power) );
                double yi = map1.getY( (iy<<map1.power) );
                glColor3f( 0, 0, 0 );
                Draw2D::z_layer = 0.0f;
                Draw2D::drawRectangle_d( {xi,yi}, {xi+map1.xStep*map1.nsub,yi+map1.yStep*map1.nsub}, true );
                Draw2D::z_layer = 1.0f;
                for( int jy=0; jy<map1.nsub;  jy++ ){
                    double y = map1.getY( (iy<<map1.power) + jy );
                    for( int jx=0; jx<map1.nsub; jx++ ){
                        double x = map1.getX( (ix<<map1.power) + jx );
                        CELL_TYPE* pcell= ptile->getPointer( jx, jy );
                        if( *pcell != 0 ){
                            float c = *pcell/255.0f;
                            glColor3f( c, c, c );
                            Draw2D::drawRectangle_d( {x,y}, {x+map1.xStep,y+map1.yStep}, true );
                        }
                    }
                }
            }
		}
    }
}

void TestAppTileTree2D::set_speed( int ntry, double xmin, double ymin, double xmax, double ymax ){
    double sum = 0;
    double xsc = ( xmax - xmin )/65536.0d;
    double ysc = ( ymax - ymin )/65536.0d;
    int h = rand_hash2( 2147483647 );
    for( int i=0; i<ntry; i++ ){
        h = rand_hash( h );
        double x = ( ( h&0xFFFF          ) ) * xsc + xmin;
        double y = ( ((h&0xFFFF0000)>>16 ) ) * ysc + ymin;
        int val  =   (h&0xFF0000)>>16 ;

        CELL_TYPE* pcell = map1.getValidPointer_d( x, y, 0 );
        (*pcell) = val;

        sum += x + y + val;
    }
    printf( " sum %e \n", sum );
}

double TestAppTileTree2D::get_speed( int ntry, double xmin, double ymin, double xmax, double ymax ){
    double sum = 0 ;
    double xsc = ( xmax - xmin )/65536.0d;
    double ysc = ( ymax - ymin )/65536.0d;
    int h = rand_hash2( 2147483647 );
    for( int i=0; i<ntry; i++ ){
        h = rand_hash( h );
        double x = ( ( h&0xFFFF          ) ) * xsc + xmin;
        double y = ( ((h&0xFFFF0000)>>16 ) ) * ysc + ymin;
        int val  =   (h&0xFF0000)>>16 ;

        CELL_TYPE* pcell = map1.getPointer_d( x, y );
        if( pcell != NULL ){
            sum += *pcell;
        }

        sum += x + y + val;
    }
    printf( " sum %e \n", sum );
    return sum;
}

void TestAppTileTree2D::set_speed_ref( int ntry, double xmin, double ymin, double xmax, double ymax ){
    double sum = 0;
    double xsc = ( xmax - xmin )/65536.0d;
    double ysc = ( ymax - ymin )/65536.0d;
    int h = rand_hash2( 2147483647 );
    for( int i=0; i<ntry; i++ ){
        h = rand_hash( h );
        double x = ( ( h&0xFFFF          ) ) * xsc + xmin;
        double y = ( ((h&0xFFFF0000)>>16 ) ) * ysc + ymin;
        int val  =   (h&0xFF0000)>>16 ;

        //CELL_TYPE* pcell = map2.getPointer( x, y );
        //(*pcell) = val;

        map2.cells[ map2.isup2D( map2.getIx(x), map2.getIy(y) ) ] = val;

        sum += x + y + val;
    }
    printf( " sum %e \n", sum );
}

double TestAppTileTree2D::get_speed_ref( int ntry, double xmin, double ymin, double xmax, double ymax ){
    double sum = 0 ;
    double xsc = ( xmax - xmin )/65536.0d;
    double ysc = ( ymax - ymin )/65536.0d;
    int h = rand_hash2( 2147483647 );
    for( int i=0; i<ntry; i++ ){
        h = rand_hash( h );
        double x = ( ( h&0xFFFF          ) ) * xsc + xmin;
        double y = ( ((h&0xFFFF0000)>>16 ) ) * ysc + ymin;
        int val  =   (h&0xFF0000)>>16 ;

        //CELL_TYPE* pcell = map2.getPointer( x, y );
        //if( pcell != NULL ){
        //    sum += *pcell;
        //}

        sum += map2.cells[ map2.isup2D( map2.getIx(x), map2.getIy(y) ) ];

        sum += x + y + val;
    }
    printf( " sum %e \n", sum );
    return sum;
}

void TestAppTileTree2D::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:   mouse_left  = true; break;
                case SDL_BUTTON_RIGHT:  mouse_right = true; break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:   mouse_left  = false; break;
                case SDL_BUTTON_RIGHT:  mouse_right = false; break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
};

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
