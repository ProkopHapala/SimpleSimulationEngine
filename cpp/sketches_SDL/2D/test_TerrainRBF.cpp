/// @file @brief This demo showcases `TerrainRBF.h`, a technique for generating or interpolating 2D terrain using Radial Basis Functions. It visualizes the smoothly interpolated terrain surface, allowing the user to interactively define or adjust the RBF control points that shape the landscape.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "TerrainRBF.h"
#include "TiledView.h"

// ======================  TestApp

class TestAppTerrainCubic : public AppSDL2OGL, public TiledView {
	public:
	int npoints;
	RBF2D * rbfs;
	TerrainRBF terrain;
	int viewList;

	bool clicked = false;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	virtual int tileToList( float x0, float y0, float x1, float y1 );

	TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTerrainCubic::TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    terrain.init( 2.0d, 16, 256 );
    terrain.generateRandom(  -10.0,-10.0, 10.0,10.0,  1.0,2.0,   0.2,1.0, 15 );

    //viewList = tileToList( -10, -10, 10, 10  );

    TiledView::init( 6, 6 ); // initialize parent
    //x0 -= 1000;
    //y0 -= 1000;
    //tiles    = new int[ nxy ];

    TiledView::renderAll( -10, -10, 10, 10 );
    //exit(0);

}

void TestAppTerrainCubic::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    //terrain.renderRect( 3, 3, 8, 8, 30 );
    //if( viewList != 0 ) glCallList( viewList );

    //exit(0);

    //TiledView::render(  3, 3, 8, 8 );

    float szview = 4.0f;
    TiledView::draw(  mouse_begin_x-szview, mouse_begin_y-szview, mouse_begin_x+szview, mouse_begin_y+szview );

    glColor3f(0.2f,0.2f,0.9f); Draw2D::drawRectangle( -10, -10, 10, 10, false );
    glColor3f(0.9f,0.2f,0.9f); Draw2D::drawRectangle( x0, y0, x0+nx*step, y0+ny*step, false );
    //exit(0);


    for( int i=0; i<terrain.nrbfs; i++ ){
        //printf( " %i  (%3.3f,%3.3f) \n", i, terrain.rbfs[i].pos.x, terrain.rbfs[i].pos.y  );
        Draw2D::drawPointCross_d( terrain.rbfs[i].pos, 0.1 );
        Draw2D::drawCircle_d( terrain.rbfs[i].pos, terrain.rbfs[i].rmax, 16, false  );
    }

    if(clicked){
        double val = terrain.getVal( mouse_begin_x, mouse_begin_y );
        printf( " terrain( %3.3f, %3.3f ) : %f \n", mouse_begin_x, mouse_begin_y, val );
        clicked = false;
    }
};


int TestAppTerrainCubic::tileToList( float x0, float y0, float x1, float y1 ){
    int ilist=glGenLists(1);
    glNewList( ilist, GL_COMPILE );
        terrain.renderRect( x0, y0, x1, y1, 9 );
        glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
    glEndList();
    return ilist;
}

void TestAppTerrainCubic::drawHUD(){
}


void TestAppTerrainCubic::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    clicked = true;
                break;
            }
            /*
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
            */
    };
    AppSDL2OGL::eventHandling( event );
};

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
