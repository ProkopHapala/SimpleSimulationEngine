
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

#include "TerrainCubic.h"
#include "TiledView.h"

// ======================  TestApp

class TestAppTerrainCubic : public AppSDL2OGL, public TiledView {
	public:
	int npoints;
	Vec2d * points;
	TerrainCubic terrain;
	int viewList;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    //virtual void eventHandling( const SDL_Event& event );
	virtual int tileToList( float x0, float y0, float x1, float y1 );

	TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTerrainCubic::TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    terrain.init( 100, 100, 1.0 );
    terrain.x0 -= 50;
    terrain.y0 -= 50;
    terrain.allocate( );
    terrain.generateRandom( 0.0, 1.0 );

    //viewList = tileToList( 3, 3, 8, 8  );

    TiledView::init( 6, 6 ); // initialize parent
    //x0 -= 1000;
    //y0 -= 1000;
    tiles    = new int[ nxy ];

    TiledView::renderAll( -10, -10, 10, 10 );
    //exit(0);

}

void TestAppTerrainCubic::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    //terrain.renderRect( 3, 3, 8, 8, 30 );
    //if( viewList != 0 ) glCallList( viewList );

    //TiledView::render(  3, 3, 8, 8 );

    float szview = 4.0f;
    TiledView::draw(  mouse_begin_x-szview, mouse_begin_y-szview, mouse_begin_x+szview, mouse_begin_y+szview );

    glColor3f(0.2f,0.2f,0.9f); Draw2D::drawRectangle( -10, -10, 10, 10, false );
    glColor3f(0.9f,0.2f,0.9f); Draw2D::drawRectangle( x0, y0, x0+nx*step, y0+ny*step, false );
    //exit(0);

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

/*
void TestAppTerrainCubic::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
                case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
                case SDLK_KP_6:  TiledView::shiftRender( +1,  0 ); break;
                case SDLK_KP_4:  TiledView::shiftRender( -1,  0 ); break;
                case SDLK_KP_8:  TiledView::shiftRender(  0, +1 ); break;
                case SDLK_KP_2:  TiledView::shiftRender(  0, -1 ); break;
            }
            break;
        case SDL_QUIT: quit(); break;
    };
};
*/

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
















