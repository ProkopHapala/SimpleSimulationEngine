
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

#include "AppSDL2OGL.h"
//#include "testUtils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "Formation.h"
#include "FormationWorld.h"

class FormationTacticsApp : public AppSDL2OGL, public TiledView {
	public:
    FormationWorld world;

    Formation * currentFormation = NULL;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	//void pickParticle( Particle2D*& picked );

	virtual int tileToList( float x0, float y0, float x1, float y1 );

	FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ );

};

FormationTacticsApp::FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    world.init();
    currentFormation = world.formations[0];

    TiledView::init( 6, 6 );
    tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );
}

void FormationTacticsApp::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	world.update( );

    float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );


    glColor3f( 0.8f, 0.1f, 0.1f );
    for( Formation* f : world.formations ){
        //printf( " f %i \n", f  );
        if( f != NULL ) f->render( );
    }

    //if(frameCount > 10) STOP = true;
};

int FormationTacticsApp::tileToList( float x0, float y0, float x1, float y1 ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
		world.terrain.renderRect( x0, y0, x1, y1, 31 );
		//glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
	glEndList();
	return ilist;
}

void FormationTacticsApp::drawHUD(){}

void FormationTacticsApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    if( currentFormation != NULL ) currentFormation->setTarget( { mouse_begin_x, mouse_begin_y } );
                break;
            }
            break;
            /*
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
}

// ===================== MAIN

FormationTacticsApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new FormationTacticsApp( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















