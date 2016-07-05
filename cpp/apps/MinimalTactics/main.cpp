
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "UnitType.h"
#include "Unit.h"
#include "Faction.h"
#include "TacWorld.h"


class FormationTacticsApp : public AppSDL2OGL, public TiledView {
	public:
    TacWorld world;

    int formation_view_mode = 0;
    Unit    * currentUnit      = NULL;
    Faction * currentFaction   = NULL;
    int       ifaction = 0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );


	void debug_buffinsert( );

	//void pickParticle( Particle2D*& picked );

	virtual int tileToList( float x0, float y0, float x1, float y1 );

	FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ );

};

FormationTacticsApp::FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    world.init();

    currentFaction = world.factions[0];  printf( "currentFaction: %s\n", currentFaction->name );
    currentUnit    = currentFaction->units[0];

    TiledView::init( 6, 6 );
    tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );
}

void FormationTacticsApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );
    //long tComp = getCPUticks();
    world.update( );
    //tComp = getCPUticks() - tComp;

    long tDraw = getCPUticks();
    for( Unit* u : world.units ){
        if( (u!= NULL) ){
            // TODO : check if on screen
            if  ( u == currentUnit ){ u->render( u->faction->color );   }
            else                    { u->render( u->faction->color );   }
        }
    }
    tDraw = getCPUticks() - tDraw;

    if( currentUnit != 0 ){
        //glColor3f(1.0,0.0,1.0);
        glColor3f(0.0,1.0,0.0);
        Draw2D::drawCircle_d( currentUnit->pos, 0.5, 16, false );
        currentUnit->renderJob( currentUnit->faction->color );
    }


    //printf( " frame %i : %i %i   %4.3f %4.3f %4.3f\n", frameCount, world.nSoldiers, world.nSoldierInteractions,
    //                             tComp*1.0e-6, tComp/(double)world.nSoldiers, tComp/(double)world.nSoldierInteractions );

    //tTot = getCPUticks() - tTot;
    //printf( " frame %i : %i %i   %4.3f %4.3f %4.3f \n", frameCount, world.nSoldiers, world.nSoldierInteractions,
    //       tComp*1.0e-6, tDraw*1.0e-6, tTot*1.0e-6 );
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
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                //case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                //case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                //case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                //case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
                case   SDLK_n: ifaction++; if(ifaction>=world.factions.size()) ifaction=0; currentFaction = world.factions[ifaction]; printf("ifaction %i\n",ifaction); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    if( currentFaction != NULL ) currentUnit = currentFaction->getUnitAt( { mouse_begin_x, mouse_begin_y } );
                break;
                case SDL_BUTTON_RIGHT:
                    //printf( "left button pressed !!!! " );
                    if( currentUnit != NULL ){
                        int imin = world.getUnitAt( { mouse_begin_x, mouse_begin_y }, currentFaction );
                        if( imin > -1 ) {
                            printf( "target selected %i %i\n", imin, world.units[imin] );
                            currentUnit->setOpponent( world.units[imin] );
                        }else{
                            printf( "goal selected (%3.3f,%3.3f)\n", mouse_begin_x, mouse_begin_y );
                            currentUnit->setGoal  ( { mouse_begin_x, mouse_begin_y } );
                        }
                    }
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
    camStep = zoom*0.05;
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
















