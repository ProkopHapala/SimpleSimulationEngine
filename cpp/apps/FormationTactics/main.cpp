
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

#include "Formation.h"
#include "FormationWorld.h"

/*

Performance:

Method  soldiers    interactions    Mticks  ticks/soldier   ticks/interaction

noBuff  512         131072          2.269   4431.172        17.309
noBuff  1012        606568          10.398  10274.514       17.14
noBuff  8188        6121488         79.428  9700.575        12.975
noBuff  16384       12419072        210.705 12860.431       16.966
noBuff  32374       48482256        787.849 24335.853       16.250

Buff    512         12902           0.613   1196.438        47.479
Buff    1012        39874           1.799   1778.087        45.128
Buff    8164        252302          11.865  1453.294        47.02
Buff    16377       527631          24.428  1491.580        46.29
Buff    32393       1731907         63.651  1964.948        36.752

*/

class FormationTacticsApp : public AppSDL2OGL, public TiledView {
	public:
    FormationWorld world;

    int formation_view_mode = 0;
    Formation * currentFormation = NULL;
    Faction   * currentFaction   = NULL;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );


	void debug_buffinsert( );

	//void pickParticle( Particle2D*& picked );

	virtual int tileToList( float x0, float y0, float x1, float y1 );

	FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ );

};

void FormationTacticsApp::debug_buffinsert( ){
    if( currentFormation == NULL ) return;
    Formation * fi = currentFormation;
    double r = world.RmaxInteract;
    world.colruler.setup(  {fi->bbox.x0-3*r,fi->bbox.y0-3*r}, {2*r+0.1,2*r+0.1} );
    //world.colruler.setup(  {fi->bbox.x0-3*r,fi->bbox.y0-3*r}, {3*r,3*r} );
    double rf = r/world.colruler.step.x;
    //printf( " pos0 (%3.3f,%3.3f) step (%3.3f,%3.3f) \n", world.colruler.pos0, world.colruler.pos0.y, world.colruler.step.y, world.colruler.step.y );
    world.colbuf.clear( );
    Vec2d pos,dipos; Vec2i ipos;
    pos.set( mouse_begin_x, mouse_begin_y );
    Draw2D::drawPointCross_d( pos, r );
    Draw2D::drawCircle_d( pos, r, 64, false );
    world.colruler.pos2index( pos, dipos, ipos );
    if( (ipos.x<1) || (ipos.x>=world.colbuf.NX-1) || (ipos.y<1) || (ipos.y>=world.colbuf.NY-1) ) return;
    world.colbuf.insert( NULL, ipos, dipos, rf );
    for( int iy=0; iy<world.colbuf.NY; iy++ ){
        double y = world.colruler.i2y( iy );
        for( int ix=0; ix<world.colbuf.NX; ix++ ){
            double x = world.colruler.i2x( ix );
            int ixy =  world.colbuf.xy2i( ix, iy );
            int ncount = world.colbuf.counts[ixy];
            for( int im=0; im<ncount; im++ ){
                Draw2D::drawRectangle_d( {x,y}, {x+world.colruler.step.x,y+world.colruler.step.y}, false );
            }
        }
    }
}

FormationTacticsApp::FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    world.init();

    currentFaction   = world.factions[0];  printf( "currentFaction: %s\n", currentFaction->name );
    currentFormation = currentFaction->formations[0];

    TiledView::init( 6, 6 );
    tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );
}

void FormationTacticsApp::draw(){
    long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );



    long tComp = getCPUticks();
    world.update( );
    tComp = getCPUticks() - tComp;


    // DEBUG grid insert
    debug_buffinsert( );

    long tDraw = getCPUticks();
    for( Formation* fm : world.formations ){
        if( fm->bbox.notOverlaps( {camXmin,camYmin, camXmax,camYmax} ) ) continue;
        //printf( " f %i \n", f  );
        //glColor3f( fm->faction->color.x, fm->faction->color.y, fm->faction->color.z );
        if( (fm!= NULL) ){
            if  ( fm == currentFormation ){ fm->render( fm->faction->color, formation_view_mode );   }
            else                          { fm->render( fm->faction->color, 0                   );   }
        }
    }
    tDraw = getCPUticks() - tDraw;

    if( currentFormation != 0 ){
        //glColor3f(1.0,0.0,1.0);
        glColor3f(0.0,1.0,0.0);
        Draw2D::drawCircle_d( currentFormation->center, 0.5, 16, true );
    }

    /*
    for( Faction* fa : world.factions ){
        glColor3f( fa->color.x, fa->color.y, fa->color.z );
        for( Formation* fm : fa->formations ){
            //printf( " f %i \n", f  );
            if( fm != NULL ) fm->render( );
        }
    }
    */

    //if(frameCount > 10) STOP = true;

    //printf( " frame %i soldiers %i interactions %i time %i  \n" );

    printf( " frame %i : %i %i   %4.3f %4.3f %4.3f\n", frameCount, world.nSoldiers, world.nSoldierInteractions,
                                 tComp*1.0e-6, tComp/(double)world.nSoldiers, tComp/(double)world.nSoldierInteractions );

    tTot = getCPUticks() - tTot;
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
                case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    if( currentFaction != NULL ) currentFormation = currentFaction->getFormationAt( { mouse_begin_x, mouse_begin_y } );
                break;
                case SDL_BUTTON_RIGHT:
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
    camStep = zoom*0.05;
}

// ===================== MAIN

FormationTacticsApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
    SDL_DisplayMode dm;
    SDL_GetDesktopDisplayMode(0, &dm);
	thisApp = new FormationTacticsApp( junk , dm.w-150, dm.h-100 );
	SDL_SetWindowPosition(thisApp->window, 100, 0 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















