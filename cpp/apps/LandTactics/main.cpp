
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_image.h>
//#include <SDL2/SDL_ttf.h>
//#include "Texture.h"

#include "Draw2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"
#include "SDL_utils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "LTUnitType.h"
#include "LTUnit.h"
#include "LTFaction.h"
#include "LTWorld.h"

// font rendering:
//  http://www.willusher.io/sdl2%20tutorials/2013/12/18/lesson-6-true-type-fonts-with-sdl_ttf
//  http://stackoverflow.com/questions/28880562/rendering-text-with-sdl2-and-opengl

int   default_font_texture;


void cmapHeight(double g){

    //glColor3f(0.2+0.8*g*g,0.2+0.3*(1-g)*g,0.2);
    //double snow = g*g*g*g;
    //glColor3f(g*0.8+0.2,0.5+0.5*snow,0.2+0.8*snow);

    glColor3f(g,g,g);

    //return g;
}

class FormationTacticsApp : public AppSDL2OGL {
	public:
    LTWorld world;

    int formation_view_mode = 0;
    LTUnit    * currentUnit      = NULL;
    LTFaction * currentFaction   = NULL;
    int       ifaction = 0;

    bool bDrawing = false;

    //GLuint       itex;

    // ==== function declaration
    void printASCItable( int imin, int imax  );
    //GLuint makeTexture( char * fname );
    //GLuint renderImage( GLuint itex, const Rect2d& rec );
    //void drawString( char * str, int imin, int imax, float x, float y, float sz, int itex );
    //void drawString( char * str, float x, float y, float sz, int itex );

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling( const SDL_Event& event  );
	void debug_buffinsert( );
	//void pickParticle( Particle2D*& picked );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

	FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ );

};

void FormationTacticsApp::printASCItable( int imin, int imax  ){
    int len = imax-imin;
    char str[len];
    for ( int i=0; i<len; i++ ){
        str[i] = (char)(i+imin);
    }
    printf("%s\n", str );
};

FormationTacticsApp::FormationTacticsApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    /*
    for(int i=0; i<100; i++){
        //Vec3d hs = {0.0,1.0,1.0};
        Vec3d hs = {randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};
        Vec2d p  = {0.2, 0.2};
        double h   = trinagleInterp( toBaricentric(p), hs );
        double hdx = trinagleInterp( toBaricentric({p.x+0.1,p.y }), hs );
        double hdy = trinagleInterp( toBaricentric({p.x,p.y+0.1 }), hs );
        Vec2d  dh  = trinagleDeriv ( hs );
        printf( " %i (%g,%g) (%g,%g) \n", i, (hdx-h)/0.1, (hdy-h)/0.1, dh.x, dh.y );
    }
    exit(0);
    */

    printASCItable( 33, 127  );

    world.init();

    camX0 = world.map_center.x;
    camY0 = world.map_center.y;

    currentFaction = world.factions[0];  printf( "currentFaction: %s\n", currentFaction->name );
    currentUnit    = currentFaction->units[0];

    //TiledView::init( 6, 6 );
    //tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );

    default_font_texture = makeTexture(  "common_resources/dejvu_sans_mono.bmp" );
    //default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //itex = makeTexture(  "data/tank.bmp" );
    //itex = makeTexture(  "data/nehe.bmp" );
    printf( "default_font_texture :  %i \n", default_font_texture );

}

void FormationTacticsApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	if(bDrawing){
        Vec2i ind;Vec2d dind;
        //int i = world.ruler.simplexIndex({mouse_begin_x,mouse_begin_y},ind,dind)>>1;
        world.ruler.simplexIndexBare({mouse_begin_x,mouse_begin_y},ind);
        ind = world.ruler.wrap_index(ind);
        int i = world.ruler.ip2i(ind);
        printf("i %i (%i,%i) \n",i, ind.x, ind.y );
        world.ground[i] = randf(0.0,1.0);
	}

	glPushMatrix();
	glScalef(world.ruler.step,world.ruler.step,1.0);
	Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, world.ruler.na, world.ground, 0.0, world.maxHeight );
	glPopMatrix();


    //float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    //TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );
    //long tComp = getCPUticks();
    world.update( );
    //tComp = getCPUticks() - tComp;

    long tDraw = getCPUticks();
    for( LTUnit* u : world.units ){
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

    double h = world.ruler.getValue( {mouse_begin_x,mouse_begin_y}, world.ground );
	cmapHeight( h/world.maxHeight );
    //Draw2D::drawPointCross_d( {mouse_begin_x,mouse_begin_y}, 100 );
    Draw2D::drawCircle_d( {mouse_begin_x,mouse_begin_y}, 0.35, 8, true );
    Vec2d  dh  = world.ruler.getDeriv( {mouse_begin_x    ,mouse_begin_y }, world.ground );
    glColor3f(0.0f,1.0f,0.0f); Draw2D::drawVecInPos_d( dh*-50.0, {mouse_begin_x,mouse_begin_y} );
    //double hdx = world.ruler.getValue( {mouse_begin_x+0.1,mouse_begin_y    }, world.ground );
    //double hdy = world.ruler.getValue( {mouse_begin_x    ,mouse_begin_y+0.1}, world.ground );
    //printf( "(%g,%g) (%g,%g) \n", (hdx-h)/0.1, (hdy-h)/0.1, dh.x, dh.y );

};

/*
int FormationTacticsApp::tileToList( float x0, float y0, float x1, float y1 ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
		world.terrain.renderRect( x0, y0, x1, y1, 31 );
		//glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
	glEndList();
	return ilist;
}
*/

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
                    //bDrawing=true;
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
















