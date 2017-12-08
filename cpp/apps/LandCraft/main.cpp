
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
#include "Draw.h"
#include "Draw2D.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "AppSDL2OGL.h"

#include "testUtils.h"
#include "SDL_utils.h"

#include "TerrainCubic.h"
#include "TiledView.h"

#include "SimplexRuler.h"
#include "Ruler2DFast.h"
#include "TerrainHydraulics.h"

#include "IO_utils.h"


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

class LandCraftApp : public AppSDL2OGL {
	public:

    SimplexRuler       ruler;
    Ruler2DFast        square_ruler;
    TerrainHydraulics  hydraulics;
    double * ground    = NULL;
    double * water    = NULL;

    Vec2d map_center;
    double maxHeight = 500.0;

    //int       ifaction = 0;
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

	LandCraftApp( int& id, int WIDTH_, int HEIGHT_ );

};

void LandCraftApp::printASCItable( int imin, int imax  ){
    int len = imax-imin;
    char str[len];
    for ( int i=0; i<len; i++ ){
        str[i] = (char)(i+imin);
    }
    printf("%s\n", str );
};

LandCraftApp::LandCraftApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

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


    ruler.setSize(128,128);
    ruler.setStep(50);

    map_center = (Vec2d){ruler.na*0.75*ruler.step,ruler.nb*0.5*ruler.step};

    ground = new double[ruler.ntot];
    hydraulics.setSize(ruler.na,ruler.nb);
    hydraulics.ground = ground;

    bool newMap = false;
    if( newMap ){
        hydraulics.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );
        for( int j=0; j<500; j++ ){
            int isz = 25;
            int ix0 = rand()%(hydraulics.nx-isz);
            int iy0 = rand()%(hydraulics.ny-isz);
            hydraulics.errodeDroples( 200, 100, 0.02, 0.15, 0.5, ix0, iy0, ix0+isz, iy0+isz );
        }
        for(int i=0; i<ruler.ntot; i++){ ground[i] *= maxHeight; };
        //for(int i=0; i<ruler.ntot; i++){ ground[i] = randf(0.0,500.0); };
        saveBin( "data/ground.bin", sizeof(double)*hydraulics.ntot, (char*)ground );
    }else{
        loadBin( "data/ground.bin", sizeof(double)*hydraulics.ntot, (char*)ground );
    }

    hydraulics.allocate_outflow();
    water = new double[ruler.ntot];
    hydraulics.water = water;
    hydraulics.init_outflow( 500.0 );

    hydraulics.nContour = 1;
    int idrain = hydraulics.nx/2 + hydraulics.nx*hydraulics.ny/2;
    printf( "idrain  %i \n", idrain   );
    hydraulics.contour2[0] = idrain;
    water[idrain]          = ground[idrain];

    //TiledView::init( 6, 6 );
    //tiles    = new int[ nxy ];
    //TiledView::renderAll( -10, -10, 10, 10 );

    default_font_texture = makeTexture(  "common_resources/dejvu_sans_mono.bmp" );
    //default_font_texture = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //itex = makeTexture(  "data/tank.bmp" );
    //itex = makeTexture(  "data/nehe.bmp" );
    printf( "default_font_texture :  %i \n", default_font_texture );

}

void LandCraftApp::draw(){
    //long tTot = getCPUticks();
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	if(bDrawing){
        Vec2i ind;Vec2d dind;
        //int i = ruler.simplexIndex({mouse_begin_x,mouse_begin_y},ind,dind)>>1;
        ruler.simplexIndexBare({mouse_begin_x,mouse_begin_y},ind);
        ind = ruler.wrap_index(ind);
        int i = ruler.ip2i(ind);
        //printf("i %i (%i,%i) \n",i, ind.x, ind.y );
        ground[i] = randf(0.0,1.0);
	}

	hydraulics.outflow_step();
	//printf( " nContour %i \n", hydraulics.nContour );
	//if(frameCount>150)exit(0);

	glPushMatrix();
	glScalef(ruler.step,ruler.step,1.0);
	//Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, ruler.na, ground, 0.0, maxHeight );
	Draw2D::drawTriaglePatch<cmapHeight>( {0,0}, {128,128}, ruler.na, water, 0.0, maxHeight );
	glPopMatrix();





    //float camMargin = ( camXmax - camXmin )*0.1;
    //float camMargin = 0;
    //TiledView::draw(  camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin  );
    //printf( " camRect  %f %f %f %f \n", camXmin-camMargin, camYmin-camMargin, camXmax+camMargin, camYmax+camMargin );
    //long tComp = getCPUticks();
    //update( );
    //tComp = getCPUticks() - tComp;



    printf("===== frameCount %i \n", frameCount);



};

/*
int LandCraftApp::tileToList( float x0, float y0, float x1, float y1 ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
		terrain.renderRect( x0, y0, x1, y1, 31 );
		//glColor3f(0.9f,0.2f,0.2f); Draw2D::drawRectangle( x0+0.1, y0+0.1, x1-0.1, y1-0.1, false );
	glEndList();
	return ilist;
}
*/

void LandCraftApp::drawHUD(){}

void LandCraftApp::eventHandling ( const SDL_Event& event  ){
    //printf( "NBodyWorldApp::eventHandling() \n" );
    int ihex;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_0:  formation_view_mode = 0;            printf( "view : default\n" ); break;
                //case SDLK_1:  formation_view_mode = VIEW_INJURY;  printf( "view : injury\n"  ); break;
                //case SDLK_2:  formation_view_mode = VIEW_STAMINA; printf( "view : stamina\n" ); break;
                //case SDLK_3:  formation_view_mode = VIEW_CHARGE;  printf( "view : charge\n"  ); break;
                //case SDLK_4:  formation_view_mode = VIEW_MORAL;   printf( "view : moral\n"   ); break;
                //case   SDLK_n: ifaction++; if(ifaction>=factions.size()) ifaction=0; currentFaction = factions[ifaction]; printf("ifaction %i\n",ifaction); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ihex = ruler.hexIndex( {mouse_begin_x, mouse_begin_y} );
                    hydraulics.water[ihex] = 1000.0;
                    //printf( "left button pressed !!!! " );
                    //if( currentFaction != NULL ) currentUnit = currentFaction->getUnitAt( { mouse_begin_x, mouse_begin_y } );
                    //bDrawing=true;
                    break;
                case SDL_BUTTON_RIGHT:
                    printf( "left button pressed !!!! " );
                    break;
            }
            break;
            /*
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    picked = NULL;
                    break;
            }
            break;
            */
    };
    AppSDL2OGL::eventHandling( event );
    camStep = zoom*0.05;
}

// ===================== MAIN

LandCraftApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new LandCraftApp( junk , 800, 600 );
	thisApp->zoom  = 3000;
	thisApp->camX0 = 4000;
	thisApp->camY0 = 2000;
	thisApp->loop( 1000000 );
	return 0;
}
















