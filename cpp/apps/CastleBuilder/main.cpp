
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

//#include "TerrainSimplex.h"
#include "CastleWorld.h"

#include "testUtils.h"

// ======================  TestApp

class CastleBuilderSingle : public AppSDL2OGL{
	public:
    int job_type = 1;
    int perframe = 3;
    CastleWorld world;
    int shape;
    bool running = true;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

    void renderMapContent ( float x0, float y0, float scale, float csc, float hsc );
    double world_color( int i );
	//void drawSimplexGrid( int n, float step );

	CastleBuilderSingle( int& id, int WIDTH_, int HEIGHT_ );

};

double CastleBuilderSingle::world_color( int i ){
    float g   = world.ground[i];
    float w   = world.water [i];
    if( w > g ){
        float c = (1-20*(w-g)); if(c<0) c=0;
        glColor3f(g*g*c,0.2+0.8*g*(1-g)*c,0.5);
    }else{
        glColor3f(g*g,0.2+0.8*(1-g)*g,0);
    }
    return g;
}

void CastleBuilderSingle::renderMapContent( float x0, float y0, float scale, float csc, float hsc ){
    //glColor3f( 0.1f,0.1f,0.1f );
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); b.mul(scale);
    //glDisable(GL_SMOOTH);
    int ii = 0;
    for (int iy=0; iy<world.ny-1; iy+=1){
        glBegin( GL_TRIANGLE_STRIP );
        for (int ix=0; ix<world.nx; ix+=1){
            p.set( ix*a.x+iy*b.x + x0, ix*a.y+iy*b.y + y0 );
            //val = world.ground[ii              ]; w = world.water[ii              ];    c = csc * val;  glColor3f(c,c,w);  glVertex3f( p.x    , p.y    , val*hsc );
            //val = world.ground[ii + world.nx ]; w = world.water[ii + world.nx ];    c = csc * val;  glColor3f(c,c,w);  glVertex3f( p.x+b.x, p.y+b.y, val*hsc );

            world_color( ii               ); glVertex3f( p.x    , p.y    , 0 );
            world_color( ii + world.nx  ); glVertex3f( p.x+b.x, p.y+b.y, 0 );

            //val = sin(p.x    )*sin(p.y    ); c = cscale * val;  glColor3f(c,c,c);  glVertex3f( p.x    , p.y    , val*hscale );
            //val = sin(p.x+b.x)*sin(p.y+b.y); c = cscale * val;  glColor3f(c,c,c);  glVertex3f( p.x+b.x, p.y+b.y, val*hscale );

            //glColor3f(randf(),randf(),randf()); glVertex3f( p.x    , p.y    , 0 );
            //glColor3f(randf(),randf(),randf()); glVertex3f( p.x+b.x, p.y+b.y, 0 );
            //printf( "%i, %i %f \n", iy, ix, val );
            ii++;
        }
        glEnd();
    }
}

CastleBuilderSingle::CastleBuilderSingle( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    world.allocate( 512, 512 );
    //world.genTerrainNoise( 14, 0.3, 0.7, 0.6, 45454, {1000.0,1000.0} );
    //world.genTerrainNoise( 14, 0.5, 1.0,  0.7, 1.2, 45454, {100.0,100.0} );
    world.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );

/*
    shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        renderMapContent( -0.1*world.nx, -0.1*world.ny,  0.2, 1.0, 1.0 );
	glEndList();
*/

    //for (int i=0; i<world.ntot; i++){  world.water[i] = world.ground[i]; }

    //for (int i=0; i<world.ntot; i++){  world.water[i] = 1.0; world.water_[i] = 1.0; }
    world.initErrosion( 1.0 );

    shape=glGenLists(1);
}

void CastleBuilderSingle::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );
    long t0;
    if( running ){
        t0 = getCPUticks();
        //perframe = 1;
        for( int i=0; i<perframe; i++ ){
            switch(job_type){
                case 1: {
                    int npix = world.flow_errosion_step_noRain( );
                    if( npix < world.nx ){
                        job_type = 2;
                        world.init_outflow();
                        for ( int ii=0; ii<50; ii++){
                            int i = rand()%world.ntot;
                            world.contour2[ii] = i;
                            world.nContour++;
                            world.water[i] = world.ground[i];
                        }
                        //running = false;
                    }
                    //world.errodeDroples( 1000, 100, 0.01 );
                    break; }
                case 2:{
                    world.outflow_step();
                    if( world.nContour == 0 ){
                        //job_type = 2;
                        running=false; break;
                    }
                    break;}
            }
        }
        long tcomp = getCPUticks() - t0;
        t0    = getCPUticks();
        if( running ){
            renderMapContent( -0.1*world.nx, -0.1*world.ny,  0.2, 1.0, 1.0 );
        }else{
            glNewList( shape, GL_COMPILE );
            renderMapContent( -0.1*world.nx, -0.1*world.ny,  0.2, 1.0, 1.0 );
            glEndList();
        }
        long tplot = getCPUticks() - t0;
        printf( " tplot %3.3f Mtick tcomp %3.3f Mtick ( %3.3f ticks/pix ) \n", tplot*1e-6, tcomp*1e-6, ((double)tcomp)/(world.ntot*perframe) );
    }else{
        glCallList( shape );
    }


    //renderMapContent( );

};

void CastleBuilderSingle::drawHUD(){}


void CastleBuilderSingle::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_r:  world.initErrosion( 0.8 ); running=true;  break;
            }
            break;
         /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //paintSimplex( mouse_begin_x, mouse_begin_y );
                    mouse_left  = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = true;
                    //eraseSimplex( mouse_begin_x, mouse_begin_y );
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    mouse_left = false;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = false;
                    break;
            }
            break;
        */
    };
    AppSDL2OGL::eventHandling( event );
};

// ===================== MAIN

CastleBuilderSingle * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new CastleBuilderSingle( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















