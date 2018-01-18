
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

#include "TerrainHydraulics.h"
#include "testUtils.h"


#include "GridIndex2D.h"
#include "Grid2DAlgs.h"
#include "Grid2DAlgs.cpp" // FIXME
#include "SquareRuler.h"

//#include "Grid.h"

// ======================  TestApp

class TestAppTerrainHydraulics : public AppSDL2OGL{
	public:
    int job_type = 1;
    int perframe = 3;
    TerrainHydraulics terrain;
    int shape;
    bool running = true;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

    void renderMapContent ( float x0, float y0, float scale, float csc, float hsc );
    double terrain_color( int i );
	//void drawSimplexGrid( int n, float step );

	TestAppTerrainHydraulics( int& id, int WIDTH_, int HEIGHT_ );

};

double TestAppTerrainHydraulics::terrain_color( int i ){
    float g   = terrain.ground[i];
    float w   = terrain.water [i];
    if( w > g ){
        float c = (1-20*(w-g)); if(c<0) c=0;
        glColor3f(g*g*c,0.2+0.8*g*(1-g)*c,0.5);
    }else{
        glColor3f(g*g,0.2+0.8*(1-g)*g,0);
    }
    return g;
}

void TestAppTerrainHydraulics::renderMapContent( float x0, float y0, float scale, float csc, float hsc ){
    //glColor3f( 0.1f,0.1f,0.1f );
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); b.mul(scale);
    //b.set( 0.0d, 1.0 ); b.mul(scale);
    //glDisable(GL_SMOOTH);
    int ii = 0;
    for (int iy=0; iy<terrain.ny-1; iy+=1){
        glBegin( GL_TRIANGLE_STRIP );
        for (int ix=0; ix<terrain.nx; ix+=1){
            p.set( ix*a.x+iy*b.x + x0, ix*a.y+iy*b.y + y0 );
            //val = terrain.ground[ii              ]; w = terrain.water[ii              ];    c = csc * val;  glColor3f(c,c,w);  glVertex3f( p.x    , p.y    , val*hsc );
            //val = terrain.ground[ii + terrain.nx ]; w = terrain.water[ii + terrain.nx ];    c = csc * val;  glColor3f(c,c,w);  glVertex3f( p.x+b.x, p.y+b.y, val*hsc );

            terrain_color( ii               ); glVertex3f( p.x    , p.y    , 0 );
            terrain_color( ii + terrain.nx  ); glVertex3f( p.x+b.x, p.y+b.y, 0 );

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

TestAppTerrainHydraulics::TestAppTerrainHydraulics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    srand(548);
    terrain.allocate( 512, 512 ); bisecNoise( 9, terrain.ground, -1.0/512, 1.0/512 );
    //terrain.allocate( 16, 16 ); bisecNoise( 4, terrain.ground, -0.5/16, 0.5/16 );
    //terrain.genTerrainNoise( 14, 0.3, 0.7, 0.6, 45454, {1000.0,1000.0} );
    //terrain.genTerrainNoise( 14, 0.5, 1.0,  0.7, 1.2, 45454, {100.0,100.0} );
    //terrain.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );


/*
    shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        renderMapContent( -0.1*terrain.nx, -0.1*terrain.ny,  0.2, 1.0, 1.0 );
	glEndList();
*/

    //for (int i=0; i<terrain.ntot; i++){  terrain.water[i] = terrain.ground[i]; }

    //for (int i=0; i<terrain.ntot; i++){  terrain.water[i] = 1.0; terrain.water_[i] = 1.0; }
    //terrain.initErrosion( 1.0 );

    shape=glGenLists(1);
}

void TestAppTerrainHydraulics::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );
    long t0;
    if( running ){
        t0 = getCPUticks();
        //perframe = 1;
        for( int i=0; i<perframe; i++ ){

            //terrain.errodeDroples( 10000, 100, 0.00, 0.2, 0.5 );
            for( int j=0; j<20; j++ ){
                int isz = 25;
                int ix0 = rand()%(terrain.nx-isz);
                int iy0 = rand()%(terrain.ny-isz);
                terrain.errodeDroples( 200, 100, 0.02, 0.15, 0.5, ix0, iy0, ix0+isz, iy0+isz );
            }
            /*
            switch(job_type){
                case 1: {
                    int npix = terrain.flow_errosion_step_noRain( );
                    if( npix < terrain.nx ){
                        job_type = 2;
                        terrain.init_outflow();
                        for ( int ii=0; ii<50; ii++){
                            int i = rand()%terrain.ntot;
                            terrain.contour2[ii] = i;
                            terrain.nContour++;
                            terrain.water[i] = terrain.ground[i];
                        }
                        //running = false;
                    }
                    //terrain.errodeDroples( 1000, 100, 0.01 );
                    break; }
                case 2:{
                    terrain.outflow_step();
                    if( terrain.nContour == 0 ){
                        //job_type = 2;
                        running=false; break;
                    }
                    break;}
            }
            */
        }
        long tcomp = getCPUticks() - t0;
        t0    = getCPUticks();
        if( running ){
            renderMapContent( -0.1*terrain.nx, -0.1*terrain.ny,  0.2, 1.0, 1.0 );
        }else{
            glDeleteLists(shape,1);
            glNewList( shape, GL_COMPILE );
            renderMapContent( -0.1*terrain.nx, -0.1*terrain.ny,  0.2, 1.0, 1.0 );
            glEndList();
        }
        long tplot = getCPUticks() - t0;
        //printf( " tplot %3.3f Mtick tcomp %3.3f Mtick ( %3.3f ticks/pix ) \n", tplot*1e-6, tcomp*1e-6, ((double)tcomp)/(terrain.ntot*perframe) );
    }else{
        glCallList( shape );
    }


    //renderMapContent( );

};

void TestAppTerrainHydraulics::drawHUD(){}


void TestAppTerrainHydraulics::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_r:  terrain.initErrosion( 0.8 ); running=true;  break;
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

TestAppTerrainHydraulics * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppTerrainHydraulics( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















