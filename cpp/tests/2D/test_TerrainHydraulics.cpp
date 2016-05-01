
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

// ======================  TestApp

class TestAppTerrainHydraulics : public AppSDL2OGL{
	public:
    int perframe = 20;
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
        glColor3f(g*c,0.7*(1-g)*c,1.0);
    }else{
        glColor3f(g,0.7*(1-g),0);
    }
    return g;
}

void TestAppTerrainHydraulics::renderMapContent( float x0, float y0, float scale, float csc, float hsc ){
    //glColor3f( 0.1f,0.1f,0.1f );
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); b.mul(scale);
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

    terrain.allocate( 512, 512 );
    //terrain.genTerrainNoise( 14, 0.3, 0.7, 0.6, 45454, {1000.0,1000.0} );
    //terrain.genTerrainNoise( 14, 0.5, 1.0,  0.7, 1.2, 45454, {100.0,100.0} );
    terrain.genTerrainNoise( 8, 2.0, 1.0,  0.5, 0.8, 45454, {100.0,100.0} );

/*
    shape=glGenLists(1);
	glNewList( shape, GL_COMPILE );
        renderMapContent( -0.1*terrain.nx, -0.1*terrain.ny,  0.2, 1.0, 1.0 );
	glEndList();
*/

    terrain.init_outflow();
    for ( int ii=0; ii<50; ii++){
        int i = rand()%terrain.ntot;
        terrain.contour2[ii] = i;
        terrain.nContour++;
        terrain.water[i] = terrain.ground[i];

    }
    //for (int i=0; i<terrain.ntot; i++){  terrain.water[i] = terrain.ground[i]; }

    shape=glGenLists(1);
}

void TestAppTerrainHydraulics::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    if( running ){
        for( int i=0; i<perframe; i++ ){
            terrain.outflow_step();
            if( terrain.nContour == 0 ){  running=false; break; }
        }
        if( running ){
            renderMapContent( -0.1*terrain.nx, -0.1*terrain.ny,  0.2, 1.0, 1.0 );
        }else{
            glNewList( shape, GL_COMPILE );
            renderMapContent( -0.1*terrain.nx, -0.1*terrain.ny,  0.2, 1.0, 1.0 );
            glEndList();
        }
    }else{
        glCallList( shape );
    }

    //renderMapContent( );

};

void TestAppTerrainHydraulics::drawHUD(){}


void TestAppTerrainHydraulics::eventHandling( const SDL_Event& event ){
    /*
    switch( event.type ){
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
    };
    */
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
















