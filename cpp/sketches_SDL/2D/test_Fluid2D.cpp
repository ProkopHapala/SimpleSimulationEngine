
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "MechGrid2D.h"
#include "MechMesh2D.h"

#include "MechPIC2D.h"
#include "MechPIC2D_Temperature.h"

CompressibleMaterial materials[] = {
    {1.,1.}
};



class TestAppMech2D : public AppSDL2OGL { public:

    //MechGrid2D mgrid;
    //MechMesh2D mmesh;
    //MechPIC2D    mpic;
    MechPIC2D_T  mpic;

    //bool bRun = false;
    bool bRun = true;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMech2D::TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    mpic.realloc( 8, {16,16} );
    mpic.materials = materials;

    mpic.setStep( 1e-2 );

    for(int i=0; i<mpic.np; i++){
        mpic.imats [i] = 0;          // materials
        mpic.pmoles[i] = 0.1;       //
        mpic.pos[i]    = { (i+2.2)*mpic.step,4.3*mpic.step };   // position (global, or local within the box)
        //mpic.vel[i]    = {0.,1000.0};
        mpic.vel[i]    = Vec2dZero;  // velocity
    }

}

void TestAppMech2D::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glDisable(GL_LIGHTING);

    glColor3f(1.0,0.0,0.0);

    double dt = 1e-7;
    if(bRun){

        mpic.particlesToCells();
        mpic.updateCellThermodynamics();
        //mpic.moveMD( dt );
        //mpic.moveMD_selfCorrect( dt );
        mpic.moveMD( dt );
        //mpic.update(dt);
    }

    int ixy=0;
    for(int iy=0; iy<mpic.nc.y; iy++){
        for(int ix=0; ix<mpic.nc.x; ix++){
            //double c = mpic.moles[ixy]/0.1;
            double c = mpic.temperature[ixy]*0.0001;
            //double c = mpic.Umol[ixy];
            //glColor3b( 1+c,1-c*c, 1-c );
            glColor3f( 1-c, 1-c*4., 1-c*16. );
            //Vec2d vert = {ix*mpic.step,iy*mpic.step};
            //Draw2D::drawRectangle( (ix-0.5)*mpic.step, (iy-0.5)*mpic.step, (ix+0.5)*mpic.step, (iy+0.5)*mpic.step, true );
            Draw2D::drawRectangle( (ix-0.5), (iy-0.5), (ix+0.5), (iy+0.5), true );
            ixy++;
        }
    }

    glColor3f(0.,0.,0.);
    double vsc   = 1e-3;
    double molsc = 5.1;
    for(int i=0; i<mpic.np; i++){
        //Draw2D::drawPointCross_d( mpic.pos[i], mpic.pmoles[i]*molsc );
        //Draw2D::drawVecInPos_d  ( mpic.vel[i]*vsc, mpic.pos[i] );

        Draw2D::drawPointCross_d( mpic.pos[i]*mpic.invStep, mpic.pmoles[i]*molsc );
        Draw2D::drawVecInPos_d  ( mpic.vel[i]*vsc, mpic.pos[i]*mpic.invStep );
    }

};

void TestAppMech2D::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
}

void TestAppMech2D::eventHandling ( const SDL_Event& event  ){
    //printf( "TestAppMech2D::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //pickParticle( world.picked );
                    //ipick = pline1.nearestPoint( {mouse_begin_x,mouse_begin_y} );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //world.picked = NULL;
                    //ipick = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppMech2D::drawHUD(){

}

// ===================== MAIN

TestAppMech2D * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new TestAppMech2D( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
