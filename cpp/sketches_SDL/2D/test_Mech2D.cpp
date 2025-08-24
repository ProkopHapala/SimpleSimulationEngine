/// @file @brief This demo simulates a 2D mechanical system using `Mech2D.h`, which likely involves particles and forces like springs or gravity. It visualizes the dynamic behavior of the system, allowing the user to pause/resume the simulation with the `SPACE` key and interactively apply forces by dragging particles with the mouse.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "globals.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "MechGrid2D.h"
#include "MechMesh2D.h"


bool bDEBUG = true;

#include "MechPIC2D.h"
#include "MechPIC2D_Temperature.h"


/*

ToDo : MechPIC2D_T needs better than linear interpolation in order to make smooth forces


*/

CompressibleMaterial materials[] = {
    {1.,1.}
};




class TestAppMech2D : public AppSDL2OGL { public:

    //MechGrid2D mgrid;
    //MechMesh2D mmesh;
    //MechPIC2D    mpic;
    MechPIC2D_T  mpic;

    int nproj=0;
    //bool bRun = false;
    bool bRun = true;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMech2D::TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {


    //mpic.realloc( 8, {16,16} );

    /*
    int nside = 5;
    nproj = 0;
    mpic.realloc( (nside-1)*(nside-1) + nproj, {nside,nside} );

    int nrow = nside-1;
    // --- Target
    for(int i=0; i<mpic.np-nproj; i++){
        int ix = i%nrow;
        int iy = i/nrow;
        mpic.imats [i] = 0;          // materials
        mpic.pmoles[i] = 0.1;       //
        //mpic.pos   [i] = { (ix+16.2)*mpic.step,(iy+24.3)*mpic.step };
        mpic.pos   [i] = { (ix+1+randf(-0.5,0.5))*mpic.step,(iy+1+randf(-0.5,0.5))*mpic.step };
        //mpic.vel[i]  = {0.,1000.0};
        mpic.vel   [i] = Vec2dZero;  // velocity
        //printf( "target[%i] p(%g,%g) v(%g,%g) \n", i, mpic.pos[i].x, mpic.pos[i].y,  mpic.vel[i].x, mpic.vel[i].y );
    }
    // --- Projectile
    int i0=mpic.np-nproj;
    for(int i=i0; i<mpic.np; i++){
        mpic.imats [i] = 0;
        mpic.pmoles[i] = 0.1;
        mpic.pos[i]    = { ((i-i0)+16.2)*mpic.step, 10.5*mpic.step };
        mpic.vel[i]  = Vec2dY * 10.0;  // velocity
        //mpic.vel[i]    = Vec2dZero;  // velocity
        //printf( "projectile[%i] p(%g,%g) v(%g,%g) \n", i, mpic.pos[i].x, mpic.pos[i].y,  mpic.vel[i].x, mpic.vel[i].y );
    }
    */

    mpic.realloc( 9, {5,5} );
    mpic.setStep( 1e-2 );
    int nrow=3;
    for(int i=0; i<mpic.np; i++){
        int ix = i%nrow;
        int iy = i/nrow;
        mpic.imats [i] = 0;
        mpic.pmoles[i] = 0.1;
        mpic.pos   [i] = { (ix+1)*mpic.step,(iy+1)*mpic.step };
        mpic.vel   [i] = Vec2dZero;
    }

    mpic.materials = materials;


    camX0 = mpic.nc.x/2;
    camY0 = mpic.nc.y/2;

    mpic.particlesToCells();
    mpic.updateCellTN ( 1.0 );
    mpic.moveMD( 0.0001 );
    bDEBUG=false;
}

void TestAppMech2D::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

    glColor3f(1.0,0.0,0.0);

    int ixy=0;
    for(int iy=0; iy<mpic.nc.y; iy++){
        for(int ix=0; ix<mpic.nc.x; ix++){
            double c = mpic.moles[ixy]/0.5;
            //double c = mpic.temperature[ixy]*0.0001;
            //double c = mpic.Umol[ixy];
            //glColor3b( 1+c,1-c*c, 1-c );
            glColor3f( 1-c, 1-c*4., 1-c*16. );
            //Vec2d vert = {ix*mpic.step,iy*mpic.step};
            //Draw2D::drawRectangle( (ix-0.5)*mpic.step, (iy-0.5)*mpic.step, (ix+0.5)*mpic.step, (iy+0.5)*mpic.step, true );
            Draw2D::drawRectangle( (ix-0.5), (iy-0.5), (ix+0.5), (iy+0.5), true );
            ixy++;
        }
    }


    double dt = 1e-7;
    //double dt_debug = 0.0001;
    double dt_debug = 0.0001;
    if(bRun){

        Vec2d u; u.fromAngle(frameCount*0.1); u.mul(0.00001);
        mpic.pos[4].add( u );
        mpic.pmoles[4]=0;

        mpic.particlesToCells();
        //mpic.updateCellThermodynamics();
        //mpic.setCellT( 10.0 );
        mpic.updateCellTN ( 1.0 );
        //mpic.moveMD( dt );
        //mpic.moveMD_selfCorrect( dt );
        //mpic.moveMD( dt );

        mpic.moveMD( dt_debug );
        //mpic.boundParticleMirror();
    }


    //if(frameCount<100){ printf( "========= %i \n", frameCount ); }

    Draw2D::drawRectangle( (Vec2f)(mpic.pmin*mpic.invStep), (Vec2f)(mpic.pmax*mpic.invStep), false );

    glColor3f(0.,0.,0.);
    //double vsc   = 1e-3;
    double vsc   = 0.1;
    double molsc = 5.1;
    for(int i=0; i<mpic.np; i++){
        //Draw2D::drawPointCross_d( mpic.pos[i], mpic.pmoles[i]*molsc );
        //Draw2D::drawVecInPos_d  ( mpic.vel[i]*vsc, mpic.pos[i] );
        //Draw2D::drawPointCross_d( mpic.pos[i]*mpic.invStep, mpic.pmoles[i]*molsc );
        Vec2d p = mpic.pos[i]*mpic.invStep;
        //if(frameCount<100){ printf("p[%i] p(%g,%g) \n", i, p.x, p.y ); };
        Draw2D::drawPointCross_d( p, 0.1 );
        Draw2D::drawVecInPos_d  ( mpic.vel[i]*vsc, mpic.pos[i]*mpic.invStep );
        //if( i>=(mpic.np-nproj) ){ printf( "projectile[%i] p(%g,%g) v(%g,%g) \n", i, mpic.pos[i].x, mpic.pos[i].y,  mpic.vel[i].x, mpic.vel[i].y ); }
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
