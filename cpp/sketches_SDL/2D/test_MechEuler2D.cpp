
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

#include "MechEuler2D.h"

class TestAppMech2D : public AppSDL2OGL { public:

    //MechGrid2D mgrid;
    //MechMesh2D mmesh;
    //MechPIC2D    mpic;
    MechEuler2D  solver;

    //bool bRun = false;
    bool bRun = true;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMech2D::TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {


    /*
    // Mini test
    solver.realloc({16,16}, 2 );
    solver.clearVelocitis();
    solver.clearMater    ();
    solver.matNs[0][solver.nc.x*8 + 7] = 0.3;
    solver.matNs[0][solver.nc.x*8 + 8] = 5.0;
    solver.matNs[0][solver.nc.x*8 + 9] = 0.3;
    //solver.matNs[1][solver.nc.x*8 + 6] = 1.0;
    */

    solver.realloc({256,256}, 2 );
    solver.clearVelocitis();
    solver.clearMater    ();

    solver.matNs[0][solver.nc.x*100 + 50] = 5.0;

    //for(int iy=110; iy<140; iy++){ for(int ix=110; ix<140; ix++){ solver.matNs[0][solver.nc.x*iy + ix] = 1.0; } }
    //for(int iy=120; iy<150; iy++){ for(int ix=120; ix<150; ix++){ solver.matNs[1][solver.nc.x*iy + ix] = 1.0; } }

    for(int iy=100; iy<120; iy++){ for(int ix=100; ix<150; ix++){ solver.matNs[0][solver.nc.x*iy + ix] = 1.0; } }
    //for(int iy=120; iy<150; iy++){ for(int ix=100; ix<150; ix++){ solver.matNs[1][solver.nc.x*iy + ix] = 1.0; } }

    for(int i=0; i<solver.nctot; i++){
        //solver.vs[i] = (Vec2d){1.0,1.0};
        solver.vs[i] = Vec2dZero;
        solver.Vs[i] = 1.0;
        solver.Ts[i] = 1.0;
        solver.ps[i] = 0.0;
        //solver.Ts[i] = 1.0;
    }



}

void TestAppMech2D::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glDisable(GL_LIGHTING);

    glColor3f(1.0,0.0,0.0);

    if(frameCount==1)zoom = 150;
    double dt = 0.003;
    solver.cohesive_pressure = 1.0;

    int perFrame = 1;
    if(bRun){
        for(int itr=0; itr<perFrame; itr++){
            //solver.update(dt);
            //solver.advec(0.001);
            solver.update(dt);
            //solver.cohesion     (dt);
            printf("#### frame %i sumN %g p[100,100]: %g \n", frameCount, solver.sumN, solver.ps[solver.nc.x*100+100] );
            //SDL_Delay(500);
        }
    }

    int ixy=0;
    for(int iy=0; iy<solver.nc.y; iy++){
        for(int ix=0; ix<solver.nc.x; ix++){
            //double c = mpic.moles[ixy]/0.1;
            double n1 = solver.matNs[0][ixy]*0;
            double n3 = solver.matNs[1][ixy]*0;
            double n2 = solver.ps[ixy];

            //double n1 = solver.vs[ixy].x*0.01+0.5;
            //double n2 = solver.Ntot[ixy]*0;
            //double n3 = solver.vs[ixy].y*0.01+0.5;

            //printf( "[%i,%i] %g %g\n", ix,iy, n1, n2 );
            glColor3f( n1, n2, n3 );
            float x=ix-0.5-solver.nc.x*0.5;
            float y=iy-0.5-solver.nc.y*0.5;
            Draw2D::drawRectangle( x, y, x+1, y+1, true );
            ixy++;
        }
    }

    ixy=0;
    for(int iy=0; iy<solver.nc.y; iy++){for(int ix=0; ix<solver.nc.x; ix++){
        glColor3f(1.0,1.0,1.0);
        float x=ix-solver.nc.x*0.5;
        float y=iy-solver.nc.y*0.5;
        Draw2D::drawVecInPos_d( solver.vs[ixy], {x+0.5,y+0.5} );
        ixy++;
    } }

    /*
    solver.cohesive_pressure = -1.0;
    int i0 = solver.nc.x*8 + 7;
    int i1 = solver.nc.x*8 + 8;
    int i2 = solver.nc.x*8 + 9;
    //solver.evalCohesion( solver.matNs[0][i1], solver.matNs[0][i0], solver.ps[i1], solver.ps[i0], dt );
    //solver.evalCohesion( solver.matNs[0][i1], solver.matNs[0][i2], solver.ps[i1], solver.ps[i2], dt );

    solver.evalCohesion( solver.matNs[0][i0], solver.matNs[0][i1], solver.ps[i0], solver.ps[i1], dt );
    solver.evalCohesion( solver.matNs[0][i1], solver.matNs[0][i2], solver.ps[i1], solver.ps[i2], dt );

    glColor3f(1.0,0.0,0.0);
    Draw2D::drawRectangle( 00, 0, 10, solver.matNs[0][i0]*10, true );
    Draw2D::drawRectangle( 10, 0, 20, solver.matNs[0][i1]*10, true );
    Draw2D::drawRectangle( 20, 0, 30, solver.matNs[0][i2]*10, true );
    glColor3f(1.0,1.0,1.0);
    Draw2D::drawRectangle(  0, 0, 10, solver.ps[i0]*10, false );
    Draw2D::drawRectangle( 10, 0, 20, solver.ps[i1]*10, false );
    Draw2D::drawRectangle( 20, 0, 30, solver.ps[i2]*10, false );
    */

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
















