/// @file @brief This program demonstrates 2D global optimization for molecular systems using `DynamicOpt.h` and `MoleculeWorld2D.h`. It visualizes the iterative search for minimum energy configurations of molecules, allowing the user to observe the system settling into stable arrangements. The optimization process can be toggled with the `SPACE` key, and the system can be perturbed with the mouse.


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

#include "Body2D.h"
#include "MoleculeWorld2D.h"
//#include <algorithm>


/*

ToDo :
    - Something wrong with charges - molecules repel each orther at long distance with non-zero charges

*/





class MoleculeWorldApp : public AppSDL2OGL {
	public:
    bool recalc = true;
    MoleculeWorld world;

    int pickedMol = -1;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	int pickMol( const Vec2d& mpos );

	MoleculeWorldApp( int& id, int WIDTH_, int HEIGHT_ );

};

MoleculeWorldApp::MoleculeWorldApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    world.allocate( 1, 5, 3 );

    static double H2O_Q  [3] = { +0.4, -0.4, 0.0 } ;
    //static double H2O_Q  [3] = { -0.2, +0.4, -0.2 } ;
    static double H2O_R  [3] = { -1.6, +1.2, +1.2 } ;
    static double H2O_pos[6] = { 0.0, 0.0,   +0.5, -1.0,  +0.5, +1.0 } ;
    world.molTypes[0].allocate( 3 );
    memcpy ( world.molTypes[0].Qs,  H2O_Q, 3*sizeof(double) );
    memcpy ( world.molTypes[0].Rs,  H2O_R, 3*sizeof(double) );
    memcpy ( world.molTypes[0].poss,  H2O_pos, 6*sizeof(double) );


/*
    static double H2O_Q  [1] = { 0.0 } ;
    static double H2O_R  [1] = { 1.0 } ;
    static double H2O_pos[2] = { 0.0, 0.0 } ;
    world.molTypes[0].allocate( 1 );
    memcpy ( world.molTypes[0].Qs,      H2O_Q, 1*sizeof(double) );
    memcpy ( world.molTypes[0].Rs,      H2O_R, 1*sizeof(double) );
    memcpy ( world.molTypes[0].poss,    H2O_pos, 2*sizeof(double) );
*/

    for(int i=0; i<world.molTypes[0].natoms; i++){
        printf( "%i xy=(%3.3f,%3.3f) R=%3.3f Q=%3.3f \n", i, world.molTypes[0].poss[i].x, world.molTypes[0].poss[i].y, world.molTypes[0].Rs[i], world.molTypes[0].Qs[i] );
    }

    srand(1545551);
    for( int i=0; i<world.nMols; i++ ){
        world.mols[i] = 0;
        world.notFixed[i] = true;
        world.rot[i]  = randf(-M_PI,M_PI);
        world.pos[i].set( randf(-10.0,10.0), randf(-10.0,10.0) );
    }

    //exit(0);
}

void MoleculeWorldApp::draw(){
    glClearColor( 0.9f, 0.9f, 0.9f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( " ==== ==== Frame %i \n", frameCount );


    if( recalc ){
        VecN::set( 3*world.nMols, 0.0, world.optimizer.force );
        world.forcesMolecules( );
        if(pickedMol>=0){
            Vec2d Fout,pmouse;
            pmouse.set(mouse_begin_x,mouse_begin_y);
            stringForce( world.pos[pickedMol], pmouse, 10.0, Fout );
            world.fpos[pickedMol].add( Fout );
            glColor3f(0.2f,0.5f,0.2f);
            Draw2D::drawLine_d( world.pos[pickedMol], pmouse );
        }
        //world.optimizer.dt = 0.01;
        //world.optimizer.move_LeapFrog();
        world.optimizer.move_MDquench();
        //world.optimizer.move_FIRE();
        //recalc=false;
    }


    for( int im=0; im<world.nMols; im++ ){
        int itype = world.mols[im];
        //printf( " mol %i itype %i \n", im, itype );
        MoleculeType2D* moli = &world.molTypes[itype];
		int npi = moli->natoms;

        Draw2D::drawPointCross_d( world.pos [im], 0.1 );
        Draw2D::drawVecInPos_d  ( world.fpos[im], world.pos [im] );

        //printf( "%i (%3.3f,%3.3f) \n", im, world.fpos[im].x, world.fpos[im].y );

        moli->rigidTransform( world.rot[im], world.pos[im], world.Tps_i );
        for( int ia =0; ia<npi; ia++ ){
            //printf( " xy=(%3.3f,%3.3f) R=%3.3f Q=%3.3f \n", world.Tps_i[ia].x, world.Tps_i[ia].y,  moli->Rs[ia], moli->Qs[ia]  );
            float qc = moli->Qs[ia];
            glColor3f( 0.5-qc, 0.5, 0.5+qc );
            Draw2D::drawCircle_d    ( world.Tps_i[ia], moli->Rs[ia]*0.8, 16, false );
        }
    }

};

int MoleculeWorldApp::pickMol( const Vec2d& mpos ){
    int    imin  = -1;
    double r2min = 1e+300;
    for(int i=0; i<world.nMols; i++){
        Vec2d d; d.set_sub( mpos, world.pos[i] );
        double r2 = d.norm2();
        if (r2<r2min){ r2min=r2; imin=i; }
    }
    return imin;
}

void MoleculeWorldApp::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
    //world.anchor.set( mouse_begin_x, mouse_begin_y );
}


void MoleculeWorldApp::eventHandling ( const SDL_Event& event  ){
    //printf( "MoleculeWorld2DApp::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //world.pos[0].set( mouse_begin_x, mouse_begin_y );
                    //recalc = true;
                    pickedMol = pickMol( {mouse_begin_x,mouse_begin_y} );
                break;
                case SDL_BUTTON_RIGHT:
                    world.pos[0].set( mouse_begin_x, mouse_begin_y );
                    //recalc = true;
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //world.picked = NULL;
                    pickedMol  = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void MoleculeWorldApp::drawHUD(){

}

// ===================== MAIN

MoleculeWorldApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new MoleculeWorldApp( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
