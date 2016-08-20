
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

class MoleculeWorldApp : public AppSDL2OGL {
	public:
    MoleculeWorld world;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	MoleculeWorldApp( int& id, int WIDTH_, int HEIGHT_ );

};

MoleculeWorldApp::MoleculeWorldApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    static double H2O_Q  [3] = { -0.2, +0.1, +0.1 } ;
    static double H2O_R  [3] = { -1.6, +1.2, +1.2 } ;
    static double H2O_pos[6] = { 0.0, 0.0,   +0.5, -1.0,  +0.5, +1.0 } ;

    world.allocate( 1, 8, 3 );

    world.molTypes[0].allocate( 3 );
    //std::copy( world.molTypes[0].Qs, world.molTypes[0].Qs+3,   H2O_Q   );
    //std::copy( world.molTypes[0].Rs, world.molTypes[0].Rs+3,   H2O_R   );
    //std::copy( world.molTypes[0].poss, world.molTypes[0].poss+3, (Vec2d*)H2O_pos );
    memcpy ( world.molTypes[0].Qs,  H2O_Q, 3*sizeof(double) );
    memcpy ( world.molTypes[0].Rs,  H2O_R, 3*sizeof(double) );
    memcpy ( world.molTypes[0].poss,  H2O_pos, 6*sizeof(double) );

    for(int i=0; i<world.molTypes[0].natoms; i++){
        printf( "%i xy=(%3.3f,%3.3f) R=%3.3f Q=%3.3f \n", i, world.molTypes[0].poss[i].x, world.molTypes[0].poss[i].y, world.molTypes[0].Rs[i], world.molTypes[0].Qs[i] );
    }

    srand(154545);
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

    printf( " ==== ==== Frame %i \n", frameCount );
    for( int im=0; im<world.nMols; im++ ){
        int itype = world.mols[im];
        //printf( " mol %i itype %i \n", im, itype );
        MoleculeType2D* moli = &world.molTypes[itype];
		int npi = moli->natoms;
        moli->rigidTransform( world.rot[im], world.pos[im], world.Tps_i );
        for( int ia =0; ia<npi; ia++ ){
            //printf( " xy=(%3.3f,%3.3f) R=%3.3f Q=%3.3f \n", world.Tps_i[ia].x, world.Tps_i[ia].y,  moli->Rs[ia], moli->Qs[ia]  );
            float qc = moli->Qs[ia];
            glColor3f( 0.5-qc, 0.5, 0.5+qc );
            Draw2D::drawCircle_d( world.Tps_i[ia], moli->Rs[ia], 16, false );
        }
    }

};

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
                    //printf( "left button pressed !!!! " );
                    //pickParticle( world.picked );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //world.picked = NULL;
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
















