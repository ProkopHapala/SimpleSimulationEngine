
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

/*
#include "Multipoles.h"
#include "PotentialFlow.h"
#include "grids3D.h"
#include "MultipoleGrid.h"
*/

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

//#include "MMFF.h"

#include "RARFF.h"


#define R2SAFE  1.0e-8f

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    RigidAtom     atom1;
    RigidAtomType type1; 

    RARFF grid;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //exit(0);

    atom1.pos  = Vec3dZero;
    atom1.qrot = Quat4dIdentity;

}

void TestAppRARFF::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


    Vec3d bhs[N_BOND_MAX];

    rotateVectors<double>(N_BOND_MAX, atom1.qrot, type1.bh0s, bhs );

    atom1.torq = (Vec3d){0.1,0.0,0.0};
    atom1.moveRotGD(0.8);
    printf( "qrot (%g,%g,%g,%g)\n", atom1.qrot.x, atom1.qrot.y, atom1.qrot.z, atom1.qrot.w );

    glColor3f(0.0,0.0,0.0);
    Draw3D::drawPointCross( atom1.pos, 0.1 );
    for(int i=0; i<N_BOND_MAX; i++){
        //printf( "%i (%g,%g,%g)\n", i, type1.bh0s[i].x, type1.bh0s[i].y, type1.bh0s[i].z );
        //printf( "%i (%g,%g,%g)\n", i, bhs[i].x, bhs[i].y, bhs[i].z );
        Draw3D::drawVecInPos( bhs[i], atom1.pos );
    }

    Draw3D::drawAxis( 1.0);

    //exit(0);
};


void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppRARFF::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppRARFF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRARFF( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















