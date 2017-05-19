
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "MMFF.h"
#include "DynamicOpt.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ==========================
// TestAppSoftMolDyn
// ==========================



constexpr int natoms=5, nbonds=4, nang=6, ntors=0;

Vec3d  apos[natoms] = {
{ 0.0, 0.0, 0.0},
{+1.0,+1.0,+1.0},
{-1.0,-1.0,+1.0},
{+1.0,-1.0,-1.0},
{-1.0,+1.0,-1.0}
};   // atomic position

Vec2i  bond2atom[nbonds] = {{0,1},{0,2},{0,3},{0,4}};
//double bond_0   [nbonds] = {1.0,1.2,1.5,1.7};  // [A]
double bond_0   [nbonds] = {1.0,1.0,1.0,1.0};  // [A]

Vec2i  ang2bond [nang]     = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};



/*
constexpr int natoms=3, nbonds=2, nang=1, ntors=0;

Vec3d  apos[natoms] = {
{ 0.0, 0.0, 0.0},
{+1.0,+1.0,+1.0},
{-1.0,-1.0,+1.0},
};   // atomic position

Vec2i  bond2atom[nbonds] = {{0,1},{0,2}};
double bond_0   [nbonds] = {0.8,1.2};  // [A]

Vec2i  ang2bond [nang]   = {{0,1}};
*/


class TestAppSoftMolDyn : public AppSDL2OGL_3D {
	public:

    MMFF       world;
    DynamicOpt opt;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSoftMolDyn::TestAppSoftMolDyn( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    world.apos      = apos;
    world.bond2atom = bond2atom;
    world.ang2bond  = ang2bond;
    world.bond_0    = bond_0;
    world.allocate( natoms, nbonds, nang, ntors );

    opt.bindArrays( 3*natoms, (double*)world.apos, new double[3*natoms], (double*)world.aforce );
    opt.setInvMass( 1.0 );

    world.ang_b2a();

    for(int i=0; i<world.nbonds; i++){
        world.bond_k[i] = 2.0;
    }

    for(int i=0; i<world.nang; i++){
        world.ang_0[i] = {1.0,0.0};
        world.ang_k[i] = 0.5;
        //Vec2i ib = world.ang2bond[i];
        //world.ang2atom [i] = (Vec3i){ world.bond2atom[ib.x].y, world.bond2atom[ib.y].y, world.bond2atom[ib.y].x };
    }

}

void TestAppSoftMolDyn::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	for(int i=0; i<natoms; i++){ world.aforce[i].set(0.0d); }
	world.eval_bonds();
	//world.eval_angles();
	world.eval_angcos();
    //exit(0);

    /*
	for(int i=0; i<world.natoms; i++){ world.aforce[i].add({0.0,-0.1,0.0}); }
	int ipivot = 0;
	world.aforce[ipivot].set(0.0);
    */

    //opt.move_LeapFrog(0.01);
    //opt.move_MDquench();
    //opt.move_FIRE();
    //exit(0);

    printf( "==== frameCount %i\n", frameCount);

    for(int i=0; i<world.natoms; i++){
        glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross(world.apos[i],0.2);
        //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(world.aforce[i],world.apos[i]);
    }
    for(int i=0; i<world.nbonds; i++){
        Vec2i ib = world.bond2atom[i];
        glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine(world.apos[ib.x],world.apos[ib.y]);
    }

};


void TestAppSoftMolDyn::eventHandling ( const SDL_Event& event  ){
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

void TestAppSoftMolDyn::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppSoftMolDyn * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppSoftMolDyn( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















