

/*

==========================================================================
    Test: Compact Linear Combination of Floating Gaussian Orbitals
==========================================================================

### DONE ###
------------
 - tested Wf projection of grid
 - tested Density projection on Aux Basis & on grid
 - tested Electrostatics (Hartree) evaluation
 - tested kinetic energy

### TO DO ###
--------------
 - How to test Fock Exchange ?
    - possibly check overlap density
 - Model for Pauli Energy
    - Perhaps Some inspire by pauli in pseudpotentials
 - Electron-ion interaction  (Electrostatics - Pauli - look how pseudopotentials are done)
 - Make derivatives


 # How test Forces on AuxDensity
 --------------------------------
 1) move wf basis function
 2) calculate some force and energy in on aux-basis
 3) make numerical derivative

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "DrawIso.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Mat3.h"
#include "Mat4.h"
#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"

//#include "GUI.h"
#include "Plot2D.h"

#include "integration.h"
#include "AOIntegrals.h"
//#include "AOrotations.h"

#include "testUtils.h"

#include "Lingebra.h"
#include "approximation.h"

#include  "Fourier.h"

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"

long timeStart;
Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;

#include "Grid.h"
#include "GaussianBasis.h"
#include "CLCFGO.h"
#include "CLCFGO_tests.h"

void testColorOfHash(){
    for(int i=0; i<10; i++){
        Draw::color_of_hash(i);
        Draw2D::drawSimplex( 0,i*1.0, 1, 1.0);
    }
}

void drawSolver(const CLCFGO& solver, float fsc ){
    glColor3f(0.,0.,0.);
    for(int i=0; i<solver.natom; i++){
        Vec3d p = solver.apos[i];
        Draw3D::drawPointCross( p, solver.aPsize[i] );
        Draw3D::drawVecInPos( solver.aforce[i]*fsc, p );
    }
    for(int io=0; io<solver.nOrb; io++){
        Draw::color_of_hash(io*15446+7545);
        for(int j=0; j<solver.perOrb; j++){
            int i = io*solver.perOrb+j;
            Vec3d p = solver.epos[i];
            Draw3D::drawSphereOctLines(16, solver.esize[i], p, Mat3dIdentity, false );
            Draw3D::drawVecInPos( solver.efpos[i]*fsc, p );
        }
    }
}


// =========================================================================
///       class   TestAppCLCFSF
// =========================================================================

class TestAppCLCFSF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;
    double dt;

    CLCFGO solver;

    Plot2D plot1;

    int  ogl=0;
    int  fontTex=0;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );


    void test_RhoDeriv();

    TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ );

};


TestAppCLCFSF::TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex     = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    /*
    // ---- 2 atoms, 2 orbs, 2 basis per orb
    //      natom  nOrb perOrb natypes
    solver.realloc( 2, 2, 2, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        _.setAtom( 0,   1,  (Vec3d){  2.0, 0.,0.} );
        _.setAtom( 0,   1,  (Vec3d){ -2.5, 0.,0.} );
        _.setElectron( 0,0, (Vec3d){  0.0, +0.2, 0.}, 1.0, -0.3 );
        _.setElectron( 0,1, (Vec3d){  1.0, -0.2, 0.}, 1.0,  0.7 );
        _.setElectron( 1,0, (Vec3d){ -2.0, -1.0, 0.}, 1.0,  1.0 );
        _.setElectron( 1,1, (Vec3d){ -2.0, +1.0, 0.}, 1.0,  1.0 );
    }
    */

    /*
    // ---- 1 H atoms, 1 orbs, 1 basis per orb
    //      natom  nOrb perOrb natypes
    solver.realloc( 1, 1, 1, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        //_.setAtom( 0,   1,  (Vec3d){  0.0, 0.0, 0.0 } );
        //_.setAtom    ( 0,   (Vec3d){  0.0, 0.0, 0.0 }, 4.0, 0.5, 0.2, 1000 );
        _.setAtom    ( 0,   (Vec3d){  0.0, 0.0, 0.0 }, 1.0,  0.1, 1.0, 0. );
        _.setElectron( 0,0, (Vec3d){  1.0, 0.0, 0.0 }, 0.25, 1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- 2 basis electron in soft potential   (1a,1o,2b)
    //      natom  nOrb perOrb natypes
    solver.realloc( 1, 1, 2, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        _.setAtom    ( 0,   (Vec3d){   0.0, 0.0, 0.0 }, 1.0,  0.9, 1.0, 0. );
        _.setElectron( 0,0, (Vec3d){  +1.0, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 0,1, (Vec3d){  -1.0, 0.0, 0.0 }, 0.5, -1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- H2+ molecule ion - 2 basis electron in 2 atom potential   (2a,1o,2b)
    //      natom  nOrb perOrb natypes
    solver.realloc( 2, 1, 2, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        //_.setAtom( 0,   1,  (Vec3d){  0.0, 0.0, 0.0 } );
        //_.setAtom    ( 0,   (Vec3d){  0.0, 0.0, 0.0 }, 4.0, 0.5, 0.2, 1000 );
        _.setAtom    ( 0,   (Vec3d){  +1.0, 0.0, 0.0 }, 1.0,  0.9, 1.0, 0. );
        _.setAtom    ( 1,   (Vec3d){  -1.0, 0.0, 0.0 }, 1.0,  0.9, 1.0, 0. );
        _.setElectron( 0,0, (Vec3d){  +0.7, 0.0, 0.0 }, 0.5, +1.0 );
        //_.setElectron( 0,1, (Vec3d){  -0.7, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 0,1, (Vec3d){  -0.7, 0.0, 0.0 }, 0.5, -1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- Free 2 basis electron    (0a,1o,2b)
    //      natom  nOrb perOrb natypes
    solver.realloc( 0, 1, 2, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        _.setElectron( 0,0, (Vec3d){  +0.7, 0.0, 0.0 }, 0.5, +1.0 );
        //_.setElectron( 0,1, (Vec3d){  -0.7, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 0,1, (Vec3d){  -0.7, 0.0, 0.0 }, 0.5, -1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- Free 2 electron 1 basis each    (0a,2o,1b)
    //      natom  nOrb perOrb natypes
    solver.realloc( 0, 2, 1, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        _.setElectron( 0,0, (Vec3d){  +0.3, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 1,0, (Vec3d){  -0.3, 0.0, 0.0 }, 0.5, +1.0 );
    }
    dt = 0.001;
    */


    // ---- 2 electron 1 basis each in Soft Atom potential   (1a,2o,1b)
    //      natom  nOrb perOrb natypes
    solver.realloc( 1, 2, 1, 1 );
    solver.setRcut( 4.0 );
    solver.setDefaultValues();
    {CLCFGO& _=solver;
        _.setAtom    ( 0,   (Vec3d){  0.0, 0.0, 0.0 }, 1.0,  0.1, 1.0, 0. );
        _.setElectron( 0,0, (Vec3d){  +0.3, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 1,0, (Vec3d){  -0.3, 0.0, 0.0 }, 0.5, +1.0 );
    }
    dt = 0.001;


    //solver.loadFromFile( "data/H2.fgo", true);
    dt = 0.001;
    //exit(0);

    solver.turnAllSwitches(false);
    solver.bNormalize     = 1;
    solver.bEvalAE        = 1;
    //solver.bEvalAECoulomb = 1;
    //solver.bEvalCoulomb   = 1;
    solver.bEvalPauli     = 1;
    solver.bEvalKinetic   = 1;

    solver.bOptEPos = 1;
    solver.bOptSize = 1;

    solver.iPauliModel = 0;

    /*
    //solver.bNormalize     = false;
    solver.bEvalKinetic   = false;
    //solver.bEvalCoulomb   = false;
    solver.bEvalPauli     = false;
    solver.bEvalExchange  = false;
    //solver.bEvalAECoulomb = false;
    solver.bEvalAEPauli   = false;
    //solver.bEvalAE      = false;
    solver.bEvalAA        = false;
    solver.iPauliModel    = 1;

    solver.bOptAtom = false;
    //solver.bOptEPos = false;
    solver.bOptSize = false;
    solver.bOptCoef = false;
    */

    // ======= Test Density projection
    plot1.init();
    plot1.fontTex = fontTex;
    //GridShape grid; grid.init( 5.0, 0.2, false);
    //solver.orb2xsf( grid, 0, "temp/orb0.xsf" );

    testDerivsTotal( 100, -1.0, 0.05, solver, plot1 );
    plot1.xsharingLines(1, 100, -3.0, 0.1 );

    plot1.scaling.y=0.05;
    plot1.update();
    plot1.render();

}

void TestAppCLCFSF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable( GL_DEPTH_TEST );

    //testColorOfHash();
    double E = solver.eval();
    float F2 = solver.moveGD(dt);
    printf( "frame[%i] E %g |F| %g \n", frameCount, E, sqrt(F2) );

    drawSolver( solver, 1.0 );
    plotOrb( solver, plot1.lines[3], 0, (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, 100.0 );
    testDerivsTotal( 0, 0, 0, solver, plot1, 0 );
    plot1.update();
    plot1.render();
    //delay=500;

    //glCallList( ogl );
    //glDisable(GL_DEPTH_TEST);
    plot1.view();
};


void TestAppCLCFSF::drawHUD(){
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
	//glTranslatef( 100.0,100.0,0.0 );
	//glScalef    ( 20.0,300.00,1.0  );
	//plot1.view();

}


void TestAppCLCFSF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppCLCFSF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppCLCFSF( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}


