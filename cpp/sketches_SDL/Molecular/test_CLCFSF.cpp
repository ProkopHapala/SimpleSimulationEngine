
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "DrawIso.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Mat3.h"
#include "Mat4.h"
//#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"

#include "Plot2D.h"

#include "integration.h"
#include "AOIntegrals.h"
//#include "AOrotations.h"

#include "CLCFSF.h"

#include "testUtils.h"

#include "Lingebra.h"
#include "approximation.h"

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"

class TestAppCLCFSF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;

    CLCFSF  solver;

    std::vector<Vec3d> ps{ 1,Vec3dZero };


    int ogl=0;
    int ipicked  = -1, ibpicked = -1;
    Plot2D plot1;


    int      fontTex;
    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCLCFSF::TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    //void realloc( natoms, nOrbs, perOrb, int nsamp_, int natypes_ ){

    //double xmax = 8.0;
    //int    n    = 40;
    //double dx   = xmax/n;

    int nsamp = 24;
    solver.realloc( 2, 2, 2, nsamp_, 1 );
    solver.setRcut( 8.0 );

    // --- Make Geometry
    // initialize atomic positions
    solver.apos[0]=(Vec3d){-1.0,0.0,0.0};
    solver.apos[1]=(Vec3d){+1.0,0.0,0.0};
    // initialize electron positions
    double dy = 0.5;
    solver.epos[0]=(Vec3d){0.0,-dy,0.0};  // e[0][0]
    solver.epos[1]=(Vec3d){0.0,+dy,0.0};  // e[0][1]
    solver.epos[2]=(Vec3d){0.0,-dy,0.0};  // e[1][0]
    solver.epos[3]=(Vec3d){0.0,+dy,0.0};  // e[1][1]


    // ----  Make Basis Functions
    double  betaWf=1.8,betaPP=1.8;
    double dr   = 0.1;
    int    nr   = Rmax/dr + 3;
    // wave function
    double* PP = solver.PPs[0];
    double sumQ = 0.0;
    for(int ir=0; ir<nr; ir++){
        double r = ((ir-1)*dr);  // WARRNING : don't forget [-1] !!!!
        double wf = exp(-betaWf*r);
        solver.Wfs [ir] = wf;    // wavefunction
        double rho = wf*wf;
        solver.Wf2s[ir] = rho;   // density
        sumQ += rho;
        solver.Vfs[ir] = sumQ/r; // potential of density basis function
        olver .PPs[ir] = exp(-betaPP*r); // pseudo-potential of ion (atom)
    }
    // pseudo potential

    for(int ir=0; ir<nr; ir++){
        double r = ((ir-1)*dr); // WARRNING : don't forget [-1] !!!!
        s
    }
    // ---- MakeBasis Function Integrals
    integrateSK( nr, 0, nsamp, solver.dsamp, Rmax, dr, solver.Wfs,  solver.Wfs, solver.ISs, solver.IKs );  // Basis Overlap and Kinetic energy
    integrateS ( nr, 0, nsamp, solver.dsamp, Rmax, dr, solver.Wf2s, solver.Vfs, solver.ICs             );  // Coulomb betwween density functions
    integrateS ( nr, 0, nsamp, solver.dsamp, Rmax, dr, solver.PPs,  solver.WFs, solver.IrhoV[0]        );  // Coulomb with ion pseudopotential

}

void TestAppCLCFSF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

};


void TestAppCLCFSF::drawHUD(){
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 20.0,300.00,1.0  );
	plot1.view();

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
















