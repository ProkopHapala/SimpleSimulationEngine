
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

    DEBUG

    int nsamp = 24;
    solver.realloc( 2, 2, 2, nsamp, 1 );
    solver.setRcut( 4.0 );

    // --- Make Geometry
    // initialize atomic positions
    solver.apos[0]=(Vec3d){-1.0,0.0,0.0};
    solver.apos[1]=(Vec3d){+1.0,0.0,0.0};
    // initialize electron positions
    double dy = 0.5;
    //solver.epos[0]=(Vec3d){-0.5,-dy,0.0};  // e[0][0]
    //solver.epos[1]=(Vec3d){-0.5 ,+dy,0.0};  // e[0][1]
    //solver.epos[2]=(Vec3d){+0.5,-dy,0.0};  // e[1][0]
    //solver.epos[3]=(Vec3d){+0.5,+dy,0.0};  // e[1][1]

    solver.epos[0]=(Vec3d){-2.0,-3,0.0};  // e[0][0]
    solver.epos[1]=(Vec3d){-2.0,+3,0.0};  // e[0][1]
    solver.epos[2]=(Vec3d){+2.0,-3,0.0};  // e[1][0]
    solver.epos[3]=(Vec3d){+2.0,+3,0.0};  // e[1][1]


    for(int i=0; i<solver.nBas; i++){
        solver.ecoefs[i] = randf()-0.5;
        printf( "orb[%i].wf[%i] (%g,%g,%g ) C %g \n",  i/solver.perOrb, i%solver.perOrb,  solver.epos[i].x, solver.epos[i].y, solver.epos[i].z, solver.ecoefs[i] );
    }


    DEBUG

    // ----  Make Basis Functions
    double  betaWf=1.8,betaPP=1.8;
    // wave function
    double sumQ = 0.0;
    for(int ir=0; ir<solver.nsampMem; ir++){
        double r    = ((ir-1)*solver.dsamp);  // WARRNING : don't forget [-1] !!!!
        double invr = 1/r;
        double wf  = exp(-betaWf*r);
        //double rho = wf*wf;
        //sumQ      += rho;
        //printf( "" );
        solver.Wfs   [ir] = wf;     // wavefunction
        //solver.Wf2s  [ir] = rho;    // density
        //solver.Vfs   [ir] = sumQ*invr; // potential of density basis function
        solver.PPs[0][ir] = (1-exp(-betaPP*r))*invr; // pseudo-potential of ion (atom)
    }

    plot1.init();
    DataLine2D* line_ref    = new DataLine2D(nsamp, 0.0, solver.dsamp, 0xFFFF00FF, "ref"    ); plot1.add(line_ref   );
    DataLine2D* line_spline = new DataLine2D(nsamp, 0.0, solver.dsamp, 0xFF0000FF, "spline" ); plot1.add(line_spline);
    Spline_Hermite::Sampler<double> spline;
    for(int i=0; i<solver.nsamp; i++){
        double r    = i*solver.dsamp + 0.06;
        spline.prepare( r*solver.idsamp, solver.Wfs );
        line_ref   ->ys[i] = solver.Wfs[i+1];
        line_spline->ys[i] = spline.y();
        printf( "[%i] wf(%g)-> %g | %g \n", i, r, line_spline->ys[i], line_ref->ys[i]  );
    }
    plot1.update();
    plot1.render();


    DEBUG

    //integrateSK( solver.nsamp, 0, solver.nsampI, solver.dsamp, solver.Rcut, dr, solver.Wfs,     solver.Wfs,  solver.ISs, solver.ITs );  // Basis Overlap and Kinetic energy
    //integrateS ( solver.nsamp, 0, solver.nsampI, solver.dsamp, solver.Rcut, dr, solver.Wf2s,    solver.Vfs,  solver.ICs             );  // Coulomb betwween density functions
    //integrateS ( solver.nsamp, 0, solver.nsampI, solver.dsamp, solver.Rcut, dr, solver.PPs[0],  solver.Wf2s, solver.IrhoV[0]        );  // Coulomb with ion pseudopotential

    solver.prepareIntegralTables();

    DEBUG

    // ===========  test project orbital

    GridShape grid;
    //grid.n    = {5,5,5};
    grid.n    = {100,100,100};
    grid.cell = (Mat3d){   10.0,0.0,0.0,     0.0,10.0,0.0,    0.0,0.0,10.0 };
    //grid.pos0 = (Vec3d){-5.0 ,-5.0 ,-5.0  };
    grid.pos0 = (Vec3d){-5.0 ,-5.0 ,-5.0  };
    grid.updateCell();

    DEBUG

    int ng = grid.n.totprod();
    double * orbOnGrid = new double[ ng ];
    //solver.orb2grid( 0, grid, orbOnGrid ); saveXSF( "temp/orb_0.xsf", grid, orbOnGrid, solver.natom, solver.apos, solver.atype );
    //solver.orb2grid( 1, grid, orbOnGrid ); saveXSF( "temp/orb_1.xsf", grid, orbOnGrid, solver.natom, solver.apos, solver.atype );
    solver.orb2grid( 0, grid, orbOnGrid ); saveXSF( "temp/orb_0.xsf", grid, orbOnGrid, solver.perOrb, solver.epos+solver.perOrb*0, 0 );
    solver.orb2grid( 1, grid, orbOnGrid ); saveXSF( "temp/orb_1.xsf", grid, orbOnGrid, solver.perOrb, solver.epos+solver.perOrb*1, 0 );
    delete [] orbOnGrid;

    DEBUG
}

void TestAppCLCFSF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

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
















