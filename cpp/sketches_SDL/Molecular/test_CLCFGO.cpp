
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
#include "VecN.h"



#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"

//#include "GUI.h"
#include "Plot2D.h"

#include "integration.h"
#include "AOIntegrals.h"
//#include "AOrotations.h"

#include "Grid.h"
#include "CLCFGO.h"

#include "testUtils.h"

#include "Lingebra.h"
#include "approximation.h"

#include  "Fourier.h"

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"

long timeStart;

// ===================================================
///        test   Wave Fucntion Overlap
// ===================================================

void test_WfOverlap( CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Wavefunction Overlap
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_ISgrid = new DataLine2D( nint, 0, Lmax/nint, 0xFFFF0000, "IS_grid" ); plot1.add(line_ISgrid );
    DataLine2D* line_ISana  = new DataLine2D( 100,  0, 0.1      , 0xFFFF8000, "IS_ana"  );  plot1.add(line_ISana  );
    {
    DEBUG_saveFile1="temp/wf0.xsf";
    DEBUG_saveFile2="temp/wf1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){ solver.epos[0].x=x; solver.orb2grid( 0, grid, f ); };
    auto func2 = [&](GridShape& grid, double* f, double x ){ solver.epos[2].x=x; solver.orb2grid( 1, grid, f ); };
    gridNumIntegral( nint, 0.2, 6.0, Lmax, line_ISgrid->ys, func1, func2 );
    }
    for(int i=0; i<line_ISana->n; i++){
        solver.epos[2].x=line_ISana->xs[i];
        line_ISana->ys[i] = solver.evalOverlap( 0, 1 );
    }
}

// ===================================================
///        test   Density Projection
// ===================================================

void test_ProjectDensity( CLCFGO& solver, Plot2D& plot1 ){
    GridShape grid;
    int mpow = grid.init( 6.0, 0.1, true );
    int ng = grid.n.totprod();
    std::vector<double> grho1(ng,0.0);
    std::vector<double> gorb1(ng,0.0);
    std::vector<double> dgrho(ng,0.0);
    for(int i=0; i<solver.nBas; i++){ solver.ecoef[i]=0; solver.epos[i]=Vec3dZero;  }
    solver.esize[0] = 1.2; solver.ecoef[0] =  1.0; solver.epos[0]=Vec3dZero;
    //solver.esize[1] = 0.7; solver.ecoef[1] = -0.7; solver.epos[1]=(Vec3d){1.0,0.0,0.0};
    Vec3d dip;
    solver.projectOrb( 0, dip, true  );
    solver.projectOrb( 0, dip, false );
    solver.orb2grid  ( 0, grid, gorb1.data() );
    //solver.rho2grid  ( 0, grid, grho1.data() );
    solver.hartree2grid  ( 0, grid, grho1.data() );
    double sumRho = 0;
    double sumWf2 = 0;
    for(int i=0;i<ng;i++){
        double rho = grho1[i];
        double wf  = gorb1[i];
        wf=wf*wf;
        gorb1[i]=wf;
        dgrho[i]=rho-wf;
        sumRho += rho;
        sumWf2 += wf;
    }
    double dV = grid.voxelVolume();
    printf      ( "DEBUG |rho| %g |wf^2| %g \n", sumRho*dV, sumWf2*dV );
    float DEBUG_sc = 0.1;
    int ixy = grid.n.x*grid.n.y*grid.n.z/2 + grid.n.x*grid.n.y/2;
    DataLine2D* line_wf2  = new DataLine2D( grid.n.x, 0, 0.1      , 0xFF0080FF, "wf2"  ); plot1.add(line_wf2  );
    DataLine2D* line_rho  = new DataLine2D( grid.n.x, 0, 0.1      , 0xFFFF00FF, "rho"  ); plot1.add(line_rho  );
    for(int i=0; i<grid.n.x; i++){ line_wf2->ys[i]=gorb1[ixy+i]*DEBUG_sc; line_rho->ys[i]=grho1[ixy+i]*DEBUG_sc; }
    grid.saveXSF( "temp/rho_orb.xsf", gorb1.data() );
    grid.saveXSF( "temp/rho_aux.xsf", grho1.data() );
    grid.saveXSF( "temp/dgrho.xsf"  , dgrho.data() );
}

// ===================================================
///        test   Density Overlap
// ===================================================

void test_DensityOverlap( CLCFGO& solver, Plot2D& plot1 ){
     // ======= Test Orbital Density Overlap
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_IrhoGrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "Irho_grid" ); plot1.add(line_IrhoGrid );
    DataLine2D* line_IrhoAna  = new DataLine2D( 100, 0, 0.1      , 0xFF0080FF, "Irho_ana"  ); plot1.add(line_IrhoAna  );
    //DataLine2D* line_IrhoWf   = new DataLine2D( 100, 0, 0.1      , 0xFFFF00FF, "Irho_Wf"   ); plot1.add(line_IrhoWf  );
    //for(int i=0; i<line_ISana->n; i++){ line_IrhoWf ->ys[i] = sq(line_ISana->ys[i]); }
    for(int i=0; i<line_IrhoAna->n; i++){
        solver.epos[2].x=line_IrhoAna->xs[i];
        solver.projectOrbs();
        line_IrhoAna->ys[i] = solver.DensOverlapOrbPair( 0, 1 );
    }
    { // ---- Test On Grid
    DEBUG_saveFile1="temp/rho0.xsf";
    DEBUG_saveFile2="temp/rho1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){
        solver.epos[0].x=x;
        //solver.projectOrbs();
        solver.orb2grid( 0, grid, f );
        int ntot=grid.n.totprod();
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    auto func2 = [&](GridShape& grid, double* f, double x ){
        solver.epos[2].x=x;
        //solver.projectOrbs();
        int ntot=grid.n.totprod();
        solver.orb2grid( 1, grid, f );
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    gridNumIntegral( nint, 0.2, 6.0, Lmax, line_IrhoGrid->ys, func1, func2, true );
    printf( "DEBUG rho*rho grid %g | | %g \n", line_IrhoGrid->ys[0], line_IrhoAna->ys[0], Gauss::norm3Ds(1) );
    }
}

// =========================================================================
///        test   Electrostatics   ( density * HartreePotential overlap )
// =========================================================================

void test_ElectroStatics( CLCFGO& solver, Plot2D& plot1 ){
    // ======= Test Orbital Density Overlap
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_IElGrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "IEl_grid" ); plot1.add(line_IElGrid );
    DataLine2D* line_IElAna  = new DataLine2D( 100, 0, 0.1       , 0xFF0080FF, "IEl_ana"  ); plot1.add(line_IElAna  );
    //DataLine2D* line_IrhoWf   = new DataLine2D( 100, 0, 0.1      , 0xFFFF00FF, "Irho_Wf"   ); plot1.add(line_IrhoWf  );
    for(int i=0; i<line_IElAna->n; i++){
        solver.epos[2].x=line_IElAna->xs[i];
        solver.projectOrbs();
        //line_IrhoAna->ys[i] = solver.DensOverlapOrbPair( 0, 1 );
        line_IElAna->ys[i] = solver.CoulombOrbPair( 0, 1 );
    }
    { // ---- Test On Grid
    DEBUG_saveFile1="temp/VH0.xsf";
    DEBUG_saveFile2="temp/rho1.xsf";
    auto func1 = [&](GridShape& grid, double* f, double x ){     //    ---- Potencial of Orb_1
        int ntot=grid.n.totprod();
        solver.epos[0].x=x;
        solver.projectOrbs();
        //solver.orb2grid( 0, grid, f );
        //solver.rho2grid( 0, grid, f );
        solver.hartree2grid( 0, grid, f );
        //for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    auto func2 = [&](GridShape& grid, double* f, double x ){     //    ---- Density of Orb_2
        int ntot=grid.n.totprod();
        solver.epos[2].x=x;
        //solver.projectOrbs();
        solver.orb2grid( 1, grid, f );
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi; }
    };
    gridNumIntegral( nint, 0.2, 6.0, Lmax, line_IElGrid->ys, func1, func2, true );
    printf( "DEBUG rho*rho grid %g | | %g \n", line_IElGrid->ys[0], line_IElAna->ys[0], Gauss::norm3Ds(1) );
    }
}


// =========================================================================
///       class   TestAppCLCFSF
// =========================================================================

class TestAppCLCFSF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;

    CLCFGO solver;

    Plot2D plot1;

    int  ogl=0;
    int  fontTex=0;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCLCFSF::TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex     = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    DEBUG

    int nsamp = 40;
    solver.realloc( 2, 2, 2, nsamp, 1 );
    solver.setRcut( 4.0 );

    // --- Make Geometry
    // initialize atomic positions
    solver.apos[0]=(Vec3d){-1.0,0.0,0.0};
    solver.apos[1]=(Vec3d){+1.0,0.0,0.0};
    // initialize electron positions
    double dy = 0.5;
    solver.epos[0]=(Vec3d){-2.0,-3,0.0};  // e[0][0]
    solver.epos[1]=(Vec3d){-2.0,+3,0.0};  // e[0][1]
    solver.epos[2]=(Vec3d){+2.0,-3,0.0};  // e[1][0]
    solver.epos[3]=(Vec3d){+2.0,+3,0.0};  // e[1][1]

    for(int i=0; i<solver.nBas; i++){
        solver.esize[i] = 1.0;
        solver.ecoef[i] = randf()-0.5;
        printf( "orb[%i].wf[%i] (%g,%g,%g ) C %g \n",  i/solver.perOrb, i%solver.perOrb,  solver.epos[i].x, solver.epos[i].y, solver.epos[i].z, solver.ecoef[i] );
    }

    // ======= Test Density projection
    plot1.init();
    plot1.fontTex = fontTex;


    {auto& _=solver;
        for(int i=0; i<_.nBas; i++){ _.ecoef[i]=0; _.epos[i]=Vec3dZero;  }
        _.ecoef[0] = 1.0;
        _.ecoef[2] = 1.0;
        //_.ecoef[1] = +0.2;
        //_.ecoef[3] = +0.3;
    }

    //test_WfOverlap     ( solver, plot1 );
    //test_ProjectDensity( solver, plot1 );
    //test_DensityOverlap( solver, plot1 );   plot1.scaling.y = 30.0;
    test_ElectroStatics( solver, plot1 );

    plot1.update();
    plot1.render();


    DEBUG


}

void TestAppCLCFSF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    glCallList( ogl );

    glDisable(GL_DEPTH_TEST);

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
















