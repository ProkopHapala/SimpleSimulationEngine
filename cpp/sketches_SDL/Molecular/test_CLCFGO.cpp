
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

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"

long timeStart;


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
    //GUI_fontTex = fontTex;
    //void realloc( natoms, nOrbs, perOrb, int nsamp_, int natypes_ ){

    //double xmax = 8.0;
    //int    n    = 40;
    //double dx   = xmax/n;

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


    DEBUG
    plot1.init();
    plot1.fontTex = fontTex;



    /*
    // ======= Test Orbital Wavefunction Overlap
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_ISgrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "IS_grid" ); plot1.add(line_ISgrid );
    DataLine2D* line_ISana  = new DataLine2D( 100,  0, 0.1      , 0xFF0080FF, "IS_ana"  );  plot1.add(line_ISana  );

    {auto& _=solver;
        for(int i=0; i<_.nBas; i++){ _.ecoef[i]=0; _.epos[i]=Vec3dZero;  }
        _.ecoef[0] = 1.0;
        _.ecoef[2] = 1.0;
        _.ecoef[1] = +0.5;
        _.ecoef[3] = -0.7;
    }
    {
    auto func1 = [&](GridShape& grid, double* f, double x ){ solver.epos[0].x=x; solver.orb2grid( 0, grid, f ); };
    auto func2 = [&](GridShape& grid, double* f, double x ){ solver.epos[2].x=x; solver.orb2grid( 1, grid, f ); };
    gridNumIntegral( nint, 0.1, 6.0, Lmax, line_ISgrid->ys, func1, func2 );
    }
    for(int i=0; i<line_ISana->n; i++){
        solver.epos[2].x=line_ISana->xs[i];
        line_ISana->ys[i] = solver.evalOverlap( 0, 1 );
    }
    plot1.update();
    plot1.render();
    */


    // ======= Test Orbital Density Overlap
    int    nint = 10;
    double Lmax = 5.0;
    DataLine2D* line_IrhoSgrid = new DataLine2D( nint, 0, Lmax/nint, 0xFF0000FF, "IS_grid" ); plot1.add(line_IrhoSgrid );
    DataLine2D* line_IrhoSana  = new DataLine2D( 10,  0, 0.5       , 0xFF0080FF, "IS_ana"  ); plot1.add(line_IrhoSana  );
    {auto& _=solver;
        for(int i=0; i<_.nBas; i++){ _.ecoef[i]=0; _.epos[i]=Vec3dZero;  }
        _.ecoef[0] = 1.0;
        _.ecoef[2] = 1.0;
        //_.ecoef[1] = +0.2;
        //_.ecoef[3] = +0.3;
    }
    for(int i=0; i<line_IrhoSana->n; i++){
        solver.epos[2].x=line_IrhoSana->xs[i];
        solver.projectOrbs();
        line_IrhoSana->ys[i] = solver.DensOverlapOrbPair( 0, 1 );
    }
    {
    double DEBUG_scale = 3.0;
    auto func1 = [&](GridShape& grid, double* f, double x ){
        solver.epos[0].x=x;
        solver.projectOrbs();
        solver.orb2grid( 0, grid, f );
        int ntot=grid.n.totprod();
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi*DEBUG_scale ; }
    };
    auto func2 = [&](GridShape& grid, double* f, double x ){
        solver.epos[2].x=x;
        solver.projectOrbs();
        int ntot=grid.n.totprod();
        solver.orb2grid( 1, grid, f );
        for(int i=0; i<ntot; i++){ double fi=f[i]; f[i]=fi*fi*DEBUG_scale ; }
    };
    gridNumIntegral( nint, 0.1, 6.0, Lmax, line_IrhoSgrid->ys, func1, func2 );
    }

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
















