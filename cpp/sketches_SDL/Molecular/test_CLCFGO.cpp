
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

#include "CLCFGO.h"

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
    //solver.epos[0]=(Vec3d){-0.5,-dy,0.0};  // e[0][0]
    //solver.epos[1]=(Vec3d){-0.5 ,+dy,0.0};  // e[0][1]
    //solver.epos[2]=(Vec3d){+0.5,-dy,0.0};  // e[1][0]
    //solver.epos[3]=(Vec3d){+0.5,+dy,0.0};  // e[1][1]

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

    //ogl = glGenLists(1);
    //glNewList(ogl,GL_COMPILE);
    //glEndList();


    // ========== Brute Force Overlap


    GridShape grid;
    //grid.n    = {5,5,5};
    grid.n    = {400,100,100};
    grid.cell = (Mat3d){ 40.0,0.0,0.0,  0.0,10.0,0.0,  0.0,0.0,10.0 };
    grid.pos0 = (Vec3d){ 0.0 ,0.0 ,0.0 };
    grid.updateCell();

    solver.ecoef[0] = 1.0;
    solver.ecoef[1] = 0;
    solver.epos [0] = (Vec3d){10.0,5.0,5.0};
    solver.epos [1] = (Vec3d){15.0,5.0,5.0};
    int ng = grid.n.totprod();
    double * orbOnGrid  = new double[ ng ];
    double * orbOnGrid_ = new double[ ng ];
    solver.orb2grid( 0, grid, orbOnGrid );
    double dV = grid.voxelVolume();
    double Q  = dV * VecN::dot (ng, orbOnGrid, orbOnGrid );
    saveXSF( "temp/orb_0.xsf", grid, orbOnGrid, solver.perOrb, solver.epos+solver.perOrb*0, 0 );

    DataLine2D* line_overlap_grid = new DataLine2D(5, 0, 1.0, 0xFF0000FF, "OverlapGrid" ); plot1.add(line_overlap_grid);
    line_overlap_grid->pointStyle = '+'; line_overlap_grid->pointSize = 0.05;
    solver.ecoef[2] = 1.0;
    solver.ecoef[3] = 0;
    for(int i=0; i<line_overlap_grid->n; i++){
        Vec3d p  = solver.epos[0];
        double x = line_overlap_grid->xs[i];
        p.x+=x;
        solver.epos[2] = p;
        solver.orb2grid( 1, grid, orbOnGrid_ );
        double Q  = dV * VecN::dot(ng, orbOnGrid, orbOnGrid_ );
        line_overlap_grid->ys[i] = Q;
    }
    plot1.update();
    plot1.render();
    delete [] orbOnGrid;
    delete [] orbOnGrid_;

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
















