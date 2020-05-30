
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

    CLCFSF solver;

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
        solver.ecoefs[i] = randf()-0.5;
        printf( "orb[%i].wf[%i] (%g,%g,%g ) C %g \n",  i/solver.perOrb, i%solver.perOrb,  solver.epos[i].x, solver.epos[i].y, solver.epos[i].z, solver.ecoefs[i] );
    }


    DEBUG
    plot1.init();
    plot1.fontTex = fontTex;
    Spline_Hermite::Sampler<double> spline;

    // ----  Make Basis Functions
    double  betaWf=1.8,betaPP=1.8;
    // wave function
    double sumQ = 0.0;
    for(int i=0; i<solver.nsampMem; i++){
        double r    = ((i-1)*solver.dsamp);  // WARRNING : don't forget [-1] !!!!
        double invr = 1/r;
        double wf   = exp(-betaWf*r);
        //double wf   = 1;

        solver.Wfs   [i] = wf;     // wavefunction
        solver.PPs[0][i] = (1-exp(-betaPP*r))*invr; // pseudo-potential of ion (atom)
        printf( "[%i] r: %g wf: %g pp: %g \n", i, r, solver.Wfs[i], solver.Wfs[i] );
    }
    spline.prepare( (solver.Rcut-1e-7)*solver.idsamp, solver.Wfs );
    printf("spline_wf(Rcut=%g|%i+%g) %g \n", (spline.ix+spline.dx)*solver.dsamp, spline.ix,spline.dx, spline.y() );
    //exit(0);

    /*
    DataLine2D* line_ref    = new DataLine2D(nsamp, -solver.dsamp, solver.dsamp, 0xFFFF00FF, "ref"    ); plot1.add(line_ref   );
    DataLine2D* line_spline = new DataLine2D(nsamp,           0.0, solver.dsamp, 0xFF0000FF, "spline" ); plot1.add(line_spline);
    line_spline->pointStyle = '+'; line_spline->pointSize = 0.03;
    for(int i=0; i<solver.nsamp; i++){
        double r    = (i+0.2)*solver.dsamp ;
        spline.prepare( r*solver.idsamp, solver.Wfs );
        line_ref   ->ys[i] = solver.Wfs[i];
        line_spline->xs[i] = r;
        line_spline->ys[i] = spline.y();
        //printf( "[%i] wf(%g)-> %g | %g \n", i, r, line_spline->ys[i], line_ref->ys[i]  );
    }
    plot1.update();
    plot1.render();
    */

    /*
    // ----- test cubic interpolation
    plot1.init();
    plot1.fontTex = fontTex;
    int ntest=20;
    double dtest=1./ntest;
    DataLine2D* line_ref    = new DataLine2D(      4, -1.0, 1.0, 0xFFFF00FF, "ref"    ); plot1.add(line_ref   );
    DataLine2D* line_spline = new DataLine2D(ntest+1,  0.0, dtest, 0xFF0000FF, "spline" ); plot1.add(line_spline);
    Spline_Hermite::Sampler<double> spline;
    int i0=5;
    for(int i=0;i<ntest;i++){
        //solver.Wfs[i0+i-1] = 1.0 + i*0.5 + (i&1)*0.3;
        solver.Wfs[i0+i-1] = 1.0 - i*i*0.25;
        line_ref->ys[i]=solver.Wfs[i0+i-1];
    }
    for(int i=0; i<(ntest+1); i++){
        spline.prepare( i0+i*dtest, solver.Wfs );
        line_spline->ys[i] = spline.y();
        //printf( "[%i] wf(%g)-> %g | %g \n", i, r, line_spline->ys[i], line_ref->ys[i]  );
    }
    plot1.update();
    plot1.render();
    */

    DEBUG

    /*
    // ----- Test 1D integral with symmetry  ::  dot_shifted_sym_( )
    int nfunc  = solver.nsamp+1;
    double dwf = sqrt( 1./((nfunc-1)*2) ); // symmetric
    DataLine2D* line_I1d    = new DataLine2D( nfunc*2, 0, solver.dsamp, 0xFF0000FF, "I   " ); plot1.add(line_I1d   );
    DataLine2D* line_I1dAna = new DataLine2D( nfunc*2, 0, solver.dsamp, 0xFFFF0000, "Iana" ); plot1.add(line_I1dAna);
    double func [nfunc];
    for(int i=0; i<nfunc; i++){
        func[i]=dwf;
    }
    printf( "dwf %g %g %g \n", dwf, dwf*dwf, 1/(dwf*dwf) );
    //dot_shifted_sym_( nfunc, nfunc+5, func, func, 1, 1 );
    //dot_shifted_sym_( nfunc, 10, func, func, 1, 1 ); exit(0);
    for(int i=0; i<nfunc*2; i++){
        //line_I1d   ->ys[i] = dot_shifted_sym( nfunc, i, func, func, 1, 1 );
        line_I1d   ->ys[i] = dot_shifted_sym_( nfunc, i, func, func, 1, 1 );
        line_I1dAna->ys[i] = 1-i*(dwf*dwf);
    }
    plot1.update();
    plot1.render();
    */

    ogl = glGenLists(1);
    glNewList(ogl,GL_COMPILE);
    solver.prepareIntegralTables();
    glEndList();

    DataLine2D* line_wf      = new DataLine2D(solver.nsampMem , -solver.dsamp, solver.dsamp, 0xFFFFFF00, "Wf"     ,solver.Wfs      ); plot1.add(line_wf    );
    DataLine2D* line_overlap = new DataLine2D(solver.nsampIMem, -solver.dsamp, solver.dsamp, 0xFF00FFFF, "Overlap",solver.IOverlap ); plot1.add(line_overlap);
    plot1.update();
    plot1.render();


    // ========== Brute Force Overlap


    GridShape grid;
    //grid.n    = {5,5,5};
    grid.n    = {400,100,100};
    grid.cell = (Mat3d){ 40.0,0.0,0.0,  0.0,10.0,0.0,  0.0,0.0,10.0 };
    grid.pos0 = (Vec3d){ 0.0 ,0.0 ,0.0 };
    grid.updateCell();

    solver.ecoefs[0] = 1.0;
    solver.ecoefs[1] = 0;
    solver.epos  [0] = (Vec3d){10.0,5.0,5.0};
    solver.epos  [1] = (Vec3d){15.0,5.0,5.0};
    int ng = grid.n.totprod();
    double * orbOnGrid  = new double[ ng ];
    double * orbOnGrid_ = new double[ ng ];
    solver.orb2grid( 0, grid, orbOnGrid );
    double dV = grid.voxelVolume();
    double Q  = dV * VecN::dot (ng, orbOnGrid, orbOnGrid );
    printf( "Integrate on grid: Q= %g(%g)  dV %g  ng %i \n", Q, solver.IOverlap[1], dV, ng );
    saveXSF( "temp/orb_0.xsf", grid, orbOnGrid, solver.perOrb, solver.epos+solver.perOrb*0, 0 );

    DataLine2D* line_overlap_grid = new DataLine2D(5, 0, 1.0, 0xFF0000FF, "OverlapGrid" ); plot1.add(line_overlap_grid);
    line_overlap_grid->pointStyle = '+'; line_overlap_grid->pointSize = 0.05;
    solver.ecoefs[2] = 1.0;
    solver.ecoefs[3] = 0;
    for(int i=0; i<line_overlap_grid->n; i++){
        Vec3d p = solver.epos[0];
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

    /*
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
    */
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















