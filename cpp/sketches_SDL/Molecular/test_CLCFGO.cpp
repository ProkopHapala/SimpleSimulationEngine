

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

int i_DEBUG = 0;

#include "Grid.h"
#include "GaussianBasis.h"
#include "CLCFGO.h"
#include "CLCFGO_tests.h"

int  fontTex=0;

void testColorOfHash(){
    for(int i=0; i<10; i++){
        Draw::color_of_hash(i);
        Draw2D::drawSimplex( 0,i*1.0, 1, 1.0);
    }
}

int orbColor(int io){
    return hash_Wang( hash_Wang( io*15446+7545 ) );
    //Draw::color_of_hash(io*15446+7545,clr);
    //return clr;
}


void drawff_atoms(const CLCFGO& ff, float fsc=1.0, float asc=0.5 ){
    glEnable(GL_DEPTH_TEST);
    glColor3f(0.,0.,0.);
    //glDisable(GL_DEPTH_TEST);
    //glEnable(GL_BLEND);
    //glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    for(int i=0; i<ff.natom; i++){
        Vec3d p = ff.apos[i];
        //Draw3D::drawPointCross( p, ff.aPsize[i]*asc );
        glColor3f(0.,0.,0.);
        Draw3D::drawPointCross( p, ff.aPars[i].z*asc );
        glColor3f( 1.0f,0.0f,0.0f );
        Draw3D::drawVecInPos( ff.aforce[i]*fsc, p );
    }
}

void drawff_wfs(const CLCFGO& ff, int oglSph, float fsc=1.0, float asc=0.5, int alpha=0x15000000 ){
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    char str[256];
    for(int io=0; io<ff.nOrb; io++){
        //Vec3f clr;
        //Draw::color_of_hash(io*15446+7545,clr);
        //glColor4f(clr.x,clr.y,clr.z,0.1);
        for(int j=0; j<ff.perOrb; j++){
            int i = io*ff.perOrb+j;
            Vec3d p = ff.epos[i];
            //float alpha=0.1;
            //if(ff.espin[i]>0){ glColor4f(0.0,0.0,1.0, alpha); }else{ glColor4f(1.0,0.0,0.0, alpha); };
            //int alphaMax=200;
            //int alpha = alphaMax*fabs(ff.rhoQ[i]); if(alpha>alphaMax)alpha=alphaMax; alpha<<=24;
            int c = orbColor(io);
            Draw  ::setRGBA( (c&0x00FFFFFF)|alpha  ); Draw3D::drawShape( oglSph, ff.epos[i], Mat3dIdentity*ff.esize[i],  false );
            //Draw  ::setRGBA(  c                    ); Draw3D::drawSphereOctLines(16, ff.esize[i], p, Mat3dIdentity, false );
            Draw  ::setRGBA(  c                    ); Draw3D::drawPointCross( p, 0.01 );
            glColor3f( 1.0f,0.0f,0.0f );  Draw3D::drawVecInPos( ff.efpos[i]*fsc, p );

            Draw  ::setRGBA( orbColor(io) );
            //sprintf(str, "%02i_%02i", io, j  );
            sprintf(str, "%3.3f", ff.ecoef[i]  );
            Draw3D::drawText(str, p, fontTex, 0.02,  0 );
        }
    }
}

void drawff_rho(const CLCFGO& ff, int oglSph, float fsc=1.0, int alpha=0x15000000 ){
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    char str[256];
    for(int io=0; io<ff.nOrb; io++){
        int i0 = ff.getRhoOffset(io);
        for(int ii=0; ii<ff.onq[io]; ii++){
            int i   = ii+i0;
            Vec3d p = ff.rhoP[i];

            //int alphaMax=200;
            //int alpha = alphaMax*fabs(ff.rhoQ[i]); if(alpha>alphaMax)alpha=alphaMax; alpha<<=24;
            int c = orbColor(io);
            Draw  ::setRGBA( (c&0x00FFFFFF)| alpha ); Draw3D::drawShape( oglSph, p, Mat3dIdentity*ff.rhoS[i],  false );
            //Draw  ::setRGBA(  c                    ); Draw3D::drawSphereOctLines(16, ff.rhoS[i], p, Mat3dIdentity, false );
            Draw  ::setRGBA(  c                    ); Draw3D::drawPointCross( p, 0.01 );

            glColor3f( 1.0f,0.0f,0.0f );
            Draw3D::drawVecInPos( ff.rhofP[i]*fsc, p );

            /*
            Draw  ::setRGBA( orbColor(io) );
            //sprintf(str, "%02i_%02i_%3.3f", io, ii, ff.ecoef[i]  );
            sprintf(str, "%3.3f", ff.rhoQ[i]  );
            Draw3D::drawText(str, p, fontTex, 0.02,  0 );
            */
        }
    }
}


// =========================================================================
///       class   TestAppCLCFSF
// =========================================================================

class TestAppCLCFSF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;
    int iter=0;
    int perFrame = 10;
    bool bRun = false;
    double dt;
    //double E,dt;

    CLCFGO ff;
    Plot2D plot1;
    int nOrbPlot = 2;

    bool bDrawPlots   = true;
    bool bDrawObjects = true;
    bool bDrawAtoms   = true;
    bool bDrawWfs     = true;
    bool bDrawRho     = false;
    bool bPlotDens    = false;

    int  ogl=0,oglSph=0;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    void test_RhoDeriv();
    void initff();
    void viewPlots();

    TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ );

};


TestAppCLCFSF::TestAppCLCFSF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //ff.setDefaultValues( );
    //ff.loadFromFile( "data/e_2g_1o.fgo"            );
    //ff.loadFromFile( "data/e2_1g_2o_singlet.fgo" );
    //ff.loadFromFile( "data/e2_1g_2o_singlet.fgo" );
    //ff.loadFromFile( "data/e2_1g_2o_triplet.fgo" );
    //ff.loadFromFile( "data/H_1g.fgo"             );
    //ff.loadFromFile( "data/H_2g.fgo"             );
    //ff.loadFromFile( "data/H_2g_split.fgo"       );
    //ff.loadFromFile( "data/H_2g_problem.fgo"     );
    //ff.loadFromFile( "data/H_2g_problem_sym.fgo"   );
    //ff.loadFromFile( "data/H_3g.fgo"   );
    //ff.loadFromFile( "data/H_5g.fgo"   );
    //ff.loadFromFile( "data/Hanti_2g_anti.fgo"    );
    //ff.loadFromFile( "data/He_singlet.fgo"       );
    //ff.loadFromFile( "data/He_triplet.fgo"       );
    //ff.loadFromFile( "data/H2_1g_2o.fgo"         );
    //ff.loadFromFile( "data/H2.fgo"               );
    //ff.loadFromFile( "data/H2_2g_half.fgo"       );
    //ff.loadFromFile( "data/H2_3g_half.fgo"       );
    //ff.loadFromFile( "data/H2_4g_half.fgo"       );
    //ff.loadFromFile( "data/He_2g_triplet_sym.fgo"  );
    //ff.loadFromFile( "data/He_2g_triplet_asym.fgo" );
    //ff.loadFromFile( "data/He_2g_triplet_sym1.fgo" );
    //ff.loadFromFile( "data/Li_2g.fgo"            );
    //ff.loadFromFile( "data/Li_3g.fgo"            );
    //ff.loadFromFile( "data/Li_4g.fgo"            );
    //ff.loadFromFile( "data/B_2g_triplet.fgo"     );
    //ff.loadFromFile( "data/B_2g_triplet_asym.fgo"     );
    //ff.loadFromFile( "data/C_1g.fgo"             );
    //ff.loadFromFile( "data/C_2g_triplet.fgo"     );
    //ff.loadFromFile( "data/C_2g_triplet-.fgo"      );
    //ff.loadFromFile( "data/C_2g_o1.fgo"          );
    ff.loadFromFile( "data/C_2g_problem.fgo"      );
    //ff.loadFromFile( "data/N2.fgo"               );
    //ff.loadFromFile( "data/O2.fgo"               );
    //ff.loadFromFile( "data/O2_half.fgo"          );
    //ff.loadFromFile( "data/H2O_1g_8o.fgo"        );

    //exit(0);

    //ff.turnAllSwitches(false);
    ff.turnAllEvalSwitches(false);
    ff.bNormalize     = 1;
    ff.bNormForce     = 1;
    ff.bEvalKinetic   = 1;
    ff.bEvalCoulomb   = 1;
    ff.bEvalPauli     = 1;
    ff.bEvalAE        = 1;
    ff.bEvalAECoulomb = 1;
    ff.bEvalAEPauli   = 1;
    ff.bEvalAA        = 1;

    //ff.bNormalize     = 0;
    //ff.bNormForce     = 0;
    //ff.bEvalKinetic   = 0;
    //ff.bEvalCoulomb   = 0;
    //ff.bEvalPauli     = 0;
    //ff.bEvalAE        = 0;
    //ff.bEvalAECoulomb = 0;
    //ff.bEvalAEPauli   = 0;

    //ff.bEvalAA        = 0;

    //ff.bOptAtom = 1;
    //ff.bOptEPos = 1;
    //ff.bOptSize = 1;

    ff.bOptAtom = 0;
    //ff.bOptEPos = 0;
    //ff.bOptSize = 0;
    //ff.bOptCoef = 0;
    //ff.ofix[0] = 1;

    //ff.iPauliModel = 2;
    //ff.iPauliModel = 0;
    dt = 0.001;

    //bDrawWfs  = 0; bDrawRho  = 1;   // Plot Density blobs instead of wavefunctions


    ff.printSetup();
    ff.printAtoms();
    ff.printElectrons();

    double E = ff.eval();
    ff.forceInfo();
    printf( "E %g | Ek %g Eee,p,ex(%g,%g,%g) Eae,p(%g,%g) Eaa %g | s[0] %g \n",E, ff.Ek, ff.Eee,ff.EeePaul,ff.EeeExch, ff.Eae,ff.EaePaul, ff.Eaa, ff.esize[0] );
    //exit(0);

    printf( "epos[0].x " ); getNumDeriv( ff, ff.epos [0].x, ff.efpos [0].x, 0.0001, 0 );
    printf( "esize[0]  " ); getNumDeriv( ff, ff.esize[0]  , ff.efsize[0]  , 0.0001, 0 );
    printf( "ecoef[0]  " ); getNumDeriv( ff, ff.ecoef[0]  , ff.efcoef[0]  , 0.0001, 0 );
    //ff.eval();
    //ff.printElectrons();
    exit(0);

    // ======= Test Density projection
    plot1.init();
    plot1.fontTex = fontTex;
    //GridShape grid; grid.init( 5.0, 0.2, false);
    //ff.orb2xsf( grid, 0, "temp/orb0.xsf" );

    //testDerivsTotal( 100, -1.0, 0.05, ff, plot1 );
    //plot1.xsharingLines(1, 100, -3.0, 0.1 );
    //plot1.xsharingLines( 2, 100,   -3.0, 0.1 );


    nOrbPlot=_min(5,ff.nOrb);
    printf( " nOrbPlot %i ff.nOrb %i \n", nOrbPlot, ff.nOrb );
    plot1.add( new DataLine2D( 100, -3.0, 0.1, 0xFF0000FF, "Vatom" ) );
    char str[32];
    for(int io=0; io<nOrbPlot; io++){
        sprintf(str,"orb%03i",io);
        plot1.add( new DataLine2D( 100, -3.0, 0.1, orbColor(io)|0xFF000000, str ) );
    }
    //plot1.add( new DataLine2D( 100, -3.0, 0.1, 0xFFFF0000, "Orb1"     ) );
    //plot1.add( new DataLine2D( 100, -3.0, 0.1, 0xFFFF8000, "Orb2"     ) );


    /*
    ff.epos [0]=Vec3dZero;
    //ff.esize[0]=1.0;
    ff.esize[0]=0.5;
    //test_Poisson( ff, 6.0, 0.1 );
    test_Poisson( ff, 0, 4.0, 0.05, 0,0,  true, false, true );
    //test_ElectroStaticsBrute( ff, plot1 );
    //test_ElectroStatics( ff, plot1 );
    exit(0);
    */

    plot1.scaling.y=0.05;
    plot1.update();
    plot1.render();

    oglSph=Draw::list(oglSph);
    Draw3D::drawSphere_oct(4,1.0d,(Vec3d){0.,0.,0.});
    glEndList();

    bRun = false;

    //dt = 0.000001;
    //dt = 0.001;

}

void TestAppCLCFSF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable( GL_DEPTH_TEST );

    //dt = 0.001;
    if(bRun){
        double F2;
        for(int itr=0; itr<perFrame; itr++){
            //testColorOfHash();
            ff.eval();
            ff.forceInfo();
            //printf( "frame[%i] E %g | Ek %g Eee,p(%g,%g) Eae,p(%g,%g) Eaa %g \n", frameCount, E, ff.Ek, ff.Eee,ff.EeePaul,  ff.Eae,ff.EaePaul, ff.Eaa );
            F2 = ff.moveGD(dt);
            iter++;
        }

        printf( "frame[%i] E %g |F| %g \n", frameCount, ff.Etot, sqrt(F2) );
        //printf( "frame[%i] |F| %g E %g | Ek %g Eee,p,ex(%g,%g,%g) Eae,p(%g,%g) Eaa %g | s[0] %g \n", frameCount, sqrt(F2), E, ff.Ek, ff.Eee,ff.EeePaul,ff.EeeExch, ff.Eae,ff.EaePaul, ff.Eaa, ff.esize[0] );
        //printf( "frame[%i] E %g pa[0](%g,%g,%g) pe[0](%g,%g,%g) s %g \n", frameCount, E, ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z,  ff.esize[0] );
        //printf( "frame[%i] E %g lHH %g lH1e1 %g se1 %g \n", frameCount, E, (ff.apos[0]-ff.apos[1]).norm(),   (ff.apos[0]-ff.epos[0]).norm(), ff.esize[0] );

        if(F2<1e-6){
            bRun=false;
            ff.printAtoms();
            ff.printElectrons();
        }
    }

    if(bDrawObjects){
        float fsc=0.01;
        if(bDrawAtoms) drawff_atoms( ff,         fsc, 0.2 );
        if(bDrawWfs  ) drawff_wfs  ( ff, oglSph, fsc      );
        if(bDrawRho  ) drawff_rho  ( ff, oglSph, fsc      );
    }
    if(bDrawPlots){ viewPlots(); }

};


void TestAppCLCFSF::drawHUD(){
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
	//glTranslatef( 100.0,100.0,0.0 );
	//glScalef    ( 20.0,300.00,1.0  );
	//plot1.view();


	//gui.draw();

    glTranslatef( 10.0,HEIGHT-20.0,0.0 );
	glColor3f(0.5,0.0,0.3);

	int nstr=2048;
	char str[nstr];
	char* s=str;
	s+=ff.Eterms2str(s);
	ff.orbs2str(s);
    Draw::drawText( str, fontTex, fontSizeDef, {100,20} );


}

void TestAppCLCFSF::viewPlots(){
        plotAtomsPot( ff, plot1.lines[0],    (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, 1.0, 0.1 );
        //plotAtomsPot( ff, plot1.lines[0],    (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, 1.0, ff.esize[0]*M_SQRT1_2 );
        for(int io=0; io<nOrbPlot; io++){
            //plot1.add( new DataLine2D( 100, -3.0, 0.1, orbColor(io)|0xFF000000, str ) );
            //plotOrb     ( ff, plot1.lines[io+1], io, (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, 30.0, false );
            int i0=ff.getOrbOffset(io);
            Vec3d dir = ff.epos[i0+1] - ff.epos[i0]; dir.normalize();
            plotOrb     ( ff, plot1.lines[io+1], io, ff.epos[i0], dir, 20.0,  bPlotDens, true );
            //plotOrb     ( ff, plot1.lines[io+1], io, ff.epos[i0], dir, 100.0, bPlotDens );

        }
        //Draw3D::drawLine(  );
        //plotOrb     ( ff, plot1.lines[1], 0, (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, 100.0 );
        //plotOrb     ( ff, plot1.lines[2], 1, (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, 100.0 );
        //testDerivsTotal( 0, 0, 0, ff, plot1, 0 );
        plot1.bGrid=false;
        //plot1.bAxes=false;
        plot1.bTicks=false;
        plot1.update();
        plot1.render();
        //glCallList( ogl );
        //glDisable(GL_DEPTH_TEST);
        plot1.view();
}

void TestAppCLCFSF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    char strtmp[64];
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                case SDLK_i:
                    sprintf(strtmp, "temp/snapshot_iter_%i03.fgo", iter );
                    ff.saveToFile( strtmp );
                    break;
                case SDLK_p:     bDrawPlots   = !bDrawPlots;    break;
                case SDLK_o:     bDrawObjects = !bDrawObjects;  break;
                case SDLK_KP_1:  bDrawAtoms   = !bDrawAtoms;    break;
                case SDLK_KP_2:  bDrawWfs     = !bDrawWfs;      break;
                case SDLK_KP_3:  bDrawRho     = !bDrawRho;      break;
                case SDLK_r:     bPlotDens = !bPlotDens; break;
                //case SDLK_a:    = !perspective; break;
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


void TestAppCLCFSF::initff(){
    /*
    Vec3d pij; double sij;
    double S = Gauss::product3D_s_new( 1.0, {0.0,0.0,0.0}, 0.5, {2.0,0.0,0.0}, sij, pij );
    printf( "S %g pij.x %g sij %g \n", S, pij.x, sij );
    exit(0);
    */

    /*
    // ---- 2 atoms, 2 orbs, 2 basis per orb
    //      natom  nOrb perOrb natypes
    ff.realloc( 2, 2, 2, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
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
    ff.realloc( 1, 1, 1, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
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
    ff.realloc( 1, 1, 2, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
        _.setAtom    ( 0,   (Vec3d){   0.0, 0.0, 0.0 }, 1.0,  0.9, 1.0, 0. );
        _.setElectron( 0,0, (Vec3d){  +1.0, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 0,1, (Vec3d){  -1.0, 0.0, 0.0 }, 0.5, -1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- H2+ molecule ion - 2 basis electron in 2 atom potential   (2a,1o,2b)
    //      natom  nOrb perOrb natypes
    ff.realloc( 2, 1, 2, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
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
    ff.realloc( 0, 1, 2, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
        _.setElectron( 0,0, (Vec3d){  +0.7, 0.0, 0.0 }, 0.5, +1.0 );
        //_.setElectron( 0,1, (Vec3d){  -0.7, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 0,1, (Vec3d){  -0.7, 0.0, 0.0 }, 0.5, -1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- Free 2 electron 1 basis each    (0a,2o,1b)
    //      natom  nOrb perOrb natypes
    ff.realloc( 0, 2, 1, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
        _.setElectron( 0,0, (Vec3d){  +0.3, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 1,0, (Vec3d){  -0.3, 0.0, 0.0 }, 0.5, +1.0 );
    }
    dt = 0.001;
    */

    /*
    // ---- 2 electron 1 basis each in Soft Atom potential   (1a,2o,1b)
    //      natom  nOrb perOrb natypes
    ff.realloc( 1, 2, 1, 1 );
    ff.setRcut( 4.0 );
    ff.setDefaultValues();
    {CLCFGO& _=ff;
        _.setAtom    ( 0,   (Vec3d){  0.0, 0.0, 0.0 }, 1.0,  0.1, 1.0, 0. );
        _.setElectron( 0,0, (Vec3d){  +0.3, 0.0, 0.0 }, 0.5, +1.0 );
        _.setElectron( 1,0, (Vec3d){  -0.3, 0.0, 0.0 }, 0.5, +1.0 );
    }
    dt = 0.001;
    */
}

// ===================== MAIN

TestAppCLCFSF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppCLCFSF( junk , 1024, 800 );
	thisApp->loop( 1000000 );
	return 0;
}


