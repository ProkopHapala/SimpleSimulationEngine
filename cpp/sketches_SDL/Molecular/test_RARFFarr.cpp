
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
#include "VecN.h"

/*
#include "Multipoles.h"
#include "PotentialFlow.h"
#include "grids3D.h"
#include "MultipoleGrid.h"
*/

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

#include "Draw3D_Molecular.h"

//#include "MMFF.h"
//#include "RARFF.h"
//#include "RARFF2.h"
#include "RARFFarr.h"
#include "QEq.h"

#define R2SAFE  1.0e-8f


// ============= This is From LIB ReactiveFF.cpp
// /home/prokop/git/ProbeParticleModel_OclPy3/cpp/ReactiveFF.cpp


// ============= Global Variables

std::vector<RigidAtomType*> atomTypes;
RARFF2arr ff;
QEq       qeq;


// ============= Functions

template< typename FF>
int pickBond( FF ff, Vec3d& ray0, const Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int ia=0; ia<ff.natom; ia++){
        if(ff.ignoreAtoms)if(ff.ignoreAtoms[ia])continue;
        int nb = ff.types[ia]->nbond;
        for(int j=0; j<nb; j++){
            int i = ia*N_BOND_MAX + j;
            Vec3d pbond = ff.bondPos(i,1.0);
            double ti = raySphere( ray0, hRay, R, pbond );
            //printf( "atom %i t %g %g hRay(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, tmin, hRay.x,hRay.y,hRay.z,  ps[i].x, ps[i].y, ps[i].z );
            //printf( "atom %i t %g %g %g ray0(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, R, tmin, ray0.x,ray0.y,ray0.z,  ps[i].x, ps[i].y, ps[i].z );
            if(ti<tmin){ imin=i; tmin=ti; }
        }
    }
    return imin;
}

/*
int insertAtomType( int nbond, int ihyb, double rbond0, double aMorse, double bMorse, double c6, double R2vdW, double Epz ){
    RigidAtomType* atyp = new RigidAtomType();
    atomTypes.push_back( atyp );
    //RigidAtomType& atyp = atomTypes.back();
    atyp->nbond  =  nbond;    // number bonds
    atyp->rbond0 =  rbond0;
    atyp->aMorse =  aMorse;
    atyp->bMorse =  bMorse;
    atyp->Epz    =  Epz;
    atyp->c6     =  c6;
    atyp->R2vdW  =  R2vdW;
    printf("insertAtomType %i %i  %g %g %g %g %g ", nbond, ihyb, rbond0, aMorse, bMorse, c6, R2vdW );
    switch(ihyb){
        case 0: atyp->bh0s = (Vec3d*)sp1_hs; printf("sp1\n"); break;
        case 1: atyp->bh0s = (Vec3d*)sp2_hs; printf("sp2\n"); break;
        case 2: atyp->bh0s = (Vec3d*)sp3_hs; printf("sp3\n"); break;
    };
    return atomTypes.size()-1;
}

template<typename Func> void numDeriv( Vec3d p, double d, Vec3d& f, Func func){
    double d_=d*0.5;
    Vec3d p0=p;
    p.x=p0.x+d_; f.x = func(p); p.x-=d; f.x-=func(p); p.x+=d_;
    p.y=p0.y+d_; f.y = func(p); p.y-=d; f.y-=func(p); p.y+=d_;
    p.z=p0.z+d_; f.z = func(p); p.z-=d; f.z-=func(p); p.z+=d_;
    f.mul(1/d);
}
*/


class TestAppRARFF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    RigidAtomType type1,type2;

    int ipicked = -1;
    Vec3d ray0;
    Vec3d mouse_p0;

    int ogl_sph=0;

    Plot2D plot1;
    bool bRun = true;
    int    perFrame = 100;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    // for exp
    type1.nbond  = 3;  // number bonds
    type1.rbond0 =  0.6;
    type1.aMorse =  1.0;
    type1.bMorse = -0.7;
    type1.bh0s = (Vec3d*)sp2_hs;
    type1.print();
    //exit(0);

    type2.nbond  = 4;  // number bonds
    type2.rbond0 =  0.75;
    type2.aMorse =  1.0;
    type2.bMorse = -0.7;
    type2.bh0s = (Vec3d*)sp3_hs;
    type2.print();
    //exit(0);

    ff.bRepelCaps =  false; // ignore caps .... good e.g. for water


    /*
    srand(0);
    //srand(2);
    int nat = 12;
    ff.realloc(nat);
    for(int i=0; i<nat; i++){
        if(randf()>0.5){ ff.types[i]=&type1;  }else{ ff.types[i]=&type2; }
        ff.apos [i].fromRandomBox((Vec3d){-5.0,-5.0,-1.0},(Vec3d){5.0,5.0,1.0});
        ff.qrots[i].setRandomRotation();
    }
    ff.cleanAux();

    RigidAtomType pairType;
    pairType.combine(type1,type1);
    printf(" >>> pairType: <<<<\n");
    pairType.print();

    ff.projectBonds();
    ff.pairEF( 0, 1, 3,3, pairType );
    //exit(0);
    */

    int nat=1;
    ff.realloc( nat, 10 );
    for(int i=0; i<1; i++){
        //if(randf()>0.5){ ff.types[i]=&type1;  }else{ ff.types[i]=&type2; }
        ff.types[i]=&type2;
        ff.apos [i].fromRandomBox((Vec3d){-5.0,-5.0,-1.0},(Vec3d){5.0,5.0,1.0});
        ff.qrots[i].setRandomRotation();
    }
    ff.apos [0]=Vec3dZero;
    ff.qrots[0]=Quat4dIdentity;
    ff.cleanAux();

    //ff.tryResize( );

    /*
    ogl_sph = glGenLists(1);
    glNewList(ogl_sph, GL_COMPILE);
        //Draw3D::drawSphere_oct( 3, 1.0, (Vec3d){0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( 3, 0.2, (Vec3d){0.0,0.0,0.0} );
    glEndList();
    */

    Draw3D::makeSphereOgl( ogl_sph, 3, 0.25 );

}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable( GL_LIGHTING );

    if( ff.tryResize( ) );

    // ---------- Simulate
    //bRun = false;
    perFrame = 10;
    if(bRun){
        for(int i=0; i<perFrame; i++){
            ff.cleanAtomForce();
            ff.projectBonds();
            ff.interEF();
            if(ipicked>=0){
                Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 );
                ff.aforce[ipicked].add( f );
            }
            //for(int i=0; i<ff.natoms; i++){
            //    ff.aforce[i].add( getForceHamakerPlane( ff.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );
            //    //printf( "%g %g %g\n",  world.aforce[i].x, world.aforce[i].y, world.aforce[i].z );
            //}
            ff.applyForceHarmonic1D( Vec3dZ, 0.0, -10.1);
            ff.evalTorques();
            ff.moveMDdamp(0.05, 0.9);
        }
    }
    if( frameCount>10 ){
        //bRun=0;
        //exit(0);
    }

    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);

    // ---------- Draw
    glColor3f(0.0,0.0,0.0);
    double fsc = 0.1;
    double tsc = 0.1;
    for(int ia=0; ia<ff.natom; ia++){
        if(ff.ignoreAtoms[ia])continue;
        glColor3f(0.3,0.3,0.3);
        Draw3D::drawShape( ogl_sph, ff.apos[ia], Mat3dIdentity );
        for(int j=0; j<ff.types[ia]->nbond; j++){
            int i=ia*N_BOND_MAX+j;
            Vec3d pb = ff.bondPos( i );
            //printf( "bondCaps[%i] %i\n", i, ff.bondCaps[i] );
            if( ff.bondCaps[i]>=0 ){ glColor3f(1.0,0.0,0.0); } else{ glColor3f(0.0,0.0,0.0); }
            Draw3D::drawLine( ff.apos[ia] , pb );
            glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( ff.fbonds[i]*fsc, pb );
            //glColor3f(0.0,0.0,0.0); Draw3D::drawVecInPos( ff.hbonds[i], ff.apos[i] );
            //glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( ff.fbonds[io]*fsc, ff.apos[i]+ff.hbonds[io] );
        }

    };
    Draw3D::drawAxis( 1.0);
};

void TestAppRARFF::drawHUD(){
	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	plot1.view();
}

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
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
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    mouse_p0 = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                    printf( "picked atom %i \n", ipicked );
                    break;
                case SDL_BUTTON_RIGHT:
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    printf( "remove atom %i \n", ipicked );
                    ff.ignoreAtoms[ ipicked ] = true;
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if( ipicked==-1 ){
                        Vec3d mouse_p = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                        Vec3d dir = mouse_p-mouse_p0;
                        if(dir.norm2()<1e-6){ dir=(Vec3d)cam.rot.b; };
                        int ib = pickBond( ff, mouse_p0, (Vec3d)cam.rot.c, 0.5 );
                        if(ib>0){
                            printf( "add atom to bond %i of atom %i \n", ib%N_BOND_MAX, ib/N_BOND_MAX );
                            Vec3d p0 = ff.bondPos(ib, 2.0 );
                            //ff.inserAtom( {nbBrush,4,4}, mouse_p0, dir, (Vec3d)cam.rot.b );
                            //ff.inserAtom( &type1, (const int[]){0,0,0,0}, p0, dir, (Vec3d)cam.rot.b  );
                            int ia = ff.inserAtom( &type2, (const int[]){1,1,-1,-1}, p0, ff.hbonds[ib]*-1, dir  );
                            ff.projectAtomBons(ia);
                            bRun = 0;
                        }
                    }
                    ipicked = -1;
                    break;
                case SDL_BUTTON_RIGHT:
                    ipicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
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



