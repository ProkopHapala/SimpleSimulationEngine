
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
    //  type1,type2, typH2O;
    RigidAtomType  typeList [4];
    //RigidAtomType* typeList_[3];
    RigidAtomType* curType=0;

    CapType capTypes[2];

    Quat4i capsBrush;

    bool bBlockAddAtom=false;
    int ipicked = -1;
    Vec3d ray0;
    Vec3d mouse_p0;

    int ogl_sph=0;

    const char* workFileName="data/work.rff";

    Plot2D plot1;
    bool bRun = true;
    int    perFrame = 100;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    {RigidAtomType& typ=typeList[0]; // C-sp2
    //typeList_[0]=&typ;
    typ.id = 0;
    typ.nbond  = 3;  // number bonds
    typ.rbond0 =  1.2/2;
    typ.aMorse =  1.0;
    typ.bMorse = -0.7;
    typ.bh0s = (Vec3d*)sp2_hs;
    typ.print();
    //exit(0);
    }
    {RigidAtomType& typ=typeList[1]; // C-sp3
    //typeList_[1]=&typ;
    typ.id = 1;
    typ.nbond  = 4;  // number bonds
    typ.rbond0 =  1.5/2;
    typ.aMorse =  1.0;
    typ.bMorse = -0.7;
    typ.bh0s = (Vec3d*)sp3_hs;
    typ.print();
    }
    {RigidAtomType& typ=typeList[2];  // H2O
    //typeList_[2]=&typ;
    //typ.name = "O";
    strcpy(typ.name, "O");
    typ.id = 2;
    typ.nbond  = 4;
    typ.rbond0 = 2.9/2;
    typ.aMorse = 1.0;
    typ.bMorse = -0.7;
    typ.bh0s   = (Vec3d*)sp3_hs;
    typ.print();
    }
    {RigidAtomType& typ=typeList[3];  // H2O
    //typeList_[2]=&typ;
    //typ.name = "O";
    strcpy(typ.name, "O");
    typ.id = 3;
    typ.nbond  = 3;
    typ.rbond0 = 2.9/2;
    typ.aMorse = 1.0;
    typ.bMorse = -0.7;
    typ.bh0s   = (Vec3d*)sp2_hs;
    typ.print();
    }
    curType = &typeList[3];

    {CapType& typ=capTypes[1]; // C-sp3
        typ.id = 1;
        typ.rbond0 = 1.0;
        strcpy(typ.name, "H");
    }

    capsBrush.set(1,1,-1,-1);

    ff.capTypeList = capTypes;
    ff.typeList    = typeList;

    ff.bRepelCaps =  false; // ignore caps .... good e.g. for water
    ff.bDonorAcceptorCap = true;


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
    ff.realloc( nat, 100 );
    for(int i=0; i<1; i++){
        //if(randf()>0.5){ ff.types[i]=&type1;  }else{ ff.types[i]=&type2; }
        ff.types[i]=curType;
        ff.apos [i].fromRandomBox((Vec3d){-5.0,-5.0,-1.0},(Vec3d){5.0,5.0,1.0});
        ff.qrots[i].setRandomRotation();
    }
    ff.apos [0]=Vec3dZero;
    ff.qrots[0]=Quat4dIdentity;
    ff.cleanAux();
    ((Quat4i*)ff.bondCaps)[0]=capsBrush;

    //ff.resize(15);

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


    //if( ff.tryResize( 5, 100, 10) );

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
            //for(int i=0; i<ff.natom; i++){ ff.aforce[i].add( getForceHamakerPlane( ff.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );}
            //ff.applyForceHarmonic1D( Vec3dZ, 0.0, -30.1);
            ff.applyForceHarmonic1D( Vec3dZ, 0.0, -1.0);
            //ff.applyForceHamacker(true, -5.0, 10.0, 2.0 );
            ff.evalTorques();
            ff.moveMDdamp(0.05, 0.9);
        }
    }else{
        if(ipicked>=0){
            Vec3d dpos = ray0 - ff.apos[ipicked];
            dpos.makeOrthoU( (Vec3d)cam.rot.c );
            ff.apos[ipicked].add(dpos);
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
    //printf( "ff.natom %i \n", ff.natom );
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


void TestAppRARFF::keyStateHandling( const Uint8 *keys ){

    if( keys[ SDL_SCANCODE_R ] ){
        if(ipicked>=0){
            //Mat3d rot;
            Quat4d qrot;  qrot.f=(Vec3d)cam.rot.c*0.01; qrot.normalizeW();
            ff.qrots[ipicked].qmul(qrot);
            ff.projectAtomBons(ipicked);
        }
    }

	AppSDL2OGL_3D::keyStateHandling( keys );
};

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );

    Quat4i* caps = &capsBrush;
    if(ipicked>=0) caps=((Quat4i*)ff.bondCaps)+ipicked;

    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;

                case SDLK_KP_0: caps->array[0]*=-1; break;
                case SDLK_KP_1: caps->array[1]*=-1; break;
                case SDLK_KP_2: caps->array[2]*=-1; break;
                case SDLK_KP_3: caps->array[3]*=-1; break;
                case SDLK_KP_MULTIPLY: caps->mul(-1); printf("capsBrush(%i,%i,%i,%i)\n",  caps->x,caps->y,caps->z,caps->w ); break;


                case SDLK_l: ff.load(workFileName); bRun=false; break;
                case SDLK_k: ff.save(workFileName);  break;
                case SDLK_x: ff.saveXYZ("data/RARFF_out.xyz");  break;
                case SDLK_SPACE: bRun = !bRun;  break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:{
                    int ip = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    bBlockAddAtom=false;
                    if( ip>=0){
                        if(ipicked==ip){ ipicked=-1; bBlockAddAtom=true; printf("inv\n"); }
                        else           { ipicked=ip;                     printf("set\n"); }
                    }else{ ipicked=-1; }
                    mouse_p0 = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                    printf( "LMB DOWN picked %i/%i bblock %i \n", ipicked,ip, bBlockAddAtom );
                    }break;
                case SDL_BUTTON_RIGHT:{
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    printf( "remove atom %i \n", ipicked );
                    ff.ignoreAtoms[ ipicked ] = true;
                    }break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    printf( "LMB UP picked %i bblock %i \n", ipicked, bBlockAddAtom );
                    if( (ipicked==-1)&&(!bBlockAddAtom) ){
                        Vec3d mouse_p = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                        Vec3d dir = mouse_p-mouse_p0;
                        if(dir.norm2()<1e-6){ dir=(Vec3d)cam.rot.b; };
                        int ib = pickBond( ff, mouse_p0, (Vec3d)cam.rot.c, 0.5 );
                        if(ib>=0){
                            printf( "add atom to bond %i of atom %i \n", ib%N_BOND_MAX, ib/N_BOND_MAX );
                            Vec3d p0 = ff.bondPos(ib, 2.0 );
                            //ff.inserAtom( {nbBrush,4,4}, mouse_p0, dir, (Vec3d)cam.rot.b );
                            //ff.inserAtom( &type1, (const int[]){0,0,0,0}, p0, dir, (Vec3d)cam.rot.b  );
                            int ia = ff.inserAtom( curType, (int*)&capsBrush, p0, ff.hbonds[ib]*-1, dir  );
                            ff.projectAtomBons(ia);
                            bRun = 0;
                        }
                        //ipicked = -1;
                    }
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



