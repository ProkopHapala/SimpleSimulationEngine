
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

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

#include "Draw3D_Molecular.h"

//#include "MMFF.h"

//#include "RARFFarr.h"


int i_DEBUG = 0;
#include "DynamicOpt.h"
#include "FlexibleAtomReactiveFF.h"
#include "FlexibleAtomReactiveFF_dyn.h"

// ======= THE CLASS

class TestAppFARFF: public AppSDL2OGL_3D { public:

    bool bRun = false;
    int perFrame = 1;

    FARFF ff;
    DynamicOpt opt;

    Plot2D plot1;

    int ogl_sph=0;

    double  Emin,Emax;
    int     npoints;
    Vec3d*  points  =0;
    double* Energies=0;
    Vec3d * Forces  =0;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppFARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppFARFF::TestAppFARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //srand(15480);  int nat = 15; double sz = 4; double szH = 1;   // Test 1

    srand(1581);  int nat = 40; double sz = 5; double szH = 0.1;    // Test 2

    ff.realloc(nat);
    for(int ia=0; ia<nat; ia++){
        ff.apos[ia].fromRandomBox( (Vec3d){-sz,-sz,-szH},(Vec3d){sz,sz,szH} );
        ff.aconf[ia].set(4,4,4);
        double rnd=randf(); if(rnd>0.7){ if(rnd<0.9){ ff.aconf[ia].a=3; }else{ ff.aconf[ia].a=2; }  }
        for(int j=0; j<N_BOND_MAX; j++){
            int io = j+ia*N_BOND_MAX;
            //printf( "atom[%i] %i %i \n", ia, j, io );
            ff.opos[io].fromRandomSphereSample();
        }
    }


    //ff.atype0.Kee = -0.5;
    //ff.atype0.Kpp = 0;

    /*
    int ia,io;
    ia=0; io=ia*N_BOND_MAX;
    ff.apos [ia].set( 0.0,0.0,0.0);
    ff.aconf[ia].set(3,4,4);
    ff.opos [io+0].set(-1.0,0.0,-0.1);
    ff.opos [io+1].set(+1.0,0.0,-0.1);
    ff.opos [io+2].set(0.0,+1.0, 1.3);
    ff.opos [io+3].set(0.0,-1.0,-1.3);

    ia=1; io=ia*N_BOND_MAX;
    ff.apos [ia].set(-1.0,0.0,0.0);
    ff.aconf[ia].set(3,4,4);
    ff.opos [io+0].set(+1.0,0.0,-0.1);
    ff.opos [io+1].set(-1.0,0.0,-0.1);
    ff.opos [io+2].set(0.0,+1.0,+0.1);
    ff.opos [io+3].set(0.0,-1.0,+0.1);

    ia=2; io=ia*N_BOND_MAX;
    ff.apos [ia].set(+1.2,0.0,0.0);
    ff.aconf[ia].set(4,4,4);
    ff.opos [io+0].set(-1.0,0.0,+0.1);
    ff.opos [io+1].set( 1.0,0.0,+0.1);
    ff.opos [io+2].set(0.0,+1.0,-0.1);
    ff.opos [io+3].set(0.0,-1.0,-0.1);
    */


    // =========== check numerical derivatives

    /*
    auto check_froce_atom0 = [&](const Vec3d& p,Vec3d& f)->double{
        ff.apos[0] = p;
        ff.cleanForce();
        double E = ff.eval();
        f = ff.aforce[0] * -1;
        return E;
    };

    auto check_froce_orb0 = [&](Vec3d h,Vec3d& f)->double{
        //printf( "check_froce_orb0 h (%g,%g,%g)\n", h.x, h.y, h.z );
        ff.opos[0] = h;
        ff.cleanForce();
        double E = ff.eval();
        ff.oforce[0].makeOrthoU(ff.opos[0]);
        f = ff.oforce[0] * -1;
        return E;
    };

    Vec3d fE,f, p = {-1.0,0.2,0.3}; p.normalize();
    printf("Atom[0] "); checkDeriv(check_froce_atom0,{1.0,1.0,1.0}, 0.0001,fE, f );
    printf("Orb [0] "); checkDeriv(check_froce_orb0 ,p, 0.0001,fE, f );
    //exit(0);
    */


    //ff.cleanForce();
    //ff.eval();


    opt.bindOrAlloc( 3*ff.nDOF, (double*)ff.dofs, 0, (double*)ff.fdofs, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    //opt.initOpt( 0.05, 0.2 );
    opt.initOpt( 0.01, 0.05 );

    //exit(0);

    Draw3D::makeSphereOgl( ogl_sph, 3, 0.25 );

}

void TestAppFARFF::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    //if(bRun){
    perFrame=10;
    for(int itr=0;itr<perFrame;itr++){
        printf( " ==== frame %i i_DEBUG  %i \n", frameCount, i_DEBUG );
        double F2 = 1.0;

        ff.cleanForce();
        ff.eval();
        //ff.apos[0].set(.0);

        //if(bRun)ff.moveGD(0.001, 1, 1);

        F2 = opt.move_FIRE();

        //for(int i=0; )printf( "",  )

        printf( " |F| %g \n", sqrt(F2) );
        //if(!(F2<1000000.0))perFrame=0;

        /*
        double cosdRot,cosdPos;
        ff.getCos( cosdRot, cosdPos );
        double damp = fmax( 0.9, fmin(cosdRot,cosdPos) );
        printf( "cosdRot %g cosdRot %g damp %g \n", cosdRot, cosdPos, damp );
        ff.moveMDdamp(0.1, damp );
        */

    }
    //}

    //Vec3d bhs[N_BOND_MAX];
    //atom1.torq = (Vec3d){0.1,0.0,0.0};
    //atom1.moveRotGD(0.8);
    //printf( "qrot (%g,%g,%g,%g)\n", atom1.qrot.x, atom1.qrot.y, atom1.qrot.z, atom1.qrot.w );

    glColor3f(1.0,1.0,1.0);
    //drawRigidAtom( atom1 );

    double fsc = 0.1;
    double tsc = 0.1;
    for(int i=0; i<ff.natom; i++){
        glColor3f(0.0,0.0,0.0); Draw3D::drawPointCross(ff.apos[i], 0.1 );
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( ff.aforce[i]*fsc, ff.apos[i]  );
        for(int j=0;j<N_BOND_MAX; j++){
            int io = j + i*N_BOND_MAX;
            if( j<ff.aconf[i].a ){
                glColor3f(0.0,0.0,0.0); Draw3D::drawVecInPos( ff.opos[io],   ff.apos[i] );
                //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( ff.oforce[io], ff.apos[i]+ff.opos[io] );
            }else{
                glColor3f(0.0,0.5,0.0); Draw3D::drawLine    ( ff.apos[i]+ff.opos[io]*0.5, ff.apos[i]-ff.opos[io]*0.5 );
            }
        }
        //glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( ff.atoms[i].torq*tsc,  ff.atoms[i].pos  );
    };

    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.25, 1.0 );
    //Draw3D::shapeInPoss( ogl_sph, ff.natom, ff.apos, ff.aenergy );
    glColor3f(1.0,1.0,1.0);
    Draw3D::shapeInPoss( ogl_sph, ff.natom, ff.apos, 0 );


/*
    printf("npoints %i Emin %g Emax %g \n",npoints, Emin, Emax);
    glPointSize(5);
    drawScalarArray( npoints, points, Energies, Emin, Emax );
*/
/*
    glColor3f(0.0,1.0,0.0);
    drawVectorArray( npoints, points, Forces, 0.02 );
*/

    //Draw3D::drawAxis( 1.0);

    //exit(0);
};


void TestAppFARFF::drawHUD(){
/*
    glColor3f(1.0,1.0,1.0);
    txt.viewHUD( {100,220}, fontTex );

	gui.draw();

	glTranslatef( 10.0,300.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, fontSizeDef, {10,5} );
*/

	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	plot1.view();

}


void TestAppFARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;


                //case SDLK_LEFTBRACKET:  i_DEBUG=(i_DEBUG+1)%6; break;
                case SDLK_RIGHTBRACKET: i_DEBUG=(i_DEBUG+1)%6; printf("i_DEBUG %i\n", i_DEBUG); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppFARFF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppFARFF( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















