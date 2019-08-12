
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

//#include "MMFF.h"

//#include "RARFF.h"
//#include "RspFF.h"
//#include "FspFF.h"
#include "FspFFclean.h"

#include "DynamicOpt.h"

#define R2SAFE  1.0e-8f

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    NBFF  nff; // non-bonded forcefield
    FspFF ff;  // flexible atom sp-bonding forcefield
    DynamicOpt opt;

    int      fontTex;

    int ipick = 0;
    bool bRun = 0;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //https://stackoverflow.com/questions/2279180/does-c-have-with-keyword-like-pascal
    {auto &_ = ff;  // instead of WIDTH

        /*
        _.realloc( 2, 1 );
        _.apos[0].set(0.,0.5,0.);
        _.apos[1].set(0.,1. ,0.);
        _.bLKs[0].set(1.5 ,1.0);
        */

        /*
        _.realloc( 3, 2 );
        _.apos[0].set(1.,-1.,0.);
        _.apos[1].set(0.,0.,0.);
        _.apos[2].set(1.,1.,0.);
        _.bLKs[0].set(1.0,1.0);
        _.bLKs[1].set(1.5,1.0);
        */

        /*
        _.realloc( 4, 3 );
        _.apos[0].set(0.,0.,0.);
        _.apos[1].set(0.,1.,0.);
        _.apos[2].set(1.,1.,0.);
        _.apos[3].set(1.,0.,0.);
        */

        _.realloc( 6, 6 );
        _.apos[0].set(0.,0.,0.);
        _.apos[1].set(0.,1.,0.);
        _.apos[2].set(1.,1.,0.);
        _.apos[3].set(1.,0.5,1.);
        _.apos[4].set(1.,0.,0.);
        _.apos[5].set(-1.,0.,0.);

        Vec2i b2a[] = { {0,1}, {1,2}, {2,3}, {3,4}, {4,0}, {5,0} };
        Vec2i hps[] = { {0,0}, {1,1}, {1,1}, {0,0}, {2,0}, {3,0} };
        //Vec2i hps[] = { {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0} };
        //Vec2i hps[] = { {0,1}, {0,1}, {0,1}, {0,1}, {0,1}, {0,0} };

        _.setBondsAndHydrogens( b2a, hps );
        printf( "natom %i nbonds %i ncap %i nepair %i nporb %i \n", _.natom, _.nbond, _.ncap, _.nepair, _.nporb );
        _.guessOrbs();
    }

    nff.realloc( ff.natom + ff.ncap, ff.nbond+ff.ncap );
    nff.setREQs(0       ,ff.natom,{1.9080,sqrt(0.003729),0});
    nff.setREQs(ff.natom,nff.n   ,{1.4870,sqrt(0.000681),0});
    int npair = ff.outputParticleBonds( nff.pairMask );
    printf("found %i bond pairs\n", npair);
    for(int i=0; i<npair; i++){ printf("pair[%i]: (%i,%i) \n", i, nff.pairMask[i].a, nff.pairMask[i].b ); };

    opt.bindOrAlloc(ff.nDOF*3,(double*)ff.dofs,0,(double*)ff.fdofs,0);
    //opt.initOpt(0.5,0.1);
    opt.initOpt(0.1,0.1);

    /*
    ff.apos[0]=ff.apos[1];
    ff.apos[0].add(0.2,0.2,0.2);
    ff.hdir[3].mul(-1);
    */

    //exit(0);
}

void TestAppRARFF::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	/*
	Vec3d cog  = Vec3dZero; for(int i=0;i<ff.natom;i++){ cog .add(ff.apos  [i]); }; cog.mul(1./ff.natom);
	//Vec3d fsum = Vec3dZero; for(int i=0;i<ff.natom;i++){ fsum.add(ff.aforce[i]); };
	Vec3d fsum = Vec3dZero; for(int i=0;i<ff.nDOF;i++){ fsum.add(ff.fdofs[i]); };
	*/
	/*
	Vec3d cog,fsum,tqsum;
	for(int i=0;i<ff.natom;i++){ cog .add(ff.apos  [i]); }; cog.mul(1./ff.natom);
	Vec3d fsum = Vec3dZero; for(int i=0;i<ff.nDOF;i++){
        fsum.add(ff.fdofs[i]);
        tq.sum()
	};
	*/

    //bRun=(frameCount==0);
	if(bRun){
        int perFrame=1;
        //perFrame=1;
        for(int i=0; i<perFrame; i++){
            ff.cleanForce();

            ff.evalForces();

            //ff.projectBondCenter();
            //ff.evalAtoms();
            //ff.transferBondForce();
            //for(int i=0;i<ff.natom;i++){ ff.aforce[i].set(0.); };

            int np1 = ff.outputParticlePositions(nff.ps);
            nff.cleanForce();
            nff.evalLJQ_sortedMask();
            int np2 = ff.inputParticleForces(nff.fs);
            //printf( "%i %i %i \n", nff.n, np1, np2 );

            // --- check force
            Vec3d cog,fsum,tqsum;
            ff.checkForceTorque(cog,fsum,tqsum);
            //printf( "===== frame %i \n", frameCount );
            printf("fsum %g tqsum %g | cog(%g,%g,%g)\n", fsum.norm(), tqsum.norm(), cog.x,cog.y,cog.z );
            glColor3f(1,1,1); Draw3D::drawPointCross(cog,0.05);
            glColor3f(1,1,1); Draw3D::drawVecInPos(fsum,cog);

            double f2err=0;
            for(int i=0; i<ff.nDOF; i++){ f2err=fmax(f2err,ff.fdofs[i].norm2()); }
            printf( "|F| %g \n", sqrt(f2err) );

            //ff.moveGD( 0.1  /fmax(1,sqrt(f2err)) );
            opt.move_FIRE();

        }

        ff.projectBondCenter();

        for(int ia=0; ia<ff.natom; ia++){
            for(int io=0; io<ff.aconf[ia].b; io++ ){
                //Draw3D::drawVecInPos( ff.hforce[ia*N_BOND_MAX+io]*fsc, ff.hdir[ia*N_BOND_MAX+io] );
                char c='H';
                if(io< ff.aconf[ia].c)c='b';
                if(io>=ff.aconf[ia].a)c='e';
                //printf( "%i[%i] %c: l=%g\n", ia, io, c, (ff.hdir[ia*N_BOND_MAX+io]-ff.apos[ia]).norm()  );
            }
        }

    }

    glColor3f(1.,1.,1.);
    for(int ia=0; ia<nff.n; ia++){
        Draw3D::drawPointCross( nff.ps[ia], 0.1 );
    }
    glColor3f(0.,0.,0.);
    for(int ib=0; ib<nff.nmask; ib++){
        //Draw3D::drawLine( nff.ps[nff.pairMask[ib].a], nff.ps[nff.pairMask[ib].b] );
    }
    glColor3f(1.,0.,0.);
    for(int ia=0; ia<nff.n; ia++){
        //Draw3D::drawVecInPos( nff.fs[ia], nff.ps[ia] );
    }

    //return;

    double fsc = 10.0;

    // ================ view froces

    // atom forces
    glColor3f(1.,0.,0.);
    for(int ia=0; ia<ff.natom; ia++){
        Draw3D::drawVecInPos( ff.aforce[ia]*fsc, ff.apos[ia] );
    }

    // sigma force (bond,epair,hydrogen)
    glColor3f(1.,1.,0.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=0; io<ff.aconf[ia].b; io++ ){
            Draw3D::drawVecInPos( ff.hforce[ia*N_BOND_MAX+io]*fsc, ff.hdir[ia*N_BOND_MAX+io] );
        }
    }

    // ================ view positions

    // atom centers
    glColor3f(0.,0.,0.);
    for(int i=0; i<ff.natom; i++){
        //Draw3D::drawPointCross( ff.apos[i], 0.1 );
    }

    // bonds
    glColor3f(0.,0.,0.);
    for(int i=0; i<ff.nbond; i++){
        //Draw3D::drawLine( ff.apos[ ff.bondIs[i].a>>2 ], ff.apos[ ff.bondIs[i].b>>2 ] );
    }

    // bond electrons
    glColor3f(0.5,0.,0.5);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=0; io<ff.aconf[ia].c; io++ ){
            Draw3D::drawLine( ff.hdir[ia*N_BOND_MAX+io], ff.apos[ia] );
        }
    }

    //return;

    // cap hydrogens
    glColor3f(0.,0.,0.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].c; io<ff.aconf[ia].a; io++ ){
            Draw3D::drawLine( ff.hdir[ia*N_BOND_MAX+io],  ff.apos[ia] );
        }
    }

    // e-pair
    glColor3f(0.,0.7,1.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].a; io<ff.aconf[ia].b; io++ ){
            Draw3D::drawLine(  ff.hdir[ia*N_BOND_MAX+io], ff.apos[ia] );
        }
    }

    // pi-bond
    glColor3f(0.4,.8,0.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].b; io<N_BOND_MAX; io++ ){
            Draw3D::drawVecInPos( ff.hdir[ia*N_BOND_MAX+io]*0.5, ff.apos[ia] );
        }
    }

    //Draw3D::drawAxis( 1.0);

};


void TestAppRARFF::drawHUD(){
/*
    glColor3f(1.0,1.0,1.0);
    txt.viewHUD( {100,220}, fontTex );

	gui.draw();

	glTranslatef( 10.0,300.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    Draw::drawText( "abcdefghijklmnopqrstuvwxyz \n0123456789 \nABCDEFGHIJKLMNOPQRTSTUVWXYZ \nxvfgfgdfgdfgdfgdfgdfg", fontTex, fontSizeDef, {10,5} );
*/

/*
	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 10.0,10.00,1.0  );
	plot1.view();
	*/

}

void TestAppRARFF::keyStateHandling( const Uint8 *keys ){

    double dstep = 0.1;
	if( keys[ SDL_SCANCODE_A ] ){ ff.apos[ipick].x += 0.1; }
	if( keys[ SDL_SCANCODE_D ] ){ ff.apos[ipick].x -= 0.1; }
    if( keys[ SDL_SCANCODE_W ] ){ ff.apos[ipick].y += 0.1; }
	if( keys[ SDL_SCANCODE_S ] ){ ff.apos[ipick].y -= 0.1; }
    if( keys[ SDL_SCANCODE_Q ] ){ ff.apos[ipick].z += 0.1; }
	if( keys[ SDL_SCANCODE_E ] ){ ff.apos[ipick].z -= 0.1; }


/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }

	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b,  cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_Q ] ){ cam.pos.add_mul( cam.rot.c, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_E ] ){ cam.pos.add_mul( cam.rot.c,  cameraMoveSpeed ); }
*/
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
















