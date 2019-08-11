
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
#include "FspFF.h"

#define R2SAFE  1.0e-8f

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:

    FspFF ff;

    int      fontTex;

    int ipick = 0;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    ff.realloc( 4, 3, 0 );

    ff.apos[0].set(0.,0.,0.);
    ff.apos[1].set(0.,1.,0.);
    ff.apos[2].set(1.,1.,0.);
    ff.apos[3].set(1.,0.,0.);

    //ff.bondIs[0].set((0<<2)+1,1<<2);
    //ff.bondIs[1].set((1<<2)+1,2<<2);
    //ff.bondIs[2].set((2<<2)+1,3<<2);
    //ff.bondIs[3].set((3<<2)+1,0<<2);

    ff.bLKs[0].set(1.5,1.0);
    ff.bLKs[1].set(1.0,1.0);
    ff.bLKs[2].set(1.0,1.0);
    //ff.bLKs[3].set(1.5,1.0);

    Vec2i b2a[] = { {0,1}, {1,2}, {2,3}, {3,0} };
    Vec2i hps[] = { {0,0}, {0,0}, {0,0}, {0,0} };

    ff.setBondsAndHydrogens( b2a, hps );
    ff.guessOrbs();

}

void TestAppRARFF::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	printf( "===== frame %i \n", frameCount );
	ff.evalForces();
	ff.moveGD(0.1);

	ff.projectBondCenter();

    // atom centers
    glColor3f(0.,0.,0.);
    for(int i=0; i<ff.natom; i++){
        Draw3D::drawPointCross( ff.apos[i], 0.1 );
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
    glColor3f(0.,.7,0.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].b; io<N_BOND_MAX; io++ ){
            Draw3D::drawVecInPos( ff.hdir[ia*N_BOND_MAX+io]*0.5, ff.apos[ia] );
        }
    }


    /*
    // bond dirs
    glColor3f(0.5,0.,0.5);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=0; io<ff.aconf[ia].c; io++ ){
            Draw3D::drawVecInPos( ff.hdir[ia*N_BOND_MAX+io]*0.3, ff.apos[ia] );
        }
    }

    //return;

    // cap hydrogens
    glColor3f(0.,0.,0.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].c; io<ff.aconf[ia].a; io++ ){
            Draw3D::drawVecInPos( ff.hdir[ia*N_BOND_MAX+io],  ff.apos[ia] );
        }
    }

    // e-pair
    glColor3f(0.,0.7,1.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].a; io<ff.aconf[ia].b; io++ ){
            Draw3D::drawVecInPos(  ff.hdir[ia*N_BOND_MAX+io]*0.5, ff.apos[ia] );
        }
    }

    // pi-bond
    glColor3f(0.,.7,0.);
    for(int ia=0; ia<ff.natom; ia++){
        for(int io=ff.aconf[ia].b; io<N_BOND_MAX; io++ ){
            Draw3D::drawVecInPos( ff.hdir[ia*N_BOND_MAX+io]*0.5, ff.apos[ia] );
        }
    }
    */
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
//                case SDLK_SPACE: bRun = !bRun;
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
















