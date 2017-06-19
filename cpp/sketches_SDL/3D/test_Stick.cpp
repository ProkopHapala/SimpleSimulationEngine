
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
#include "MMFF.h"
#include "DynamicOpt.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ==========================
// TestAppStick
// ==========================

double R_target = 1.0;
Vec2d p_target = (Vec2d){6.0,0.0};


inline Vec3d hitForce( const Vec3d& pos, const Vec3d& vel ){
    Vec2d dp; dp.set_sub(pos.xy(), p_target);
    double r2 = dp.norm2();
    if(r2<(R_target*R_target)){
        return (Vec3d){-dp.x, -dp.y, 0.0};
    }
    return (Vec3d){0.0,0.0,0.0};
}

class TestAppStick : public AppSDL2OGL_3D {
	public:
	MMFFparams  params;
    MMFF        world;
    DynamicOpt  opt;

    int     fontTex;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1;
    int perFrame =  50;

    double drndv =  10.0;
    double drndp =  0.5;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppStick( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppStick::TestAppStick( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    params.loadBondTypes("common_resources/BondTypes.dat");

    world.allocate( 10, 9, 8, 0 );


    //for(int i=0; i<world.natoms; i++){
    //    world.apos[i].set( i, i, 0.0  );
    //}

    for(int i=0; i<world.nbonds; i++){
        world.bond2atom[i] = {i,i+1};
        world.bond_0[i] = 1.0;
        world.bond_k[i] = 2.0;
    }

    printf( "DEBUG 3 \n" );

    for(int i=0; i<world.nang; i++){
        world.ang2bond[i] = {i,i+1};
        world.ang_0   [i] = {1.0,0.0};
        world.ang_k   [i] = 0.5;
        //Vec2i ib = world.ang2bond[i];
        //world.ang2atom [i] = (Vec3i){ world.bond2atom[ib.x].y, world.bond2atom[ib.y].y, world.bond2atom[ib.y].x };
    }

    world.ang_b2a();

    //exit(0);

    world.printBondParams();
    //exit(0);

    opt.bindArrays( 3*world.natoms, (double*)world.apos, new double[3*world.natoms], (double*)world.aforce );
    opt.setInvMass( 1.0 );
    opt.cleanVel( );
    opt.initOpt( 0.01, 0.9 );

    for(int i=0; i<world.natoms; i++){
        ((Vec3d*)opt.pos)[i].set( i,  i, 0.0  );
        ((Vec3d*)opt.vel)[i].set( 0.0  );
        //((Vec3d*)opt.vel)[i].set( i*0.001, -i*0.001, 0.0  );
    }


	for(int itr=0; itr<1000; itr++){
        for(int i=0; i<world.natoms; i++){ world.aforce[i].set(0.0d); }
        world.eval_bonds(false);
        world.eval_angcos();
        //for(int i=0; i<opt.n; i++){ opt.vel[i] *=0.9999; }
        //opt.move_LeapFrog(0.01);
        opt.move_FIRE();
        world.apos[0].set(0.0);
    }

    for(int i=0; i<world.natoms; i++){
        //((Vec3d*)opt.pos)[i].set( i,  i, 0.0  );
        ((Vec3d*)opt.vel)[i].set( i*0.003, -i*0.003, 0.0  );
    }

}

void TestAppStick::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    ray0 = camMat.a*mouse_begin_x + camMat.b*mouse_begin_y;
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( world.apos[ipicked], ray0);

	double F2;
	for(int itr=0; itr<perFrame; itr++){

        for(int i=0; i<world.natoms; i++){
            //world.aforce[i].set(0.0d);
            world.aforce[i] = hitForce( world.apos[i], ((Vec3d*)opt.vel)[i] ) * -0.1;
        }

        world.eval_bonds(false);
        world.eval_angcos();
        //for(int i=0; i<opt.n; i++){ opt.vel[i] *=0.9999; }
        opt.move_LeapFrog(0.01);
        // F2 = opt.move_FIRE();
        world.apos[0].set(0.0);
    }

    for(int i=0; i<world.nbonds; i++){
        Vec2i ib = world.bond2atom[i];
        glColor3f(0.0f,0.0f,0.0f);
        if(i==ibpicked) glColor3f(1.0f,0.0f,0.0f); ;
        Draw3D::drawLine(world.apos[ib.x],world.apos[ib.y]);
        sprintf(str,"%i\0",i);
        Draw3D::drawText(str, (world.apos[ib.x]+world.apos[ib.y])*0.5, fontTex, 0.02, 0,0);
    }

    Draw2D::drawCircle_d( p_target, R_target, 64, false );

    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    for(int i=0; i<world.natoms; i++){
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross(world.apos[i],0.2);
        glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(world.aforce[i]*30.0,world.apos[i]);

    }
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);

};


void TestAppStick::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                case SDLK_w: world.apos[1].mul( 1.1 ); break;
                case SDLK_s: world.apos[1].mul( 0.9 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = pickParticle( world.natoms, world.apos, ray0, camMat.c , 0.5 );
                    break;
                case SDL_BUTTON_RIGHT:
                    ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
                    printf("ibpicked %i \n", ibpicked);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = -1;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppStick::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppStick * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppStick( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















