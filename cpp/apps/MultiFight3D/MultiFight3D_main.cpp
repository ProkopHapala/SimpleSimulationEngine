
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "MultiFight3DWorld.h"
#include "Warrior3D.h"
#include "Projectile3D.h"

class MultiFight3D_single : public AppSDL2OGL_3D {
	public:
    MultiFight3DWorld world;
    double dvel = 10.0;

    Warrior3D *warrior1;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );
    virtual void mouseHandling( );

	MultiFight3D_single( int& id, int WIDTH_, int HEIGHT_ );

};

MultiFight3D_single::MultiFight3D_single( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    world.init_world();
    printf( " world.defaultObjectShape %i \n", world.defaultObjectShape );

    warrior1 = world.makeWarrior( {0.0d,0.0d,0.0d}, {0.0d,0.0d,1.0d}, {0.0d,1.0d,0.0d}, world.defaultWarriorShape );

    zoom = 5.0;
    first_person = true;
    perspective  = true;

}

void MultiFight3D_single::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );
	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

    warrior1->gun_rot.set( (Vec3d)cam.rot.c );
    world.update_world( );



    for( auto o : world.objects ) {
        float glMat[16];
        glPushMatrix();
        Draw3D::toGLMat( o->lpos, o->lrot, glMat );
        glMultMatrixf( glMat );

        double t = raySphere( (Vec3d)cam.pos, (Vec3d)cam.rot.c, world.objR, o->lpos );
        if( ( t>0 ) && (t < 1000.0 ) ){
            //printf( " t %f  pos (%3.3f,%3.3f,%3.3f) \n", t, o->lpos.x, o->lpos.y, o->lpos.z );
            glCallList( world.defaultObjectHitShape );
        }

        glCallList( world.defaultObjectShape );
        glPopMatrix();
    }

    for( auto p : world.projectiles ) {
        float glMat[16];
        glPushMatrix();
        glTranslatef( p->pos.x, p->pos.y, p->pos.z );
        //printf( " %3.3f %3.3f %3.3f \n", p->pos.x, p->pos.y, p->pos.z);
        glCallList  ( world.defaultProjectileShape );
        glPopMatrix();
    }

};

/*
void MultiFight3D_single::mouseHandling( ){
    int mx,my;
    //Uint32 buttons = SDL_GetMouseState( &mx, &my );
    //SDL_WarpMouse(0, 0);
    SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    Quat4d q; q.fromTrackball( 0, 0, mx*0.001, my*0.001 );
    //qCamera.qmul( q );
    qCamera.qmul_T( q );
    camMat( );

    //defaultMouseHandling( mouseX, mouseY );
    //world.anchor.set( mouse_begin_x, mouse_begin_y );
}
*/


void MultiFight3D_single::keyStateHandling( const Uint8 *keys ){

    if( warrior1 != NULL ){
        //if( keys[ SDL_SCANCODE_W ] ){ warrior1->pos.add_mul( camMat.c, +0.1 ); }
        //if( keys[ SDL_SCANCODE_S ] ){ warrior1->pos.add_mul( camMat.c, -0.1 ); }
        //if( keys[ SDL_SCANCODE_A ] ){ warrior1->pos.add_mul( camMat.a, -0.1 ); }
        //if( keys[ SDL_SCANCODE_D ] ){ warrior1->pos.add_mul( camMat.a, +0.1 ); }

        RigidBody* rb = warrior1->asRigidBody();

        if( keys[ SDL_SCANCODE_W ] ){ rb->vel.add_mul( (Vec3d)cam.rot.c, +0.1 ); }
        if( keys[ SDL_SCANCODE_S ] ){ rb->vel.add_mul( (Vec3d)cam.rot.c, -0.1 ); }
        if( keys[ SDL_SCANCODE_A ] ){ rb->vel.add_mul( (Vec3d)cam.rot.a, -0.1 ); }
        if( keys[ SDL_SCANCODE_D ] ){ rb->vel.add_mul( (Vec3d)cam.rot.a, +0.1 ); }
        if( keys[ SDL_SCANCODE_SPACE ] ){ rb->vel.mul( 0.9 ); }

        cam.pos.set( (Vec3f)rb->pos );
    }

    //if( keys[ SDL_SCANCODE_W ] ){ camPos.add_mul( camMat.c, +0.1 ); }
	//if( keys[ SDL_SCANCODE_S ] ){ camPos.add_mul( camMat.c, -0.1 ); }
	//if( keys[ SDL_SCANCODE_A ] ){ camPos.add_mul( camMat.a, -0.1 ); }
	//if( keys[ SDL_SCANCODE_D ] ){ camPos.add_mul( camMat.a, +0.1 ); }
	//if( keys[ SDL_SCANCODE_W ] ){ camPos.z += +0.1; }
	//if( keys[ SDL_SCANCODE_S ] ){ camPos.z += -0.1; }
	//if( keys[ SDL_SCANCODE_A ] ){ camPos.x += +0.1; }
	//if( keys[ SDL_SCANCODE_D ] ){ camPos.x += -0.1; }
    if( keys[ SDL_SCANCODE_Q ] ){ qCamera.droll2(  +0.01 ); }
	if( keys[ SDL_SCANCODE_E ] ){ qCamera.droll2(  -0.01 ); }

};

void MultiFight3D_single::mouseHandling( ){
    int mx,my;
    //SDL_GetMouseState( &mouseX, &mouseY );
    Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    //if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {

        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        qCamera.qmul_T( q );
    //}
    //qCamera.qmul( q );
}



void MultiFight3D_single::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f:  warrior1->tryJump(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    warrior1->trigger = true;
                    break;
            }
            break;

        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    warrior1->trigger = false;
                    break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}




void MultiFight3D_single::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );
    drawCrosshair( 10 );
}

// ===================== MAIN

MultiFight3D_single * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new MultiFight3D_single( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















