
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
	//virtual void keyStateHandling( const Uint8 *keys );

	MultiFight3D_single( int& id, int WIDTH_, int HEIGHT_ );

};

MultiFight3D_single::MultiFight3D_single( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    world.init_world();
    printf( " world.defaultObjectShape %i \n", world.defaultObjectShape );
    zoom = 16.0;

}

void MultiFight3D_single::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );
	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

    for( auto o : world.objects ) {
        float glMat[16];
        glPushMatrix();
        Draw3D::toGLMat( o->lpos, o->lrot, glMat );
        glMultMatrixf( glMat );

        double t = raySphere( {0.0d,0.0d,0.0d}, camMat.c, 2.0, o->lpos );
        if( ( t>0 ) && (t < 1000.0 ) ){
            //printf( " t %f  pos (%3.3f,%3.3f,%3.3f) \n", t, o->lpos.x, o->lpos.y, o->lpos.z );
            glCallList( world.defaultObjectHitShape );
        }

        glCallList( world.defaultObjectShape );
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

/*
void MultiFight3D_single::keyStateHandling( const Uint8 *keys ){
};
*/


void MultiFight3D_single::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f:  warrior1->tryJump(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    pickParticle( world.picked );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
        */
    };
    AppSDL2OGL::eventHandling( event );
}




void MultiFight3D_single::drawHUD(){

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
















