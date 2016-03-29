
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

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

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //glEnable( GL_LIGHTING );
    //glCallList( world.defaultObjectShape );
	glDisable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

    for( auto o : world.objects ) {
        float glMat[16];
        glPushMatrix();
        Draw3D::toGLMat( o->lpos, o->lrot, glMat );
        glMultMatrixf( glMat );
        glCallList( world.defaultObjectShape );
        glPopMatrix();
    }

};

/*
void NonInert_seats::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
    world.anchor.set( mouse_begin_x, mouse_begin_y );
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
	int junk;
	thisApp = new MultiFight3D_single( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















