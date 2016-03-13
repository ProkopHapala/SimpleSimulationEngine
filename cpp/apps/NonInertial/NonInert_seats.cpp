
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "drawMath2D.h"

#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "NonInertWorld.h"
#include "Warrior2D.h"
#include "Projectile2D.h"

class NonInert_seats : public AppSDL2OGL {
	public:
    NonInertWorld world;
    double dvel = 10.0;

    Warrior2D *warrior1,*warrior2;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

	NonInert_seats( int& id, int WIDTH_, int HEIGHT_ );

};

NonInert_seats::NonInert_seats( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    world.init_world();

    world.makeWarrior( {-5.0d,0.0d}, 0.0d, "", world.defaultWarriorShape );
    world.makeWarrior( {+5.0d,0.0d}, 0.0d, "", world.defaultWarriorShape );

    warrior1 = world.warriors[0];
    warrior2 = world.warriors[1];

    zoom = 16.0;
}

void NonInert_seats::draw(){

    glClearColor( 0.1f, 0.1f, 0.1f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glDisable( GL_DEPTH_TEST);

    //printf( "zoom %f \n", zoom );
    world.update_world();

    //warrior1->pos.set( 0.0,  -5.0 );    warrior1->vel.set( 0.0,  +2.0 );
    //warrior2->pos.set( 0.0, -10.0 );    warrior2->vel.set( 0.0,  -2.0 );

    //warrior2->pos.set( mouse_begin_x, mouse_begin_y );    warrior2->vel.set( 0.0, -1.0 );
    //warrior2->clean_temp ( );
    //world.addEnviroForces(  warrior2->pos,  warrior2->vel,  warrior2->force );

    if( world.projectiles.size() > 1 ) printf( " frame %i : %i projectiles in air \n", frameCount, world.projectiles.size() );

    glPushMatrix();
    glRotatef( -world.phi * 180 / M_PI, 0.0f, 0.0f, 1.0f );

    //glColor3f( 0.8, 0.8, 0.8 ); Draw2D::drawLine_d( {0.0d,0.0d}, { world.rot.x*world.rmax, world.rot.y*world.rmax } );
    glColor3f( 0.8, 0.8, 0.8 ); Draw2D::drawRectangle( -20.0, -3.0, 20.0, 3.0, true );
    glColor3f( 1.0, 1.0, 1.0 );	Draw2D::drawCircle_d( world.center, world.rmax+1.0, 64, true );
	//glColor3f( 0.8, 0.8, 0.8 ); Draw2D::drawPointCross_d( world.center, 0.5 );
	//glColor3f( 0.8, 0.8, 0.8 ); Draw2D::drawLine_d( {0.0d,world.rmax+1.0}, { 0, world.rmax+1.5 } );

    for( auto w : world.warriors ) {
        glColor3f( 0.8f, 0.8f, 0.8f ); 	w->draw_shape( );

        //w->clean_temp ( );
        //world.addEnviroForces( w->pos, w->vel, w->force );
        glColor3f( 0.1, 0.1, 0.1 ); Draw2D::drawPointCross_d( w->pos, 0.1 );
        glColor3f( 0.0, 0.0, 0.8 ); Draw2D::drawVecInPos_d  ( w->vel     * 1.0,  w->pos );
        glColor3f( 0.8, 0.0, 0.0 ); Draw2D::drawVecInPos_d  ( w->force   * 10.0, w->pos );
        glColor3f( 0.0, 0.5, 0.0 ); Draw2D::drawVecInPos_d  ( w->gun_rot * 1.0, w->pos );
    }

    for( auto p : world.projectiles ) {
        glColor3f( 0.0f, 0.5f, 0.0f ); 	Draw2D::drawVecInPos_d( p->vel*0.02, p->pos );
    }

    glPopMatrix();
};

/*
void NonInert_seats::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
    world.anchor.set( mouse_begin_x, mouse_begin_y );
}
*/

void NonInert_seats::keyStateHandling( const Uint8 *keys ){

    Vec2d dvx,dvy;
    dvx.set_mul ( world.rot, world.dt*10 );
    dvy.set_perp( dvx );
    if( keys[ SDL_SCANCODE_W ] ){ warrior1->vel.add( dvy ); }
	if( keys[ SDL_SCANCODE_S ] ){ warrior1->vel.sub( dvy ); }
	if( keys[ SDL_SCANCODE_D ] ){ warrior1->vel.add( dvx ); }
    if( keys[ SDL_SCANCODE_A ] ){ warrior1->vel.sub( dvx ); }
    if( keys[ SDL_SCANCODE_Q ] ){ warrior1->rotate_gun( -0.03 ); }
    if( keys[ SDL_SCANCODE_E ] ){ warrior1->rotate_gun( +0.03 ); }
    if( keys[ SDL_SCANCODE_R ] ){ world.fireProjectile( warrior1 ); }

    if( keys[ SDL_SCANCODE_I ] ){ warrior2->vel.add( dvy ); }
	if( keys[ SDL_SCANCODE_K ] ){ warrior2->vel.sub( dvy ); }
	if( keys[ SDL_SCANCODE_L ] ){ warrior2->vel.add( dvx ); }
    if( keys[ SDL_SCANCODE_J ] ){ warrior2->vel.sub( dvx ); }
    if( keys[ SDL_SCANCODE_U ] ){ warrior2->rotate_gun( -0.03 ); }
    if( keys[ SDL_SCANCODE_O ] ){ warrior2->rotate_gun( +0.03 ); }
    if( keys[ SDL_SCANCODE_P ] ){ world.fireProjectile( warrior2 ); }

    //case SDLK_r:  world.fireProjectile( warrior1 ); break;

};


void NonInert_seats::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_f:  warrior1->tryJump(); break;
                case SDLK_h:  warrior1->tryJump(); break;
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




void NonInert_seats::drawHUD(){

}

// ===================== MAIN

NonInert_seats * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new NonInert_seats( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















