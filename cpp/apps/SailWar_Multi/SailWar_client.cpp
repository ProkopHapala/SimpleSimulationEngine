
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <SDL2/SDL_net.h>

// ========= include from common
#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

//#include "drawUtils.h"

#include "Body.h"
#include "Body2D.h"
#include "AeroSurf2D.h"
#include "geom2D.h"

// ========= include 3rd party library dependent

#include  "AppSDL2OGL.h"
#include  "UDPNode.h"

#include "GridMap2D.h"
#include "drawMath.h"
#include "drawMath2D.h"

// ========= include from local app

#include "SailWarWorld.h"
#include "Projectile.h"
#include "Gun.h"
#include "Yacht2D.h"
#include "Frigate2D.h"

class SailWar_client : public AppSDL2OGL, public UDPNode, public SailWarWorld {
    public:
    Frigate2D* thisShip;


// ==== overide AppSDL2OGL
    virtual void draw            (                                 );
    virtual void drawHUD         (                                 );
    virtual void eventHandling   ( const SDL_Event& event          );
	virtual void keyStateHandling( const Uint8 *keys               );

	SailWar_client( int& id, int WIDTH_, int HEIGHT_ );

// ==== overide UDPNode
	virtual void onRecieve( );
	virtual bool onSend   ( );

// ==== overide SailWarWorld

};

//////////////////////////////////
//      overide AppSDL2OGL
//////////////////////////////////

SailWar_client::SailWar_client( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

	init_UDP      ( 1, 2000, 512      );
	tryConnect_UDP( "localhost", 2001 );

    printf( " ==== main.setup \n" );
    init_world();
    thisShip = ships.front();;       // FIXME - this is just temporary, later this should be specified by server
    //thisScreen->zoom = 100;
    printf( " ==== world.init DONE \n" );

}

void SailWar_client::draw(){

    glDisable  (GL_LIGHTING);
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glShadeModel(GL_FLAT);

    tryReceive();     // HERE WE SHOULD READ INPUTS FROM CLIENTS


	//update_world();
	//trySend();      // HERE WE SHOULD SEND WORLD STATE TO CLIENTS

	for( auto ship : ships ) {
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->drawHitBox( );
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->draw_shape( );
        glColor3f( 0.2f, 0.2f, 0.2f );  ship->draw( );
	}
	for( auto p : projectiles ) {
		p->draw();
	}

};

void SailWar_client::drawHUD(){
    /*
    float bar_y =  10.0f;
    float bar_x = 100.0f;

    Vec2d compass_pos; compass_pos.set( WIDTH*0.9, HEIGHT*0.9 );
	glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( compass_pos, 10 );
	glColor3f( 0.2f, 0.5f, 0.2f );  Draw2D::drawVecInPos_d  ( ( *(Vec2d*)&(wind_speed))*10.1, compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  Draw2D::drawVecInPos_d  ( watter_speed*10.1, compass_pos );
    */
};


void SailWar_client::keyStateHandling( const Uint8 *keys ){

    //if( keys[ SDL_SCANCODE_LEFT  ] ){ camX0 -= camStep; }
	//if( keys[ SDL_SCANCODE_RIGHT ] ){ camX0 += camStep; }
	//if( keys[ SDL_SCANCODE_UP    ] ){ camY0 += camStep; }
	//if( keys[ SDL_SCANCODE_DOWN  ] ){ camY0 -= camStep; }

    if( keys[ SDL_SCANCODE_LEFT  ] ){ thisShip->rudder.setAngle( thisShip->rudder.phi + 0.01 );  }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ thisShip->rudder.setAngle( thisShip->rudder.phi - 0.01 );  }
	if( keys[ SDL_SCANCODE_UP    ] ){ thisShip->mast.setAngle  ( thisShip->mast.phi   + 0.01 );  }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ thisShip->mast.setAngle  ( thisShip->mast.phi   - 0.01 );  }

};

void SailWar_client::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
				case SDLK_KP_1:     thisShip->fire_left ( &projectiles ); break;
				case SDLK_KP_2:     thisShip->fire_right( &projectiles ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
};


//////////////////////////////////
//      overide UDPNode
//////////////////////////////////

void SailWar_client::onRecieve( ){

};

bool SailWar_client::onSend   ( ){

};

// ============== main

SailWar_client * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new SailWar_client( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}


