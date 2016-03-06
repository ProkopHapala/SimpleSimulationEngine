
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

class SailWar_server : public AppSDL2OGL, public UDPNode, public SailWarWorld {
    public:
// ==== overide AppSDL2OGL
    virtual void draw    ( );
    virtual void drawHUD ( );
	SailWar_server( int& id, int WIDTH_, int HEIGHT_ );

// ==== overide UDPNode
	virtual void onRecieve( );
	virtual bool onSend   ( );

// ==== overide SailWarWorld

};

//////////////////////////////////
//      overide AppSDL2OGL
//////////////////////////////////

SailWar_server::SailWar_server( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

	init_UDP      ( 1, 2000, 512      );
	tryConnect_UDP( "localhost", 2001 );

    printf( " ==== main.setup \n" );
    init_world();
    //thisScreen->zoom = 100;
    printf( " ==== world.init DONE \n" );

}

void SailWar_server::draw(){

    glDisable  (GL_LIGHTING);
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glShadeModel(GL_FLAT);

    //tryReceive();     // HERE WE SHOULD READ INPUTS FROM CLIENTS
	update_world();
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

void SailWar_server::drawHUD(){
    /*
    float bar_y =  10.0f;
    float bar_x = 100.0f;

    Vec2d compass_pos; compass_pos.set( WIDTH*0.9, HEIGHT*0.9 );
	glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( compass_pos, 10 );
	glColor3f( 0.2f, 0.5f, 0.2f );  Draw2D::drawVecInPos_d  ( ( *(Vec2d*)&(wind_speed))*10.1, compass_pos );
	glColor3f( 0.2f, 0.2f, 0.8f );  Draw2D::drawVecInPos_d  ( watter_speed*10.1, compass_pos );
    */
};

//////////////////////////////////
//      overide UDPNode
//////////////////////////////////

void SailWar_server::onRecieve( ){

};

bool SailWar_server::onSend   ( ){

};

// ============== main

SailWar_server * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new SailWar_server( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}


