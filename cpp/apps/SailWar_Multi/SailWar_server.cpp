
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
//#include  "UDPNode.h"
#include  "UDPServer.h"

//#include "GridMap2D.h"
#include "Draw2D.h"

// ========= include from local app

#include "SailWarWorld.h"
#include "Projectile.h"
#include "Gun.h"
#include "Yacht2D.h"
#include "Frigate2D.h"

class SailWar_server : public AppSDL2OGL, public UDPServer, public SailWarWorld {
    public:
// ==== overide AppSDL2OGL
    virtual void draw    ( );
    virtual void drawHUD ( );
	SailWar_server( int& id, int WIDTH_, int HEIGHT_ );

// ==== overide UDPNode
	virtual void onRecieve( int iClient );
	virtual bool onSend   ( int iClient );
	virtual int  onConnect( IPaddress client );

// ==== overide SailWarWorld

};

//////////////////////////////////
//      overide AppSDL2OGL
//////////////////////////////////

SailWar_server::SailWar_server( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

	init_UDP      ( 2000, 512, 8 );

    printf( " ==== main.setup \n"    );
    init_world();
    //makeShip( { randf(-5.0,5.0), randf(-5.0,5.0)}, M_PI*0.6, "data/FrigateType.txt", defaultShipShape, defaultCollisionShape );
    //ships[ 0 ]->fire_left ( &projectiles );
    //thisScreen->zoom = 100;
    printf( " ==== world.init DONE \n" );

}

void SailWar_server::draw(){

    //printf( " ======== frame %i \n", frameCount );
    //printf( "fromClients()\n" );
    fromClients();    // HERE WE SHOULD READ INPUTS FROM CLIENTS
    //printf( "update_world()\n" );
	update_world();
	//printf( "toClients();\n" );
	toClients();      // HERE WE SHOULD SEND WORLD STATE TO CLIENTS

    //printf( "\n" );

    glDisable  (GL_LIGHTING);
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glShadeModel(GL_FLAT);

	for( auto ship : ships ) {
        glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->drawHitBox( );
        //glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->draw_shape( );
        glColor3f( 0.8f, 0.8f, 0.8f ); 	Draw2D::drawShape ( ship->pos, ship->rot, ship->shape );
        glColor3f( 0.2f, 0.2f, 0.2f );  ship->draw( );
	}
	for( auto p : projectiles ) {
		p->draw();
		//printf( " draw (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f \n", p->pos.x,p->pos.y,p->pos.z, p->vel.x,p->vel.y,p->vel.z, p->mass );
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

void SailWar_server::onRecieve( int iClient ){
    //printPacketInfo();
    Frigate2D * thisShip = ships[ iClient ];
    bool      * keys     = (bool*) packet->data;
    if( keys[ 0 ] ){ thisShip->rudder.setAngle( thisShip->rudder.phi + 0.01 );  }
	if( keys[ 1 ] ){ thisShip->rudder.setAngle( thisShip->rudder.phi - 0.01 );  }
	if( keys[ 2 ] ){ thisShip->mast.setAngle  ( thisShip->mast.phi   + 0.01 );  }
	if( keys[ 3 ] ){ thisShip->mast.setAngle  ( thisShip->mast.phi   - 0.01 );  }
    if( keys[ 4 ] ){ thisShip->fire_left ( &projectiles );                      }
    if( keys[ 5 ] ){ thisShip->fire_right( &projectiles );                      }
};

bool SailWar_server::onSend   ( int iClient ){

    //printf( "onSend : iClinet %i \n", iClient  );
    char * buff = (char*)packet->data;

    (*(int*)buff) = iClient;       buff += sizeof( int );

    (*(int*)buff) = ships.size();  buff += sizeof( int );
    /*
    auto it_ship    = ships.begin();
    while( it_ship != ships.end  ( ) ) {
        Frigate2D * ship = *it_ship;
        buff = ship->toBytes( buff );
        ++it_ship;
    }
    */

    for( auto ship : ships ) {
        //printf( " onSend (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f \n", p->pos.x,p->pos.y,p->pos.z, p->vel.x,p->vel.y,p->vel.z, p->mass );
        buff = ship->toBytes( buff );
	}

    printf( " =============== frame %i \n", frameCount );
    (*(int*)buff) =  projectiles.size();  buff += sizeof( int );
    for( auto p : projectiles ) {
        //printf( " onSend (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f \n", p->pos.x,p->pos.y,p->pos.z, p->vel.x,p->vel.y,p->vel.z, p->mass );
        buff = p->toBytes( buff );
	}

    packet->len = buff - (char*)packet->data;
};

int SailWar_server::onConnect( IPaddress client ){
    //makeShip( { 3.0, -3.0}, M_PI*0.6, "data/FrigateType.txt", defaultShipShape, defaultCollisionShape );
    makeShip( { randf(-5.0,5.0), randf(-5.0,5.0)}, M_PI*0.6, "data/FrigateType.txt", defaultShipShape, defaultCollisionShape );
    printf( " new Player [ %i ] initialized \n", ships.size() );
    return UDPServer::onConnect( client );
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


