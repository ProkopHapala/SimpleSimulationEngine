
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
#include "drawMath2D.h"

// ========= include from local app

#include "SailWarWorld.h"
#include "Projectile.h"
#include "Gun.h"
#include "Yacht2D.h"
#include "Frigate2D.h"

class SailWar_client : public AppSDL2OGL, public UDPNode, public SailWarWorld {
    public:
    Frigate2D* thisShip = NULL;

    int nProj = 0;

    bool keys_to_send[ 6 ];


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

	//init_UDP      ( 1, 2001, 512      );
	init_UDP      ( 1, 0, 512      );   // random port

	tryConnect_UDP( "localhost", 2000 );

    printf( " ==== main.setup \n" );
    init_world();
    //if( ) thisShip = ships.front();;       // FIXME - this is just temporary, later this should be specified by server
    //thisScreen->zoom = 100;
    printf( " ==== world.init DONE \n" );

}

void SailWar_client::draw(){

    glDisable  (GL_LIGHTING);
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glShadeModel(GL_FLAT);

    trySend();               // HERE WE SHOULD SEND WORLD STATE TO CLIENTS
    receiveQuedPackets( );   // HERE WE SHOULD READ INPUTS FROM CLIENTS

    if( thisShip != NULL ){

        printf( "net_id %i   thisShip %i \n", net_id, thisShip );
        printf( "pos (%3.3f,%3.3f) phi %3.3f \n", thisShip->pos, thisShip->phi );

        glPushMatrix();
        glRotatef   ( (float)-thisShip->phi*180/M_PI + 90, 0.0f, 0.0f, 1.0f );
        glTranslatef( (float)-thisShip->pos.x, (float)-thisShip->pos.y, 0 );
        glColor3f( 0.55f, 0.55f, 0.55f ); Draw2D::drawGrid( -50.0, -50.0, +50.0, +50.0, 5.0, 5.0 );

        for( auto ship : ships ) {
            glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->drawHitBox( );
            glColor3f( 0.8f, 0.8f, 0.8f ); 	ship->draw_shape( );
            glColor3f( 0.2f, 0.2f, 0.2f );  ship->draw( );
        }

        for( int i=0; i<nProj; i++ ){
            Projectile * p = projectiles[ i ];
            //printf( " draw %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f \n", i, p->pos.x,p->pos.y,p->pos.z, p->vel.x,p->vel.y,p->vel.z, p->mass );
            p->draw();
        }

        //glPushMatrix();
        //glRotatef   ( (float)thisShip->phi*180/M_PI, 0.0f, 0.0f, 1.0f );
        //Draw2D::drawGrid( -5.0, -5.0, +5.0, +5.0, 1.0, 1.0 );
        glPopMatrix();

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

    if( keys[ SDL_SCANCODE_LEFT  ] ){ keys_to_send[0] = true; }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ keys_to_send[1] = true; }
	if( keys[ SDL_SCANCODE_UP    ] ){ keys_to_send[2] = true; }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ keys_to_send[3] = true; }

};

void SailWar_client::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
				//case SDLK_KP_1:     thisShip->fire_left ( &projectiles ); break;
				//case SDLK_KP_2:     thisShip->fire_right( &projectiles ); break;
                case SDLK_KP_1:     keys_to_send[4] = true; break;
				case SDLK_KP_2:     keys_to_send[5] = true; break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
};


//////////////////////////////////
//      overide UDPNode
//////////////////////////////////

void SailWar_client::onRecieve( ){

    //printPacketInfo();

    char * buff = ( char* ) packet->data;

    net_id     = (*(int*)buff);   buff += sizeof( int );
    int nShips = (*(int*)buff);   buff += sizeof( int );
    int nShipsLeft = nShips - ships.size();
    if( nShipsLeft > 0 ){
        for( int i=0; i<nShipsLeft; i++ ){
            SailWarWorld::makeShip( { 0.0, -0.0}, M_PI*0.6, "data/FrigateType.txt", defaultShipShape, defaultCollisionShape );
        }
    }

    auto it_ship  = ships.begin();
    while( it_ship != ships.end  ( ) ) {
        Frigate2D * ship = *it_ship;
        buff = ship->fromBytes( buff );
        ++it_ship;
    }

    nProj       = (*(int*)buff);   buff += sizeof( int );
    int nProjSz = projectiles.size();
    for( int i=0; i<nProj; i++ ){
        Projectile * p;
        if( i >= nProjSz ){
            p = new Projectile( );
            projectiles.push_back( p );
            //printf( "new projectile i %i nProjSz %i ", i, nProjSz );
        }else{
            p = projectiles[ i ];
        }
        buff = p->fromBytes( buff );
        //printf( " onRecieve %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) %3.3f \n", i, p->pos.x,p->pos.y,p->pos.z, p->vel.x,p->vel.y,p->vel.z, p->mass );
    }

    thisShip   = ships[net_id];

};

bool SailWar_client::onSend   ( ){
    bool * keys = (bool*) packet->data;
    for( int i=0; i<6; i++ ){
        keys[i] = keys_to_send[i];
        keys_to_send[i]  = false;
    }
    packet->len = 6;
    return true;
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


