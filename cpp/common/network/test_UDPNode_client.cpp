
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <SDL2/SDL.h> 
#include <SDL2/SDL_opengl.h>
#include "SDL_net.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "UDPNode.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Vec2.h"
#include "geom2D.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"


// ==================== Declaration

class UDPClientApp : public AppSDL2OGL, public UDPNode {
	public:
	int frameCount = 0;
	bool loopEnd,STOP;
	double t=0,x=0,y=0,vx=0,vy=0;

	// ---- function declarations 

	//virtual void loop( int nframes );
	virtual void draw   ();

	virtual void onRecieve();
	virtual bool onSend();

	virtual void eventHandling( const SDL_Event& event );

	void add_char_to_buff( char ch );

	UDPClientApp( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

UDPClientApp::UDPClientApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
	init_UDP      ( 1, 2001, 512      );
	tryConnect_UDP( "localhost", 2000 );
}

void UDPClientApp::onRecieve(){
	//printPacketInfo( );
	int ibyte = 0;
	double t_ = ((double *)(packet->data+ibyte))[0];   ibyte+=sizeof(double);
	if( t_ > t ){
		t = t_;
		x  = ((double *)(packet->data+ibyte))[0];      ibyte+=sizeof(double);
 		y  = ((double *)(packet->data+ibyte))[0];      ibyte+=sizeof(double);
		vx = ((double *)(packet->data+ibyte))[0];      ibyte+=sizeof(double);
 		vy = ((double *)(packet->data+ibyte))[0];      ibyte+=sizeof(double);
		//((int    *)(packet->data+ibyte))[0] = iframe;  ibyte+=sizeof(int);
	}
}

bool UDPClientApp::onSend(){
	//printf( "sending: \n" );

	if( packet->len == 0 ) return false;
	packet->data[ packet->len ] = 0;
	packet->len++;
	printf( "sending: %s \n", (char*)packet->data );
	return true;
/*
	sprintf( (char*)packet->data, "client frame %06i\n", frameCount );
	packet->len = 19;
	printf("sending: %s", (char*)packet->data );
*/

};

void UDPClientApp::add_char_to_buff( char ch ){
	((char*)packet->data)[ packet->len ] = ch;
	packet->len++;
};

void UDPClientApp::eventHandling( const SDL_Event& event ){
	packet->len = 0;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_a: add_char_to_buff( 'a' ); break;
                case SDLK_d: add_char_to_buff( 'd' ); break;
            }
            break;
        case SDL_QUIT: quit(); break;
    };
	trySend( );
	AppSDL2OGL::eventHandling( event );
}

void UDPClientApp::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//trySend( );
	receiveQuedPackets( );

    glColor3f( 0.9f, 0.2f, 0.9f ); Draw2D::drawPointCross_d( { x, y }, 0.3 );
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine_d      ( { x, y }, { x+vx*0.1, y+vy*0.1 } );
	glColor3f( 0.2f, 0.2f, 0.2f ); Draw2D::drawRectangle   ( -8.0f, -8.0f, 8.0f, 8.0f, false );

	frameCount++;
};

// ===================== MAIN

UDPClientApp * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new UDPClientApp( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}

