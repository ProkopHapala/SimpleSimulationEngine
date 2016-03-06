
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include <SDL2/SDL.h>
#include "SDL_net.h"

#include "UDPClient.h"

// super simple example
/*
UDPClient client;
int main(int argc, char **argv){
	client.init( "localhost", 2000, 0, 512 );
	while (true){
		//printf( "client try\n" );
		client.send( );
		SDL_Delay( 5 );
	}
}
*/

class UDPClient_mod : public UDPClient {
	public:
	int iframe=0;
	double t=0,x=0,y=0,vx=1,vy=0;

	double move( double dt ){
		double G = 9.81; 
		vy -= G * dt;
		x  += vx * dt;
		y  += vy * dt;
		if( ( x >  10 ) && ( vx > 0 ) ){ vx = -vx; }
		if( ( x < -10 ) && ( vx < 0 ) ){ vx = -vx; }
		if( ( y >  10 ) && ( vy > 0 ) ){ vy = -vy; }
		if( ( y < -10 ) && ( vy < 0 ) ){ vy = -vy; }
		iframe++;
		t+=dt;
	}

	virtual bool onSend(){
		move( 0.01 );
		const int expected_length = sizeof(double)*3 + sizeof(int)*1;
		((double *)packet->data)[0] = x;
 		((double *)packet->data)[1] = y;
		((double *)packet->data)[2] = t;
		((int *)&(((double *)packet->data)[3]))[0] = iframe;
		packet->len  = expected_length;
	}

};

UDPClient_mod client;

int main(int argc, char **argv){
	int frameCount = 0;
	//client.init( "localhost", 2000, 512 );

	client.init_UDP     ( 1, 0, 512 );
	client.tryConnect_UDP( "localhost", 2000 );

	while (true){
		printf( "test_UDPClient %04i: try\n", frameCount );
		//printf( "client try\n" );
		client.send( );
		SDL_Delay( 20 );
		frameCount++;
	}
}


