
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
	int iframe;
	double t,x,y,vx,vy;

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
		((double *)p->data)[0] = x;
 		((double *)p->data)[1] = y;
		((double *)p->data)[2] = t;
		((int *)&(((double *)p->data)[3]))[0] = iframe;
		p->len  = expected_length;
	}

};

UDPClient_mod client;

int main(int argc, char **argv){
	client.init( "localhost", 2000, 0, 512 );

	client.iframe;
	client.t=0;
	client.x=0;
	client.y=0;
	client.vx=1;
	client.vy=0;

	while (true){
		//printf( "client try\n" );
		client.send( );
		SDL_Delay( 20 );
	}
}


