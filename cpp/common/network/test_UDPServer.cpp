
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <SDL2/SDL.h> 
#include "SDL_net.h"

#include "UDPServer.h"

// super simple example
/*
UDPServer server;
int main(int argc, char **argv){
	server.init( 2000, 512 );
	while (true){
		printf( "server try\n" );
		server.tryReceive( );
		SDL_Delay( 20 );
	}
}
*/

class UDPServer_mod : public UDPServer {
	public:
	virtual void onRecieve(){
		const int expected_length = sizeof(double)*3 + sizeof(int)*1;
		if( p->len >= expected_length ){
			double x   = ((double *)p->data)[0];
			double y   = ((double *)p->data)[1];
			double t   = ((double *)p->data)[2];
			int iframe = ((int *)&(((double *)p->data)[3]))[0]; 
			printf( " %f %f %f %i \n", x, y, t, iframe );
		}
	}
};

UDPServer_mod server;

int main(int argc, char **argv){
	server.init( 2000, 512 );
	while (true){
		printf( "server try\n" );
		//server.tryReceive( );
		server.receiveQuedPackets( );
		SDL_Delay( 20 );
	}
}
