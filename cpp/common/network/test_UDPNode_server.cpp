
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include <SDL2/SDL.h>
#include "SDL_net.h"

#include "UDPNode.h"

class TestServer : public UDPNode {
	public:
	int iframe;
	double t=0,x=0,y=0,vx=-1,vy=0;
	double friction    = -0.5;
	double restitution =  0.5;
	double gravity     = -9.81;

	double move( double dt ){
		 
		vx += ( vx * friction           ) * dt;
		vy += ( vy * friction + gravity ) * dt;

		if( ( x >  10 ) && ( vx > 0 ) ){ vx = -vx; }
		if( ( x < -10 ) && ( vx < 0 ) ){ vx = -vx; }
		if( ( y >  10 ) && ( vy > 0 ) ){ vy = -vy; }
		if( ( y < -10 ) && ( vy < 0 ) ){ vy = -vy; }

		x  += vx * dt;
		y  += vy * dt;

		iframe++;
		t+=dt;
	}

	virtual bool onSend(){
		int ibyte = 0;	
		((double *)(packet->data+ibyte))[0] = t;       ibyte+=sizeof(double);	
		((double *)(packet->data+ibyte))[0] = x;       ibyte+=sizeof(double);
 		((double *)(packet->data+ibyte))[0] = y;       ibyte+=sizeof(double);
		((double *)(packet->data+ibyte))[0] = vx;      ibyte+=sizeof(double);
 		((double *)(packet->data+ibyte))[0] = vy;      ibyte+=sizeof(double);
		//((int    *)(packet->data+ibyte))[0] = iframe;  ibyte+=sizeof(int);
		packet->len = ibyte;
		return true;
	}

	virtual void onRecieve(){ 
		printPacketInfo( );
		//printf( "recieved: %s \n", (char*)packet->data );
		for( int i=0; i<packet->len; i++ ){
			char key = ((char*)packet->data)[i];
            switch( key ){
                case 'a': vy+=3.0; vx+=3.0; break;
                case 'd': vy+=3.0; vx-=3.0; break;
            }
		}
	}

};

TestServer server;

int main(int argc, char **argv){

	server.init_UDP   ( 0, 2000, 512       );
	
	while (true){
		//printf( "server try %i :\n", server.iframe );
		server.receiveQuedPackets( );
		server.move( 0.01 );

		if( !server.connected ){ 
			server.tryConnect_UDP( "localhost", 2001  );
		}

		server.trySend(   );
		SDL_Delay( 10 );
	}
}


