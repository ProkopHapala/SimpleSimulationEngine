#ifndef  UDPClient_h
#define  UDPClient_h

#include "SDL_net.h"

class UDPClient{
	public:
	UDPsocket sd;
	IPaddress srvadd;
	UDPpacket *p;
 
	virtual int init( const char *host, Uint16 port, int socket, int nbytes ){
		if (        SDLNet_Init() < 0                         ){ printf( "SDLNet_Init:        %s\n", SDLNet_GetError());                     return 1; }
		if ( !(sd = SDLNet_UDP_Open(socket) )                 ){ printf( "SDLNet_UDP_Open:    %s\n", SDLNet_GetError());                     return 2; }
		if (        SDLNet_ResolveHost( &srvadd, host, port ) ){ printf( "SDLNet_ResolveHost(%s %d): %s\n", host, port, SDLNet_GetError());	 return 3; }
		if ( !( p = SDLNet_AllocPacket(512))                  ){ printf( "SDLNet_AllocPacket: %s\n", SDLNet_GetError());                     return 4; }
		return 0;
	}

	void send_buff( ){
		p->address.host = srvadd.host;	// Set the destination host
		p->address.port = srvadd.port;	// And destination port		
		SDLNet_UDP_Send(sd, -1, p);     // This sets the p->channel
	}

	virtual bool onSend(){
		printf(" Fill the buffer\n>");
		scanf("%s", (char *)p->data);
		p->len = strlen((char *)p->data) + 1;
		if( p->len > 1 ) return true; 
		printf( " empty \n" );
		return false;
	}

	virtual void send( ){ if ( onSend() ){ send_buff( ); } }
  
	virtual void close(){
		SDLNet_FreePacket(p);
		SDLNet_Quit();
	}

};

#endif
