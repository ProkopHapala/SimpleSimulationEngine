#ifndef  UDPClient_h
#define  UDPClient_h

#include "SDL_net.h"

class UDPClient{
	public:
	UDPsocket socket;
	IPaddress receiver;
	UDPpacket *packet;
 
	virtual int init( const char *host_str, Uint16 port, int nbytes ){
		if (             SDLNet_Init() < 0                               ){ printf( "SDLNet_Init:        %s\n", SDLNet_GetError());                        return -1; }
		if ( !( socket = SDLNet_UDP_Open(0) )                            ){ printf( "SDLNet_UDP_Open:    %s\n", SDLNet_GetError());                        return -2; }
		if (             SDLNet_ResolveHost( &receiver, host_str, port ) ){ printf( "SDLNet_ResolveHost(%s %d): %s\n", host_str, port, SDLNet_GetError()); return -3; }
		if ( !( packet = SDLNet_AllocPacket(nbytes))                     ){ printf( "SDLNet_AllocPacket: %s\n", SDLNet_GetError());                        return -4; }
		printf( "connected to (%i %i) \n", receiver.host, receiver.port );
		return 0;
	}

	void send_buff( ){
		packet->address.host = receiver.host;	// Set the destination host
		packet->address.port = receiver.port;	// And destination port		
		SDLNet_UDP_Send( socket, -1, packet );     // This sets the p->channel
	}

	virtual bool onSend(){
		printf(" Fill the buffer\n>");
		scanf("%s", (char *)packet->data);
		packet->len = strlen((char *)packet->data) + 1;
		if( packet->len > 1 ) return true; 
		printf( " empty \n" );
		return false;
	}

	virtual void send( ){ if ( onSend() ){ send_buff( ); } }
  
	virtual void close(){
		SDLNet_FreePacket( packet );
		SDLNet_Quit();
	}

};

#endif
