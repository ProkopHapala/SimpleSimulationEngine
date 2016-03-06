#ifndef  UDPServer_h
#define  UDPServer_h

#include "SDL_net.h"

class UDPServer{
	public:
	UDPsocket socket;			// Socket descriptor 
	UDPpacket *packet = NULL;	// Pointer to packet memory 

	// ---- incoming

	void printPacketInfo( ){
		printf("UDP Packet incoming\n");
		printf("\tChan:    %d\n", packet->channel);
		printf("\tData:    %s\n", (char *)packet->data);
		printf("\tLen:     %d\n", packet->len);
		printf("\tMaxlen:  %d\n", packet->maxlen);
		printf("\tStatus:  %d\n", packet->status);
		printf("\tAddress: %x %x\n", packet->address.host, packet->address.port);
	}

	virtual void onRecieve(){
		printPacketInfo( );
	}

	bool tryReceive( ){
		if ( SDLNet_UDP_Recv( socket, packet ) )	{
			onRecieve();
		}
	}

	bool receiveQuedPackets( ){
		while ( SDLNet_UDP_Recv( socket, packet ) )	{
			onRecieve();
		}
	}
	virtual int init( int port, int nbytes ){   
		if (   SDLNet_Init() < 0                     ){ printf("SDLNet_Init: %s\n", SDLNet_GetError() );           return -1; }
		if ( !(socket = SDLNet_UDP_Open(port)      ) ){ printf("SDLNet_UDP_Open: %s\n", SDLNet_GetError() );	   return -2; }
		if ( !(packet = SDLNet_AllocPacket(nbytes) ) ){ printf("SDLNet_AllocPacket: %s\n", SDLNet_GetError() );    return -3; }
		return 0;
	}

	virtual void close(){  // Clean and exit
		SDLNet_FreePacket(packet);
		SDLNet_Quit();
	}

};

#endif
