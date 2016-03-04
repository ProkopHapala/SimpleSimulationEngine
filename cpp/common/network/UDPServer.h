#ifndef  UDPServer_h
#define  UDPServer_h

#include "SDL_net.h"

class UDPServer{
	public:
	UDPsocket sd;			// Socket descriptor 
	UDPpacket *p = NULL;	// Pointer to packet memory 

	virtual void onRecieve(){
		printf("UDP Packet incoming\n");
		printf("\tChan:    %d\n", p->channel);
		printf("\tData:    %s\n", (char *)p->data);
		printf("\tLen:     %d\n", p->len);
		printf("\tMaxlen:  %d\n", p->maxlen);
		printf("\tStatus:  %d\n", p->status);
		printf("\tAddress: %x %x\n", p->address.host, p->address.port);
	}

	bool tryReceive( ){
		if ( SDLNet_UDP_Recv( sd, p ) )	{
			onRecieve();
		}
	}

	bool receiveQuedPackets( ){
		while ( SDLNet_UDP_Recv( sd, p ) )	{
			onRecieve();
		}
	}

	virtual int init( int socket, int nbytes ){   
		if ( SDLNet_Init() < 0 )                        { printf("SDLNet_Init: %s\n", SDLNet_GetError() );        return 1; }
		sd = SDLNet_UDP_Open(2000);	   if ( sd == NULL ){ printf("SDLNet_UDP_Open: %s\n", SDLNet_GetError() );	  return 2; }
		p  = SDLNet_AllocPacket(512);  if ( p == NULL  ){ printf("SDLNet_AllocPacket: %s\n", SDLNet_GetError() ); return 3; }
		return 0;
	}

	virtual void close(){  // Clean and exit
		SDLNet_FreePacket(p);
		SDLNet_Quit();
	}

};

#endif
