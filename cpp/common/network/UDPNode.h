#ifndef  UDPNode_h
#define  UDPNode_h

#include "SDL_net.h"

class UDPNode{
	public:
	UDPsocket  in_socket;	
	UDPsocket  out_socket;	 
	UDPpacket *packet = NULL;
	IPaddress  server_address;
	bool       connected = false;

	Uint16  net_id; 
	//char*   host_IP;
	//Uint16  in_port, out_port; 
 
	// ---- outgoing

	void send_buff( ){
		packet->address.host = server_address.host;
		packet->address.port = server_address.port;
		SDLNet_UDP_Send( out_socket, -1, packet );
	}

	virtual bool onSend(){
		printf(" Fill the buffer\n>");
		scanf("%s", (char *)packet->data);
		packet->len = strlen((char *)packet->data) + 1;
		if( packet->len > 1 ) return true; 
		printf( " empty \n" );
		return false;
	}

	virtual void trySend( ){ 
		if ( onSend() ){ 
			//printf( "sending bytes %i \n", packet->len );
			send_buff( ); 
		} 
	}

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
		if ( SDLNet_UDP_Recv( in_socket, packet ) )	{
			onRecieve();
		}
	}

	bool receiveQuedPackets( ){
		while ( SDLNet_UDP_Recv( in_socket, packet ) )	{
			onRecieve();
		}
	}

	// ---- other

	virtual int init_UDP( Uint16 net_id_, Uint16 in_port, Uint16 nbuff ){
		if (                 SDLNet_Init()             < 0   ){  printf( "SDLNet_Init:        %s\n", SDLNet_GetError());       return -1; }
		if ( !( packet     = SDLNet_AllocPacket( nbuff   ) ) ){  printf( "SDLNet_AllocPacket: %s\n", SDLNet_GetError() );      return -2; }
		if ( !( in_socket  = SDLNet_UDP_Open(   in_port  ) ) ){  printf( "SDLNet_UDP_Open  in_port: %s\n", SDLNet_GetError()); return -3; }			
		return 0;
	}

	virtual int tryConnect_UDP( const char *host, Uint16 out_port ){
		printf( "tryConnect_UDP( %s, %i )\n", host, out_port );
		//if ( !( in_socket  = SDLNet_UDP_Open    ( 0        )                      ) ){  printf( "SDLNet_UDP_Open  port 0 : %s\n", SDLNet_GetError());                    return -1; }
		if (            0 != SDLNet_ResolveHost ( &server_address, host, out_port ) ){  printf( "SDLNet_ResolveHost(%s %d): %s\n", host, out_port, SDLNet_GetError());	 return -2; }
		//if ( !( out_socket = SDLNet_UDP_Open    ( out_port )                      ) ){  printf( "SDLNet_UDP_Open out_port: %s\n", SDLNet_GetError() );                   return -3; }	
		connected = true;
		return 0;
	}
  
	virtual void close_UDP(){
		if( packet  != NULL ) SDLNet_FreePacket(packet);
		SDLNet_Quit();
	}

};

#endif
