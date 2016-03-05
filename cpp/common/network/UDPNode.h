#ifndef  UDPNode_h
#define  UDPNode_h

#include "SDL_net.h"

class UDPNode{
	public:
	//UDPsocket  in_socket;	
	//UDPsocket  out_socket;	
 	UDPsocket  socket;
	UDPpacket *packet = NULL;
	IPaddress  receiver;
	bool       connected = false;

	Uint16  net_id; 
	//char*   host_IP;
	//Uint16  in_port, out_port; 
 
	// ---- outgoing

	void send_buff( ){
		packet->address.host = receiver.host;
		packet->address.port = receiver.port;
		//SDLNet_UDP_Send( in_socket, -1, packet );
		SDLNet_UDP_Send( socket, -1, packet );
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
		//printf( " tryReceive: \n" );
		if ( SDLNet_UDP_Recv( socket, packet ) )	{
			onRecieve();
		}
	}

	bool receiveQuedPackets( ){
		//printf( " receiveQuedPackets: \n" );
		while ( SDLNet_UDP_Recv( socket, packet ) )	{
			onRecieve();
		}
	}

	// ---- other

	virtual int init_UDP( Uint16 net_id_, Uint16 in_port, Uint16 nbuff ){
		net_id = net_id_;
		if (              SDLNet_Init()             < 0   ){  printf( "SDLNet_Init:        %s\n", SDLNet_GetError());       return -1; }
		if ( !( packet  = SDLNet_AllocPacket( nbuff   ) ) ){  printf( "SDLNet_AllocPacket: %s\n", SDLNet_GetError() );      return -2; }
		if ( !( socket  = SDLNet_UDP_Open(   in_port  ) ) ){  printf( "SDLNet_UDP_Open  in_port: %s\n", SDLNet_GetError()); return -3; }		
		//if ( !( in_socket  = SDLNet_UDP_Open(   in_port  ) ) ){  printf( "SDLNet_UDP_Open  in_port: %s\n", SDLNet_GetError()); return -3; }		
		return 0;
	}

	virtual int tryConnect_UDP( const char *host_str, Uint16 port ){
		printf( "tryConnect_UDP( %s, %i )\n", host_str, port );
		if (             -1 == SDLNet_ResolveHost( &receiver, host_str, port )   ){  printf( "SDLNet_ResolveHost(%s %d): %s\n", host_str, port, SDLNet_GetError()); return -2; }
		connected = true;
		printf( "connected to (%i %i) \n", receiver.host, receiver.port );
		return 0;
	}
  
	virtual void close_UDP(){
		if( packet  != NULL ) SDLNet_FreePacket(packet);
		SDLNet_Quit();
	}

};

#endif
