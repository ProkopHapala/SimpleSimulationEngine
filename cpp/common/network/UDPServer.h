#ifndef  UDPServer_h
#define  UDPServer_h

#include "SDL_net.h"


/*

TODO :

    for Server

    There will be multiple clients for one server; therefore server should has following structure:

    int       nclinents;
    IPaddress receivers[];
    UDPsocket sockets  [];

    - when client is connected virtual function onConnect() is activated;
    - server give the client "net_id" according to order of connection, it is also index in array of clients


    functions:
    ---------
    onSend()
    onReceive()
    onConnect()
    onDisconnect()  ?( at begining not required         )
    onTimeOut()     ?( can be used instead onDisconnect )

*/


class UDPServer {
    public:
	UDPsocket   socket;
	UDPpacket  *packet = NULL;

	int nMaxClients = 0;
	int nClients    = 0;
	IPaddress * clients;
	bool      * connected = false;

    virtual bool onSend   ( int iClient       );
    virtual void onRecieve( int iClient       );
    virtual int  onConnect( IPaddress client );

    //virtual void onDisconnect( int iClient );
    //virtual void onTimeOut   ( int iClient );

    int  findKnownClient ( IPaddress client );
    void toClients       ( );
    void fromClients     ( );
	void printPacketInfo ( );

    virtual int  init_UDP ( Uint16 in_port, Uint16 nbuff, int nMaxClinets_ );
	virtual void close_UDP();

};

/// ===============================
///    Implementation
/// ===============================

bool UDPServer::onSend( int iClient ){
    printf( "UDPServer::onSend( %i )\n", iClient );
    return true;
}

void UDPServer::onRecieve( int iClient ){
    printf( "UDPServer::onRecieve( %i )\n", iClient );
    printPacketInfo( );
}

int  UDPServer::onConnect( IPaddress client ){
    int iClient = nClients;
    if( iClient < nMaxClients ){
        clients[iClient] = client;
        nClients++;
        printf( " connected new client %i \n", iClient );
        return iClient;
    }else{
        printf( " server is full on full capacity ( %i clients ) \n", nMaxClients  );
        return -1;
    }
}

void UDPServer::toClients( ){
    //printf( "UDPServer::toClients %i nClients \n", nClients );
    for( int i=0; i<nClients; i++){
        //printf( "clients %i begin \n", i );
        onSend( i );
        packet->address.host = clients[ i ].host;
        packet->address.port = clients[ i ].port;
        SDLNet_UDP_Send( socket, -1, packet );
        //printf( "clients %i end \n", i );
    }
    //printf( "UDPServer::toClients DONE \n", nClients );
}

void UDPServer::fromClients( ){
    while ( SDLNet_UDP_Recv( socket, packet ) )	{
        //printPacketInfo();
        int iClient = findKnownClient ( packet->address );
        if ( iClient >= 0 ){
            onRecieve( iClient );
        }else{
            onConnect( packet->address );
            onRecieve( iClient );
        }
    }
}

int  UDPServer::findKnownClient ( IPaddress client ){
    for( int i=0; i<nClients; i++ ){
        if( client.host == clients[ i ].host ){
            if( client.port == clients[ i ].port ){
                return i;
            }
        }
    }
    return -1;
}

void UDPServer::printPacketInfo( ){
    printf("UDP Packet incoming\n");
    printf("\tChan:    %d\n", packet->channel);
    printf("\tData:    %s\n", (char *)packet->data);
    printf("\tLen:     %d\n", packet->len);
    printf("\tMaxlen:  %d\n", packet->maxlen);
    printf("\tStatus:  %d\n", packet->status);
    printf("\tAddress: %x %x\n", packet->address.host, packet->address.port);
}

int UDPServer::init_UDP( Uint16 in_port, Uint16 nbuff, int nMaxClients_ ){
    nMaxClients = nMaxClients_;
    clients = new IPaddress[nMaxClients];
    if (              SDLNet_Init()             < 0   ){  printf( "SDLNet_Init:        %s\n", SDLNet_GetError());       return -1; }
    if ( !( packet  = SDLNet_AllocPacket( nbuff   ) ) ){  printf( "SDLNet_AllocPacket: %s\n", SDLNet_GetError() );      return -2; }
    if ( !( socket  = SDLNet_UDP_Open(   in_port  ) ) ){  printf( "SDLNet_UDP_Open  in_port: %s\n", SDLNet_GetError()); return -3; }
    return 0;
}

void UDPServer::close_UDP(){
    if( packet  != NULL ) SDLNet_FreePacket(packet);
    SDLNet_Quit();
    delete clients;
}

#endif
