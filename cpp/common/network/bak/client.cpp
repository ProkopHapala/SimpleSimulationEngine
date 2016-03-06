#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "SDL_net.h"
 
int main(int argc, char **argv){

	UDPsocket sd;
	IPaddress srvadd;
	UDPpacket *p;
	int quit;
 
	// Check for parameters
	if (argc < 3){
		fprintf(stderr, "Usage: %s host port\n", argv[0]);
		exit(EXIT_FAILURE);
	}
 
	// Initialize SDL_net
	if (SDLNet_Init() < 0)	{
		fprintf(stderr, "SDLNet_Init: %s\n", SDLNet_GetError());
		exit(EXIT_FAILURE);
	}
 
	// Open a socket on random port 
	if (!(sd = SDLNet_UDP_Open(0)))	{
		fprintf(stderr, "SDLNet_UDP_Open: %s\n", SDLNet_GetError());
		exit(EXIT_FAILURE);
	}
 
	// Resolve server name  
	if (SDLNet_ResolveHost(&srvadd, argv[1], atoi(argv[2])) == -1)	{
		fprintf(stderr, "SDLNet_ResolveHost(%s %d): %s\n", argv[1], atoi(argv[2]), SDLNet_GetError());
		exit(EXIT_FAILURE);
	}
 	printf( "connected to (%i %i) \n", srvadd.host, srvadd.port );

	// Allocate memory for the packet 
	if (!(p = SDLNet_AllocPacket(512)))	{
		fprintf(stderr, "SDLNet_AllocPacket: %s\n", SDLNet_GetError());
		exit(EXIT_FAILURE);
	}
 
	// Main loop
	quit = 0;
	while (!quit)	{
		printf("Fill the buffer\n>");
		scanf("%s", (char *)p->data);
 
		p->address.host = srvadd.host;	// Set the destination host
		p->address.port = srvadd.port;	// And destination port
 
		p->len = strlen((char *)p->data) + 1;
		SDLNet_UDP_Send(sd, -1, p);     // This sets the p->channel
 
		// Quit if packet contains "quit" 
		if (!strcmp((char *)p->data, "quit")) quit = 1;
	}
 
	SDLNet_FreePacket(p);
	SDLNet_Quit();
 
	return EXIT_SUCCESS;
}
