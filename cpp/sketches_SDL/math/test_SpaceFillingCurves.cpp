
#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>

#include "spaceFillingCurves.h"
#include "pixMapFractal.h"
#include "SDLplot.h"

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

/*
SDL_Renderer*	render			= NULL;
SDL_Window*		window        	= NULL;
SDL_Surface*	screenSurface 	= NULL;
SDL_Surface*	tempSurf;
SDL_Texture*	tempTex;

SDL_Rect SrcR;
SDL_Rect DestR;
*/

SDLplot plot;

SDL_Event		event;
bool 			STOP          	= false;
int 			frameCount		=	0;

// ====================================
// ===== FUNCTION FORWARD DECLARATIONS
// ====================================

void quit(){SDL_Quit(); exit(1);};
void setup();
void draw();
void inputHanding();

// ===============================
// ===== FUNCTION IMPLEMENTATION
// ===============================

int pixelFunc(int ix, int iy){
    //return binaryPixMapFrac( 256, 2, ix, iy, patterns_1 )<<15;
    return HiblertCurve::xy2d( 65536, ix, iy);
    //return ix^iy;
}

void draw(){
	SDL_RenderPresent( plot.render );
	SDL_UpdateWindowSurface( plot.window);
}

void setup(){

    plot.setZoom( 100.0d );
    plot.init( "space filling curve", 512, 512 );

	//window          = SDL_CreateWindow( "SDL Tutorial",   SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
	//render        	= SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC );
	//SDL_UpdateWindowSurface( window );

	//tempSurf        = SDL_CreateRGBSurface(0,SCREEN_WIDTH,SCREEN_HEIGHT,32,0,0,0,0 );

    //pixelFunc(1554,6548);

	plot.setPixelsFunc2i( plot.tempSurf, 0, 0, plot.tempSurf->w,    plot.tempSurf->h, &pixelFunc );

	plot.tempTex   = SDL_CreateTextureFromSurface( plot.render, plot.tempSurf  );
	//SrcR.x  = 0; SrcR.y  = 0; SrcR.w  = tempSurf->w; SrcR.h  = tempSurf->h;
	//DestR.x = 0; DestR.y = 0; DestR.w = tempSurf->w; DestR.h = tempSurf->h;
	SDL_RenderCopy( plot.render, plot.tempTex, &plot.SrcR, &plot.DestR);


	SDL_SetRenderDrawBlendMode( plot.render, SDL_BLENDMODE_BLEND );
}

void inputHanding(){
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE    ){ STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); }
		}
		if( event.type == SDL_QUIT){ quit();  };
	}
}

int main( int argc, char* args[] ){
	setup();
	for( frameCount=0; frameCount<1000000; frameCount++ ){
		if (!STOP){ draw();  }
		inputHanding();
		SDL_Delay( 100 );
	}
	return 0;
}
