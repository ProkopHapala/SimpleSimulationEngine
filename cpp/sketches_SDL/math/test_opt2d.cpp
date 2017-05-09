
#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>

#include "fastmath.h"
#include "functions.h"
#include "warpFunction2D.h"
#include "optimizer_random.h"
#include "SDLplot.h"

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

SDL_Renderer*	render			= NULL;
SDL_Window*		window        	= NULL;
SDL_Surface*	screenSurface 	= NULL;
SDL_Surface*	tempSurf;
SDL_Texture*	tempTex;

SDL_Rect SrcR;
SDL_Rect DestR;

SDL_Event		event;
bool 			STOP          	= false;
int 			frameCount		=	0;

OptimizerRandom* opt1;

double dxMax[] = { 0.05, 0.05 };
double xBest[] = { 2.0, +1.2 };

double tolerance = 0.001;

int nEvals = 0;

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


inline double myFittnessFunc_( double x, double y ){
	//return rosenbrok( x, y );
	//return sinValey( x, y );
	//return cosValey( x, y );
	//return mandelbort( x, y );
	return spiral( x, y );
}

inline double myFittnessFunc( double * xs ){ return -myFittnessFunc_( xs[0], xs[1] ); nEvals++; };

void draw(){
	//nEvals = 0;
	int perFrame=1;
	for (int i=0; i<perFrame; i++){
		double x0 = opt1->xBest[0];
		double y0 = opt1->xBest[1];
		double dfitness = opt1->step();
		double x1 = opt1->xNew[0];
		double y1 = opt1->xNew[1];
		if( dfitness > 0 ){
			SDL_SetRenderDrawColor( render, 0, 128, 0, 255 );
			printf( " %i %f \n", nEvals, opt1->bestFitness  );
			if( -opt1->bestFitness < tolerance){
				printf( " CONVERGENCE ACHIEVED in %i iterations, fittness = %f \n", nEvals, opt1->bestFitness  );
				STOP = true;
				break;
			}
		}else  {
			SDL_SetRenderDrawColor( render, 255, 0, 0, 50  );
		}
		SDL_RenderDrawLine    ( render,  x2i(x0), y2i(y0), x2i(x1), y2i(y1) );
	}
	SDL_RenderPresent( render );
	SDL_UpdateWindowSurface(window);
}

void setup(){

    setZoom( 100.0d );

	window          = SDL_CreateWindow( "SDL Tutorial",   SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
	render        	= SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC );
	SDL_UpdateWindowSurface( window );

	tempSurf        = SDL_CreateRGBSurface(0,SCREEN_WIDTH,SCREEN_HEIGHT,32,0,0,0,0 );
	setPixelsFunction( tempSurf, 0, 0, tempSurf->w,    tempSurf->h, &myFittnessFunc_ );
	//setPixelsFunction( tempSurf, 0, 0, tempSurf->w,    tempSurf->h, -2.0,-2.0, 2.0, 2.0, &mandelbort );

	tempTex   = SDL_CreateTextureFromSurface( render, tempSurf  );

	SrcR.x  = 0; SrcR.y  = 0; SrcR.w  = tempSurf->w; SrcR.h  = tempSurf->h;
	DestR.x = 0; DestR.y = 0; DestR.w = tempSurf->w; DestR.h = tempSurf->h;

	SDL_RenderCopy( render, tempTex, &SrcR, &DestR);

	opt1 = new OptimizerRandom_3( 2, xBest, dxMax, 0.7, 1.2, &myFittnessFunc );
	SDL_SetRenderDrawBlendMode( render, SDL_BLENDMODE_BLEND );
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
