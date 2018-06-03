

#include <SDL.h>
#include <math.h>
#include <stdio.h>


const  float INV_RAND_MAX = 1.0f/RAND_MAX;
inline float randf(){ return INV_RAND_MAX*rand(); }



// ===============================
// ===== GLOBAL VARIABLES
// ===============================

#include "integrator.h"

double dx  = 0.01;

double x0  = 0.0;
double Y0  = 0.0;
double dY0 = 1.0;

double dE = 0.01;
double E  = 4.0;

double dfcomb = 0.01;
double fcomb  = 0.5;

const int nbuff = 1000;

double  V_buff[nbuff];
double  x_buff[nbuff];

double  y_buff [nbuff];
double  yA_buff[nbuff];
double  yB_buff[nbuff];
double dy_buff [nbuff];
double dyA_buff[nbuff];
double dyB_buff[nbuff];

const int nE = 1000;
double E_buff     [nbuff];
double y2out_buff [nbuff];



const int SCREEN_WIDTH  = 1024;
const int SCREEN_HEIGHT = 512;

const int IX0 = SCREEN_WIDTH/2;
const int IY0 = SCREEN_HEIGHT/2;

const double ZOOM  = 100.0;
const double IZOOM = 1.0d/ZOOM;

double i2x( int ix ){ return (ix - IX0)*IZOOM;  };
double i2y( int iy ){ return (iy - IY0)*IZOOM;  };

int x2i( double  x ){ return x*ZOOM + IX0;     };
int y2i( double  y ){ return y*ZOOM + IY0;     };

#include "../../../common/SDL2/plotFunc.h"

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


/*

 ( (1-a)*yAL + a*yBL )^2   +     ( (1-a)*yAR + a*yBR )^2   minimize
 ( yAL + a*(yBL-yAL) )^2   +     ( yAR + a*(yBR-yAR) )^2   minimize
 ( yAL + a*yABL      )^2   +     ( yAR + a*yABR      )^2   minimize

  yAL^2 + 2*a*yABL + (a*yABL)^2   +    yAL^2 + 2*a*yABR + (a*yABR)^2

  2*yABL   +   a*yABL^2            2*yABR   +   a*yABR^2 = 0 

  a * ( yABL^2 + yABR^2  )   =  -2 *( yABL + yABR )

*/


void draw(){

/*
	double yAL=0,yAR=0,yBL=0,yBR=0;
	integrate     ( nbuff, nbuff/2, dx, 1.0, 0.0, E, V_buff, yAL, yAR );
	integrate     ( nbuff, nbuff/2, dx, 0.0, 1.0, E, V_buff, yBL, yBR );

	double yABL    = yBL - yAL;
	double yABR    = yBR - yAR;
	double fb_best  =  -( yABL*yAL + yABR*yAR ) / ( yABL*yABL + yABR*yABR );
	double fa_best  =  1.0-fa_best;

	double yL      =  fa_best*yAL + fb_best*yBL;
	double yR      =  fa_best*yAR + fb_best*yBR;
	double y2out   =  yL*yL + yR*yR;
*/

	double fb_best=0, y2out=0;
	evalResidua( nbuff, nbuff/2, dx, E, V_buff, fb_best, y2out );
	double fa_best= 1.0 - fb_best;
	integrate_buff ( nbuff, nbuff/2, dx, fa_best, fb_best, E, V_buff, y_buff, dy_buff );

	double densIn = dot( nbuff/4, 3*nbuff/4, y_buff, y_buff )*dx;

	printf( " E % fA = %f fB = %f y2out = %f densIn \n", E, fa_best, fb_best, y2out, densIn );

	double Yscale  = 0.5/sqrt(densIn);
	double dYscale = 0.5;
	double Vscale  = 0.2;
	SDL_SetRenderDrawColor(render,  255, 255, 255, 255 );   SDL_RenderClear(render);
	SDL_SetRenderDrawColor( render, 128, 128, 128, 255 );   plotAxes( render );
	SDL_SetRenderDrawColor( render,  50, 128,  50, 255 );   plotHline( render, E*Vscale );
	SDL_SetRenderDrawColor( render, 0, 0, 0, 255 );         plotFunc( render, nbuff, x_buff,  V_buff,   Vscale  );
	SDL_SetRenderDrawColor( render, 255, 128, 255,   255 );   plotFunc( render, nbuff, x_buff,  y_buff, Yscale  );
	SDL_SetRenderDrawColor( render, 100, 200, 100,   255 );   plotFunc( render, nE, y2out_buff, E_buff, Vscale  );

	SDL_RenderPresent( render );
	SDL_UpdateWindowSurface(window);
	
	STOP = true;
}



void setup(){
	window          = SDL_CreateWindow( "SDL Tutorial",   SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
	render        	= SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC );
	SDL_UpdateWindowSurface( window );

	eval_potential( nbuff, dx, -5.0, x_buff, V_buff );

	scanResidua( nbuff, nbuff/2, dx,  nE,  0.0,  15.0, V_buff, E_buff, y2out_buff );

}

void inputHanding(){
	while(SDL_PollEvent(&event)){ 
		if( event.type == SDL_KEYDOWN ){ 
			if(event.key.keysym.sym == SDLK_ESCAPE) { quit(); }
			if(event.key.keysym.sym == SDLK_SPACE    ){ STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); }
			//if(event.key.keysym.sym == SDLK_PLUS     ){ E += dE; printf( " E = %f \n", E ); }
			//if(event.key.keysym.sym == SDLK_MINUS    ){ E -= dE; printf( " E = %f \n", E ); }
			if(event.key.keysym.sym == SDLK_PAGEUP     ){ E += dE; /*printf( " E = %f \n", E );*/	STOP = false; }
			if(event.key.keysym.sym == SDLK_PAGEDOWN   ){ E -= dE; /*printf( " E = %f \n", E );*/ STOP = false; }
			if(event.key.keysym.sym == SDLK_HOME       ){ fcomb += dfcomb; printf( " fA = %f fB = %f \n", 1.0-fcomb, fcomb );	STOP = false; }
			if(event.key.keysym.sym == SDLK_END        ){ fcomb -= dfcomb; printf( " fA = %f fB = %f \n", 1.0-fcomb, fcomb ); 	STOP = false; }
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
