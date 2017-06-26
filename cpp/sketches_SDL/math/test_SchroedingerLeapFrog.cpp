
#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>

#include "fastmath.h"

#include "SDLplot.h"

#include "SchroedingerLeapFrog.h"

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

double dx  = 0.01;

double x0  = 0.0;
double Y0  = 0.0;
double dY0 = 1.0;

double dE = 0.01;
double E  = 4.0;

double dfcomb = 0.01;
double fcomb  = 0.5;

const int nbuff = 2000;

double  V_buff[nbuff];
double  x_buff[nbuff];

double  y_buff [nbuff];
double  yA_buff[nbuff];
double  yB_buff[nbuff];
double  dy_buff [nbuff];
double  dyA_buff[nbuff];
double  dyB_buff[nbuff];

const int nE = 1000;
double E_buff     [nbuff];
double y2out_buff [nbuff];

// ===============================
// ===== GLOBAL VARIABLES PLOTTING
// ===============================


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

void draw(){

	double fb_best=0, y2out=0;
	evalResidua( nbuff, nbuff/2, dx, E, V_buff, fb_best, y2out );
	double fa_best= 1.0 - fb_best;
	integrate_buff ( nbuff, nbuff/2, dx, fa_best, fb_best, E, V_buff, y_buff, dy_buff );

	double densIn = dot( nbuff/4, 3*nbuff/4, y_buff, y_buff )*dx;

	printf( " E % fA = %f fB = %f y2out = %f densIn \n", E, fa_best, fb_best, y2out, densIn );

	double Yscale  = 0.5/sqrt(densIn);
	double dYscale = 0.5;
	double Vscale  = 0.2;
	SDL_SetRenderDrawColor( plot.render,  255, 255, 255, 255 );    SDL_RenderClear( plot.render);
	SDL_SetRenderDrawColor( plot.render, 128, 128, 128, 255 );     plot.plotAxes ( );
	SDL_SetRenderDrawColor( plot.render,  50, 128,  50, 255 );     plot.plotHline( E*Vscale );
	SDL_SetRenderDrawColor( plot.render,   0,   0,   0, 255 );     plot.plotFunc ( nbuff, x_buff,  V_buff,   Vscale  );
	SDL_SetRenderDrawColor( plot.render, 255, 128, 255,   255 );   plot.plotFunc ( nbuff, x_buff,  y_buff, Yscale  );
	SDL_SetRenderDrawColor( plot.render, 100, 200, 100,   255 );   plot.plotFunc ( nE, y2out_buff, E_buff, Vscale  );

	SDL_RenderPresent( plot.render );
	SDL_UpdateWindowSurface( plot.window);

	STOP = true;
}



void setup(){
	//plot.window         = SDL_CreateWindow( "SDL Tutorial",   SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN );
	//plot.render        	= SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC );

	plot.init( "Schoedinger Leap Frog 1D", 1024, 512 );

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
