
#include <SDL2/SDL.h>
#include <math.h>
#include <stdio.h>

#include "fastmath.h"
#include "functions.h"
#include "warpFunction2D.h"
#include "optimizer_random.h"
#include "SDLplot.h"

#include "TestFuncND.h"

// ===============================
// ===== GLOBAL VARIABLES
// ===============================

SDLplot plot;
SDL_Event		event;
bool 			STOP          	= false;
int 			frameCount		=	0;

OptimizerRandom* opt1;

double dxMax[] = { 0.05, 0.05 };
double xBest[] = { 2.0, +1.2 };
double freqScales[] = { 2.0, 1.0, 0.7, 0.4 };
double tolerance = 0.001;
int nEvals = 0;
TestFuncND funcND;

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
    Vec2d xs; xs.set(x,y);
    double f =  funcND.eval( ((double*)&xs) ); return 1-1/(1+f*VALEY_TIGHTNESS) + fabs(x+2)*0.1;
    //return funcND.eval( ((double*)&xs) );
	//return rosenbrok( x, y );
	//return sinValey( x, y ) + fabs(x+2)*0.1;
	//return spiral( x, y ) + sqrt(harmonic(x,y))*0.1 - lorenz(x*30,y*30);

    //return warped_function_noDiff( x, y );
	//return cross_valey(x,y) + harmonic(x,y)*0.1;
	//return mandelbort( x, y );
	//return particles(x,y) + harmonic(x,y)*0.1;
}

inline double myFittnessFunc( double * xs ){ nEvals++; return -myFittnessFunc_( xs[0], xs[1] );  };

int perform_relaxation( int nMaxIter ){

	for (int i=0; i<nMaxIter; i++){
		double x0 = opt1->xBest[0];
		double y0 = opt1->xBest[1];
		double dfitness = opt1->step();
		double x1 = opt1->xNew[0];
		double y1 = opt1->xNew[1];
		if( dfitness > 0 ){
			SDL_SetRenderDrawColor( plot.render, 0, 128, 0, 255 );
			if( -opt1->bestFitness < tolerance){
				printf( " CONVERGENCE ACHIEVED in %i iterations, fittness = %f \n", nEvals, opt1->bestFitness  );
				//STOP = true;
				return i;
			}
		}else  {
			SDL_SetRenderDrawColor( plot.render, 255, 0, 0, 50  );
		}
		//SDL_RenderDrawLine    ( render,  x2i(x0), y2i(y0), x2i(x1), y2i(y1) );
		plot.drawLine( x0, y0, x1, y1 );
	}
    return -1;
}

void restartRelaxation(){
    nEvals = 0;
    double t = randf( 2.0,4.0);
	opt1->xBest[0]= t;
	funcND.getNodeAt( t, opt1->xBest );
    //opt1->xBest[0]=randf( 2.0,4.0);
	//opt1->xBest[1]=randf(-3.0,3.0);
    opt1->restart();
}

void newTestFunc(){
    printf("generationg new test func ... \n");
    funcND.setRandomStiffness(0.5d,2.0d);
    funcND.setRandomCoefs(freqScales, 1.0, 0.0 );

    restartRelaxation();

    plot.setPixelsFunctionClamped( plot.tempSurf, 0, 0, plot.tempSurf->w,  plot.tempSurf->h, &myFittnessFunc_, 0.0, 1.0 );
    //SDL_DestroyTexture(tempTex);
    //tempTex   = SDL_CreateTextureFromSurface( render, tempSurf  );
    //SDL_RenderCopy( render, tempTex, &SrcR, &DestR );
    plot.updateSurfTex();
    printf(" ... DONE \n");
}

void draw(){
    /*
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
    */
	restartRelaxation();
	perform_relaxation( 10000 );

	SDL_RenderPresent(       plot.render );
	SDL_UpdateWindowSurface( plot.window);

}

void setup(){

    NWRAP = 2;

    funcND.allocate(2, 4);
    //funcND.setRandomStiffness(0.5d,2.0d);
    //funcND.setRandomCoefs(freqScales);

    //double f = myFittnessFunc_( 0.5, 0.4 );
    //printf( "(%f,%f) -> %f",  f);
    //exit(0);

    opt1 = new OptimizerRandom_2( 2, xBest, dxMax, 0.85, 1.2, &myFittnessFunc );
	//opt1 = new OptimizerRandom_3( 2, xBest, dxMax, 0.7, 1.2, &myFittnessFunc );

    plot.setZoom( 50.0d );

	plot.init( "Schoedinger Leap Frog 1D", 512, 512 );
	newTestFunc();

	SDL_RenderCopy( plot.render, plot.tempTex, &plot.SrcR, &plot.DestR);
	SDL_SetRenderDrawBlendMode( plot.render, SDL_BLENDMODE_BLEND );
}

void inputHanding(){
	while(SDL_PollEvent(&event)){
		if( event.type == SDL_KEYDOWN ){
			if(event.key.keysym.sym == SDLK_ESCAPE   ){ quit(); }
			if(event.key.keysym.sym == SDLK_SPACE    ){ STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); }
			if(event.key.keysym.sym == SDLK_n        ){ newTestFunc(); }
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
