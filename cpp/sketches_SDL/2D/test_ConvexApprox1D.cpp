
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "SDL_utils.h"

#include "Plot2D.h"
#include "PlotScreen2D.h"

//#include "Noise.h"
//#include "Grid2DAlgs.h"
#include "dataprocess1D.h"


int np;
double * data;



//#include "integration.h"

// ======================  TestApp



class TestAppPlotting : public AppSDL2OGL{   public:

    Plot2D plot1;
    //int fontTex;

    GLint ogl_samples=0;

	virtual void draw();
    //virtual void eventHandling( const SDL_Event& event );

	TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ );
};










TestAppPlotting::TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    plot1.init();
    plot1.fontTex = fontTex;

    srand(545454);
    int npow = 8;
    int n    = 1<<npow;
    DataLine2D * line1 = new DataLine2D(n); plot1.lines.push_back( line1 ); line1->clr = 0xFFFF00FF; line1->linspan(-3,3);
    DataLine2D * line2 = new DataLine2D(n); plot1.lines.push_back( line2 ); line2->clr = 0xFFFF0000; line2->linspan(-3,3);
    DataLine2D * line3 = new DataLine2D(n); plot1.lines.push_back( line3 ); line3->clr = 0xFF008000; line3->linspan(-3,3);

    bisecNoise1D( npow, line1->ys, -0.02, 0.02 );
    //runningMax( n, 8, line1->ys, line2->ys);

    double slope = -0.05;
    slopeSweep    ( n, slope, line1->ys, line2->ys );
    slopeSweepBack( n, slope, line1->ys, line3->ys );

    /*
    int     is[n];
    double his[n];
    int ni = convexApprox( n, 16, line1->ys, is, his);
    DataLine2D * line2 = new DataLine2D(ni); plot1.lines.push_back( line2 ); line2->clr = 0xFF008000; line2->linspan(-3,3);
    for(int i=0; i<ni; i++){
        printf( "[%i]  i %i x %g y %g \n", i, is[i], line1->xs[is[i]], his[i] );
        line2->ys[i]=line1->xs[is[i]]; line2->ys[i]=his[i];
    }
    */
    //plot1.lines.push_back( line2 );
    //plot1.lines.push_back( line3 );
    //plot1.lines.push_back( line2 );
    plot1.render();

    //bDrawSamples = true;
    //ogl_samples = glGenLists(1);
    //glNewList(ogl_samples, GL_COMPILE);
    //glEndList();
    //printf( "I(testFunc,a,b) = %g | neval %i \n", I, neval );

}

void TestAppPlotting::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	//plot1.drawAxes();
    plot1.view();

    glCallList(ogl_samples);
};

// ===================== MAIN

TestAppPlotting * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppPlotting( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















