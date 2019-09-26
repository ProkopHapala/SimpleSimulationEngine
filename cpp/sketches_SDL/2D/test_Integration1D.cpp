
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

#include "integration.h"

// ======================  TestApp

class TestAppPlotting : public AppSDL2OGL{
	public:

    Plot2D plot1;
    int fontTex;

	virtual void draw   ();
    //virtual void eventHandling( const SDL_Event& event );

	TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ );
};

TestAppPlotting::TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    Func1d myFunc = &sin;

    plot1.init();
    plot1.fontTex = fontTex;


    DataLine2D * line1 = new DataLine2D(100);
    line1->linspan(-3*M_PI,2*M_PI);
    line1->yfunc = myFunc;

    DataLine2D * line2 = new DataLine2D(400);
    line2->linspan(-M_PI,3*M_PI);
    for(int i=0; i<line2->n; i++){
        double x     = line2->xs[i];
        line2->ys[i] = sin(x*x);
    }
    line2->clr = 0xFF00FF00;

    plot1.lines.push_back( line1 );
    plot1.lines.push_back( line2 );
    plot1.render();

}

void TestAppPlotting::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	//plot1.drawAxes();
    plot1.view();
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
















