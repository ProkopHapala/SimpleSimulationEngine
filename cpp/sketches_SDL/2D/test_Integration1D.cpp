/// @file @brief This program demonstrates various numerical integration methods for 1D functions. It visualizes the function and the approximation process (e.g., Riemann sums, trapezoidal rule), allowing the user to compare the accuracy and behavior of different integration techniques.

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

const int npar = 10;
Vec3d params[npar];

int neval = 0;
bool bDrawSamples = false;

double testFunc(double x){
    neval++;
    double sum=0.5;
    for(int i=0;i<npar;i++){
        const Vec3d& par = params[i];
        double dx = (x-par.x)/par.y;
        //sum+= par.z/(1+dx*dx);
        sum+= par.z*exp( -sqrt( 0.01 + dx*dx) );
    }
    //glVertex3f(x,sum,0);
    //Draw3D::drawPointCross();
    if(bDrawSamples){ Draw2D::drawLine({x,0},{x,sum});  }
    return sum;
}

class TestAppPlotting : public AppSDL2OGL{
	public:

    Plot2D plot1;
    //int fontTex;

    GLint ogl_samples=0;

	virtual void draw();
    //virtual void eventHandling( const SDL_Event& event );

	TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ );
};










TestAppPlotting::TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    //Func1d myFunc = &sin;
    Func1d myFunc = &testFunc;

    for(int i=0;i<npar;i++){
        params[i].x=randf(-1.,1.);
        params[i].y=randf(0.1,0.3);
        params[i].z=randf(-1.,1.);
    }


    plot1.init();
    plot1.fontTex = fontTex;

    srand(545454);
    DataLine2D * line1 = new DataLine2D(500);
    line1->linspan(-3,3);
    line1->yfunc = myFunc;

    /*
    DataLine2D * line2 = new DataLine2D(400);
    line2->linspan(-M_PI,3*M_PI);
    for(int i=0; i<line2->n; i++){
        double x     = line2->xs[i];
        line2->ys[i] = sin(x*x);
    }
    line2->clr = 0xFF00FF00;
    */

    plot1.lines.push_back( line1 );
    //plot1.lines.push_back( line2 );
    plot1.render();

    bDrawSamples = true;
    ogl_samples = glGenLists(1);
    glNewList(ogl_samples, GL_COMPILE);

    double eps = 1.0;
    for(int i=0; i<7; i++){

        glColor3f(0.0,0.0,1.0);
        neval = 0;
        double I1 = AdaptiveIntegral::trapez ( &testFunc, -3., 3., eps, 6 );
        int n1=neval;

        glColor3f(1.0,0.0,0.0);
        neval = 0;
        double I2 = AdaptiveIntegral::simpson( &testFunc, -3., 3., eps, 6 );
        int n2=neval;

        printf( "%i %g trapez %g %i simpson %g %i \n", i, eps, I1, n1, I2, n2 );
        eps*=0.1;
    }

    glEndList();
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
