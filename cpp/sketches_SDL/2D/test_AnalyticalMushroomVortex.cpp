
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "Vec2.h"
#include "Plot2D.h"

#include "PotentialFlow.h"


/*

Analytical expression:
Cornu Spiral


https://en.wikipedia.org/wiki/Fresnel_integral#Euler_spiral

Doublet Vortex:
https://nbviewer.jupyter.org/github/barbagroup/AeroPython/blob/master/lessons/03_Lesson03_doublet.ipynb

https://www.ecourses.ou.edu/cgi-bin/ebook.cgi?doc=&topic=fl&chap_sec=07.4&page=theory
https://link.springer.com/chapter/10.1007/978-981-13-1678-4_4


*/


void flowField(Vec2d p, Vec2d& v){
    //v.x = sin(p.y);// + cos(p.x);
    //v.y = sin(p.x);// - sin(p.x);

    v.x =  0.0;
    v.y = -1.0;

    /*
    Vec2d d  = p-(Vec2d){0.0,0.1};
    double r = d.norm();
    v.add_mul( d, 10./(r*r*r) ); // source

    d = p-(Vec2d){0.0,-0.1};
    r = d.norm();
    v.add_mul( d, -10./(r*r*r) ); // source
    */

    Vec3d v3 = sourceDipol( {p.x,p.y,0.0}, (Quat4d){0.0,1.0,0.0,0.0} ) * 2.0;
    v.x+=v3.x;
    v.y+=v3.y;


    double sz = v.norm();
    if(sz>1.0) v.mul(1./sz); ;


}

Vec2d propagate(Vec2d p, int nt, double dt, double* xs, double* ys){
    Vec2d v;
    Vec2d rot=(Vec2d){-0.70710678118,0.70710678118};
    for(int i=0; i< nt; i++){

        //p.mul_cmplx(rot);

        flowField( p, v );
        //v.udiv_cmplx(rot);
        //p.udiv_cmplx(rot);

        p.add_mul( v, dt );
        if(xs)xs[i] = p.x;
        if(ys)ys[i] = p.y;
    }
    return p;
}


class AnalyticVortexApp : public AppSDL2OGL {
	public:
	int frameCount = 0;
	bool loopEnd,STOP;

    static const int np = 1024;
	Vec2d particles[np];

	Vec2d pmin={-3.,0.0};
	Vec2d pmax={+3.,5.5};

	QuePlot2D trj;
	Plot2D plot1;

	//int nt;
	//double xs[nt];
	//double ys[nt];

	// ---- function declarations
	//virtual void loop( int nframes );
	virtual void draw();
	//virtual void drawHUD();
	AnalyticVortexApp( int& id, int WIDTH_, int HEIGHT_ );

};


AnalyticVortexApp::AnalyticVortexApp( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    trj.init( 256 , 2 );

    plot1.init();
    plot1.xsharingLines( 2, 1024 );


    double L = M_PI*M_SQRT2;
    //Vec2d dp=(Vec2d){l/np,0.0};
    for(int i=0; i<np; i++){
        //particles[i]={-0.5*i*L/np,-L/2};
        //particles[i]=(Vec2d){randf(-M_PI*M_SQRT2,0),randf(-M_PI*M_SQRT2,0)};
        particles[i] = (Vec2d){randf(pmin.x,pmax.x),randf(pmin.y,pmax.y)};
    }

    //propagate( particles[np/2], plot1.lines[0]->n, 0.1, plot1.lines[0]->ys, plot1.lines[1]->ys );

    plot1.render();

    //.zoom = 10;

}

void AnalyticVortexApp::draw(){
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(1.,1.,1., 0.02);
    Draw2D::drawRectangle( {-100.,-100.}, {100.,100.}, true );
    glDisable(GL_BLEND);


    //glClearColor( 0.9f, 0.9f, 0.9f, 0.1f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //if(frameCount==2)zoom=1.0;

    for(int i=0;i<5; i++){
        int ip = rand()%np;
        particles[ip] = (Vec2d){randf(pmin.x,pmax.x),randf(pmin.y,pmax.y)};
    }

    glColor4f(1.,.0,.0,1.);
    Draw2D::drawLine( {0.0,0.0}, {M_PI*0.5,0.0} );


    int ipick = 0;
    glColor4f(0.,.0,.0,1.);
    glBegin(GL_POINTS);
    for(int i=0; i<np; i++){
        particles[i] = propagate( particles[i], 1, 0.02,  0,0 );
        //printf( "[%i] %g %g \n", i, particles[i].x, particles[i].y );
        glVertex2f( particles[i].x, particles[i].y );

        if(i==ipick){
            trj.next( frameCount );
            trj.set_back(0, particles[i].x );
            trj.set_back(1, particles[i].y );
        }
    }
    glEnd();

    glColor3f(1.,.0,.0);
    //Draw2D::drawRectangle({-M_PI*M_SQRT2,-M_PI*M_SQRT2},{0,0},false);


    //glColor3f(1.,1.,.0);
    //Draw2D::drawPlot2D( trj.n, trj.data[0], trj.data[1],(Vec2d){1.,1.}, (Vec2d){0.,0.} );


    particles[0] = (Vec2d){mouse_begin_x,mouse_begin_y};
    propagate( particles[0], plot1.lines[0]->n, 0.1, plot1.lines[0]->ys, plot1.lines[1]->ys );
    Draw2D::plot( plot1.lines[0]->n, plot1.lines[0]->ys, plot1.lines[1]->ys );

    //plot1.render();
    //plot1.view();

    /*
	glColor3f(0.0,.0,.0);
	for( int ix = 0; ix<5; ix++ ){

        double xoff = ix * 0.1;
        glBegin(GL_LINE_STRIP);
        double dt = 0.2;
        for(int i=0; i<256; i++){
            double t   = i*dt;
            //double t2  = t*t;
            double lor = 1/(1+t);
            double x = cos(t)*lor*(1+xoff);
            double y = sin(t)*lor*(1+xoff);
            glVertex2f( x, y );
        }
        glEnd();
	}
	*/

};

// ===================== MAIN

AnalyticVortexApp * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new AnalyticVortexApp( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















