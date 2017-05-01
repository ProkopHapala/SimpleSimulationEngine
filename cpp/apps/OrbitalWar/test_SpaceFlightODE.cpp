
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Vec3.h"
#include "Mat3.h"

#include "Draw2D.h"
#include "Draw3D.h"

//#include "DynamicOpt.h"
#include "spline_hermite.h"
#include "CubicBSpline.h"
#include "ODEintegrator.h"

#include "SDL_utils.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"


inline void addGravity( const Vec3d& dp, Vec3d& G, double GM ){
    double ir2 = 1.0d/dp.norm2();
    double ir  = sqrt( ir2 );
    double ir3 = ir2*ir;
    G.add_mul( dp, ir3*GM);
}

void getAcceleration( double t, int n, double * Ys, double * dYs ){
    Vec3d pos = *(Vec3d*) Ys;
    Vec3d vel = *(Vec3d*)(Ys+3);

    Vec3d acc; acc.set(0.0d);
    addGravity( pos, acc, -1.0d );

    ((Vec3d*) dYs   )->set( vel );
    ((Vec3d*)(dYs+3))->set( acc );
};


class TestApp_SpaceFlightODE : public AppSDL2OGL_3D {
	public:

    ODEintegrator_RKF45 odeint;

    Vec3d *pos=NULL,*vel=NULL,*acc=NULL;
    Vec3d  opos;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	//virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling( );

	TestApp_SpaceFlightODE( int& id, int WIDTH_, int HEIGHT_ );

};

TestApp_SpaceFlightODE::TestApp_SpaceFlightODE( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    odeint.reallocate( 6 );
    odeint.getDerivs = getAcceleration;
    odeint.dt_max   = 0.5;
    odeint.dt_min   = 0.01;
    odeint.dt_adapt = 0.1;

    ((Vec3d*)(odeint.invMaxYerr  ))->set(1e+8);
    ((Vec3d*)(odeint.invMaxYerr+3))->set(1e+8);


    pos = (Vec3d*)(odeint.Y   );
    vel = (Vec3d*)(odeint.Y+3 );
    acc = (Vec3d*)(odeint.dY+3);

    pos->set(1.0,0.0,0.0); // pos
    vel->set(0.0,1.3,0.0); // vel

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

void TestApp_SpaceFlightODE::draw(){
    //printf( " ==== frame %i \n", frameCount );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    //glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4f(0.5f,0.5f,0.5f, 0.02f);Draw2D::drawRectangle( {-20.0,-20.0},{20.0,20.0}, true);
	//glDisable(GL_BLEND);
	//glClear( GL_DEPTH_BUFFER_BIT );

	int nstep = 0;
	//odeint.step( 0.1d );
	//odeint.adaptive_step_RKF45( );
	nstep = odeint.integrate_adaptive( odeint.dt_adapt, odeint.t+0.2d );
	//odeint.step_RKF45( 0.1d );
	//odeint.step_euler( 0.1d );

	printf( "%i %i (%g,%g,%g) %g %g\n ", frameCount, nstep, pos->x,pos->y,pos->z, odeint.t, odeint.dt_adapt );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine    ( opos      , *pos ); Draw3D::drawPointCross( *pos, 0.1 );
    glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( (*vel)*2.0, *pos );

    opos = *pos;

    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross( {0.0,0.0,0.0}, 0.5 );

	glDisable ( GL_LIGHTING );

};

void TestApp_SpaceFlightODE::drawHUD(){
    glDisable ( GL_LIGHTING );

    //exit(0);
}

// ===================== MAIN

TestApp_SpaceFlightODE * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestApp_SpaceFlightODE( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















