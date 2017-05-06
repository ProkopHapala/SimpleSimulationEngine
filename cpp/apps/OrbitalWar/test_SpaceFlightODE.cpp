
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

// ==== Globals

constexpr int nCP  = 10;
double  thetaCPs [nCP];
double  thrustCPs[nCP];
Vec3d   dirCPs   [nCP];

double  uT    = 50.0; // duration of each CP interval
double inv_uT = 1.0/uT;

// rocket parameters - Saturn V
double  vexhaust      = 4500.0; // [km/s]
double  dm_F          = 1.0/vexhaust;
double  mass_initial  = 2.970e+6; // put there zero to switch off the engine
double  mass_empty    = 130e+3 + 140e+3 + 40e+3 + 13.5e+3;
double  thrust_full   = 35e+6;

// planet parameters
const double planet_R   = 6371e+3;
double       planet_GM  = 3.9860044189e+14; //https://en.wikipedia.org/wiki/Standard_gravitational_parameter
Vec3d        planet_pos = {0.0,-planet_R,0.0};



double view_scale = 1e-6;

// ==== Functions

inline double getSpline( double u, double* CPs ){
    double cp0 = CPs[0];
    double cp1 = CPs[1];
    double cp2 = CPs[2];
    double cp3 = CPs[3];
    return Spline_Hermite::val( u, cp1, cp2, (cp2-cp0)*0.5d, (cp3-cp1)*0.5d );
}

inline Vec3d getSpline3d( double u, Vec3d* CPs ){
    //double b0,b1,b2,b3;
    //basis( TYPE x, TYPE& c0, TYPE& c1, TYPE& d0, TYPE& d1 );
    Vec3d cp0 = CPs[0];
    Vec3d cp1 = CPs[1];
    Vec3d cp2 = CPs[2];
    Vec3d cp3 = CPs[3];
    return { // TODO Optimize this
        Spline_Hermite::val( u, cp1.x, cp2.x, (cp2.x-cp0.x)*0.5d, (cp3.x-cp1.x)*0.5d ),
        Spline_Hermite::val( u, cp1.y, cp2.y, (cp2.y-cp0.y)*0.5d, (cp3.y-cp1.y)*0.5d ),
        Spline_Hermite::val( u, cp1.z, cp2.z, (cp2.z-cp0.z)*0.5d, (cp3.z-cp1.z)*0.5d )
    };
}

void splinePixel( double t, double& x, double& y, double& z ){
    int    i = (int) t;
    double u = t - i;
    /*
    double y1 = thetaCPs[i+1];
    double y2 = thetaCPs[i+2];
    x = t+1.0;
    y = Spline_Hermite::val( u, y1, y2, (y2-thetaCPs[i])*0.5d, (thetaCPs[i+3]-y1)*0.5d ) * 5.0;
    z = 0;
    */
    Vec3d p = getSpline3d( u, dirCPs+i );
    x=p.x; y=p.y; z=p.z;
};

inline void addGravity( const Vec3d& dp, Vec3d& G, double GM ){
    double ir2 = 1.0d/dp.norm2();
    double ir  = sqrt( ir2 );
    double ir3 = ir2*ir;
    G.add_mul( dp, -ir3*GM);
}

void getAcceleration( double t, int n, double * Ys, double * dYs ){
    Vec3d pos = *(Vec3d*) Ys;
    Vec3d vel = *(Vec3d*)(Ys+3);
    double m  = Ys[6];

    Vec3d acc; acc.set(0.0d);
    addGravity( pos-planet_pos, acc, planet_GM );

    // engine thrust
    double dm = 0.0d;
    if( m > mass_empty ){   //  check if propelant exhausted
        // --- interpolate from inputs
        double tcp = t*inv_uT;
        int    icp = (int) tcp;
        if((icp+3)<nCP){
            double u   = tcp - icp;

            // TODO : this can be very much optimized if we evaluate spline basis instead !!!!
            Vec3d  dir    = getSpline3d( u, dirCPs   +icp );
            double thrust = getSpline  ( u, thrustCPs+icp );

            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawVecInPos( dir*5.0, pos*view_scale );
            glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( pos*view_scale, pos*view_scale + dir );
            //Vec3d dir_ = {1.0,dir.y,0.0};
            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( pos*view_scale, pos*view_scale + dir_ );
            //printf( "t %g tcp %g icp %i u %g (%g,%g,%g) \n", t, tcp, icp, u, dir.x,dir.y,dir.z  );

            dm =     -dm_F * thrust; // propelant consumption

            //printf( "T %g T/m %g acc %g m %g dm %g \n", thrust,thrust/m, acc.norm(), m, dm  );
            acc.add_mul( dir, thrust/m );
        }
        //exit(0);
    };

    // TODO : air drag

    ((Vec3d*) dYs   )->set( vel );
    ((Vec3d*)(dYs+3))->set( acc );
    dYs[6] = dm;
};


// ==== Class

class TestApp_SpaceFlightODE : public AppSDL2OGL_3D {
	public:
    bool RUNNING = true;
    ODEintegrator_RKF45 odeint;

    Vec3d  *pos=NULL,*vel=NULL,*acc=NULL;
    double *mass=NULL;
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

    odeint.reallocate( 7 );
    odeint.getDerivs = getAcceleration;
    odeint.dt_max   = 0.005;
    odeint.dt_min   = 0.0001;
    odeint.dt_adapt = 0.001;

    ((Vec3d*)(odeint.invMaxYerr  ))->set(1/1e-1);
    ((Vec3d*)(odeint.invMaxYerr+3))->set(1/1e-2);
    odeint.invMaxYerr[6] = 1/1e+1;

    pos  = (Vec3d*)(odeint.Y   );
    vel  = (Vec3d*)(odeint.Y+3 );
    acc  = (Vec3d*)(odeint.dY+3);
    mass = odeint.Y+6;

    pos->set(1.0,0.0,0.0); vel->set(0.0,0.0,0.0);
    //pos->set(0.0,100.0e+3,0.0); vel->set(7.9e+3,0.0,0.0); satelite
    *mass = mass_initial;

    for(int i=0; i<nCP; i++){
        double t = (i-1)*0.5;
        double   theta = M_PI_2/(1.0d + t*t*1.5 );
        if(i==0) theta = M_PI_2;
        thetaCPs [i] = theta;
        thrustCPs[i] = thrust_full;
        dirCPs   [i].set( cos(theta), sin(theta), 0.0 );
    }

    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //zoom = 0.5;
}

void TestApp_SpaceFlightODE::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	/*
	// --- fade out
	glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    //glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4f(0.5f,0.5f,0.5f, 0.05f);Draw2D::drawRectangle( {-20.0,-20.0},{20.0,20.0}, true);
	//glDisable(GL_BLEND);
	//glClear( GL_DEPTH_BUFFER_BIT );
    */

	// --- draw control spline
    glColor3f(1.0f,1.0f,1.0f); Draw3D::drawCurve( 0.0, 8.0, 100, splinePixel );
    //for(int i=0; i<nCP; i++){ Draw3D::drawPointCross( {i,thetaCPs[i]*5.0,0.0}, 0.1 ); }
    for(int i=0; i<nCP; i++){ Draw3D::drawPointCross( dirCPs[i], 0.1 ); }

    //RUNNING = 0;

    if(RUNNING){
        int nstep = 0;
        //odeint.step( 0.1d );
        //odeint.adaptive_step_RKF45( );
        //nstep = odeint.integrate_adaptive( odeint.dt_adapt, odeint.t+0.2d );
        nstep = odeint.integrate_adaptive( odeint.dt_adapt, odeint.t+5.0d );

        //odeint.step_RKF45( 0.1d );
        //odeint.step_euler( 0.1d );

        //printf( "%i %i (%g,%g,%g) %g %g\n ", frameCount, nstep, pos->x,pos->y,pos->z, odeint.t, odeint.dt_adapt );
        printf( " %g : (%g,%g,%g) %g %g \n ", odeint.t, pos->x,pos->y,pos->z, vel->norm(), *mass );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine    ( opos      , *pos ); Draw3D::drawPointCross( *pos, 0.1 );
        //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( (*vel)*2.0, *pos );

        if( (*pos-planet_pos).norm2()<(planet_R*planet_R) ){
            printf("crashed !!!\n");
            RUNNING = false;
        }
    }



    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine    ( opos*view_scale      , (*pos)*view_scale ); Draw3D::drawPointCross( *pos*view_scale, 0.1 );
    glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( (*vel)*2.0*view_scale, (*pos)*view_scale );

    opos = *pos;



    //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross( {0.0,0.0,0.0}, 0.5 );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross( {0.0,0.0,0.0}, 0.1 );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawPointCross( planet_pos*view_scale, 0.5 );

    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawCircleAxis( 64, planet_pos*view_scale, {0.0,planet_R*view_scale,0.0}, {0.0,0.0,1.0} );

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
















