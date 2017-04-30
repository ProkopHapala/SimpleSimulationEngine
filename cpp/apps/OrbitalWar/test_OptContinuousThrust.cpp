
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "Draw3D.h"
#include "SDL_utils.h"

#include "CubicBSpline.h"

#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "testUtils.h"

// ReImplementation of
// Java version
// /home/prokop/Dropbox/MyDevSW/Processing_arch/Kosmos/_Orbital/OrbitalOpt2D_cubicBsplineFIRE_ends

constexpr int nCP = 10;
double CPs[nCP] = { 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0};


class TestApp_OptContinuousThrust : public AppSDL2OGL_3D {
	public:

    //DEGUB
    //Vec3d *  dpos;
    //Vec3d * ddpos;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	//virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );
    //virtual void mouseHandling( );

	TestApp_OptContinuousThrust( int& id, int WIDTH_, int HEIGHT_ );

};

TestApp_OptContinuousThrust::TestApp_OptContinuousThrust( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

}

void TestApp_OptContinuousThrust::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


	glDisable ( GL_LIGHTING );

    int nsub = 10;
    double du =  1.0d/nsub;
    double x    = 0;
    double oy   = 0;
    double ody  = 0;
    double oddy = 0;

    double odynum  = 0;
    double oddynum = 0;
    for( int i=0; i<nCP-2; i++ ){
        double u = 0;
        double c0 = CPs[i  ];
        double c1 = CPs[i+1];
        double c2 = CPs[i+2];
        double c3 = CPs[i+3];
        // spline is in inverse order
        //double c3 = CPs[i  ];
        //double c2 = CPs[i+1];
        //double c1 = CPs[i+2];
        //double c0 = CPs[i+3];
        //printf( " >> CP[%i] \n", i );
        for(int j=0; j<nsub; j++){
            double   y = CubicBSpline::val  ( u, c0, c1, c2, c3 ); glColor3f( 0.0,0.0,0.0 ); Draw2D::drawLine_d( {x-du,  oy},{x,  y} );
            double  dy = CubicBSpline::dval ( u, c0, c1, c2, c3 ); glColor3f( 0.0,0.0,1.0 ); Draw2D::drawLine_d( {x-du, ody},{x, dy} );
            double ddy = CubicBSpline::ddval( u, c0, c1, c2, c3 ); glColor3f( 1.0,0.0,0.0 ); Draw2D::drawLine_d( {x-du,oddy},{x,ddy} );

            double  dynum = (y-  oy)/du; glColor3f( 0.0,0.7,0.7 ); Draw2D::drawLine_d( {x-du,   odynum},{x,   dynum} );
            double ddynum = (dy-ody)/du; glColor3f( 0.7,0.7,0.0 ); Draw2D::drawLine_d( {x-du,  oddynum},{x,  ddynum} );
            //printf( " %g  %g %g %g   %g \n", x, y, dy, ddy, du );
            oy=y; ody=dy; oddy=ddy;
            odynum=dynum; oddynum=ddynum;
            u += du;
            x += du;
        }
    }

    //printf( "camMat.a (%3.3f,%3.3f,%3.3f) \n", camMat.a.x, camMat.a.y, camMat.a.z );

};

void TestApp_OptContinuousThrust::drawHUD(){
    glDisable ( GL_LIGHTING );

    //exit(0);
}

//void TestApp_OptContinuousThrust::keyStateHandling( const Uint8 *keys ){ };
/*


void TestApp_OptContinuousThrust::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY ); mouseY=HEIGHT-mouseY;
    mouse_t   = (mouseX-20 )/timeScale + tstart;
    mouse_val = (mouseY-200)/valScale;
    sprintf( curCaption, "%f %f\0", mouse_t, mouse_val );
    int ipoint_ = binSearchFrom<double>(mouse_t,splines.n,splines.ts);
    if( (splines.ts[ipoint_+1]-mouse_t)<(mouse_t-splines.ts[ipoint_]) ) ipoint_++;
    if(ipoint_!=ipoint){
        ipoint=ipoint_;
        char buff[100];
        Vec3d r,v;
        r.set( splines.CPs[0][ipoint],splines.CPs[1][ipoint],splines.CPs[2][ipoint] );
        v.set( splines.getPointDeriv(ipoint,0), splines.getPointDeriv(ipoint,1), splines.getPointDeriv(ipoint,2) );
        sprintf(buff, "%i %f r(%3.3f,%3.3f,%3.3f) v(%3.3f,%3.3f,%3.3f)", ipoint, splines.ts[ipoint], r.x,r.y,r.z,   v.x,v.y,v.z );
        txtStatic.inputText = buff;
    }
    //printf("%i %i %3.3f %3.3f \n", mouseX,mouseY, mouse_t, mouse_val );
}

void TestApp_OptContinuousThrust::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_f:  warrior1->tryJump(); break;
                //case SDLK_h:  warrior1->tryJump(); break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
                case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;

                case SDLK_x:  iedit=0; break;
                case SDLK_y:  iedit=1; break;
                case SDLK_z:  iedit=2; break;
            }
            break;
    };
    //AppSDL2OGL::eventHandling( event );
}

*/



// ===================== MAIN

TestApp_OptContinuousThrust * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestApp_OptContinuousThrust( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















