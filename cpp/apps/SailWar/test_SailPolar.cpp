
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "fastmath.h"
#include "Vec2.h"
//#include "geom2D.h"
#include "drawMath2D.h"
#include "AppSDL2OGL.h"

#include "Body.h"
#include "Body2D.h"
#include "AeroSurf2D.h"
#include "Yacht2D.h"

#include "testUtils.h"

// ======================  TestApp

class TestAppSailPolar : public AppSDL2OGL {
	public:

	Yacht2D * thisShip;

    bool   reseting   = true;
	double rudder_phi = 0.02  +     M_PI/2;
    double mast_phi   = 0.7;
    double ship_phi   = M_PI/2;
    double wind_speed = 10.0;

    constexpr static int npolar = 100;
    double phis    [npolar];
    double mastCL  [npolar];
    double mastCD  [npolar];
    double rudderCL[npolar];
    double rudderCD[npolar];
    double hullCL  [npolar];
    double hullCD  [npolar];

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling   ( const SDL_Event& event          );
	virtual void keyStateHandling( const Uint8 *keys               );
	void evalComponentPolars     ( double phi_min, double phi_max  );
	void printSailSetup( );
	TestAppSailPolar( int& id, int WIDTH_, int HEIGHT_ );
};

void TestAppSailPolar::evalComponentPolars( double phi_min, double phi_max ){
    double dphi=(phi_max-phi_min)/(npolar-1);
    for( int i=0; i<npolar; i++ ){
        phis[i] = phi_min + dphi * i;
    }
    thisShip -> keel.eval_polar  ( npolar, phis, hullCD, hullCL     );
    thisShip -> rudder.eval_polar( npolar, phis, rudderCD, rudderCL );
    thisShip -> mast.eval_polar  ( npolar, phis, mastCD, mastCL     );
};

TestAppSailPolar::TestAppSailPolar( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

	const int npts = 4;
	static double poss[npts*2] = { -1.0, 0.0,   0.0, -0.1,   0.0, +0.1,   +1.0, 0.0  };
	static double mass[npts  ] = {  10.0, 50.0, 50.0, 10.0  };

	thisShip =  new Yacht2D();
		                              //  pos.x  pos.y  angle         area   CD0      dCD   dCDS   dCL     dCLS    sStall   wStall
    //thisShip -> keel  .fromString( "  0.0    0.0    1.57079632679 0.4    0.04     1.5   0.9    3.00    2.00    0.20     0.40"   );
	//thisShip -> rudder.fromString( " -1.1    0.0    1.77079632679 0.03   0.008    1.5   0.9    6.28    2.70    0.16     0.08"   );
	//thisShip -> mast  .fromString( " +0.05   0.0    0.0           3.0    0.1     -1.0   0.8    3.50    2.20    0.20     0.40"   );

	                              //  pos.x  pos.y  angle         area   CD0      dCD   dCDS   dCL     dCLS    sStall   wStall
    thisShip -> keel  .fromString( "  0.0    0.0    1.57079632679 0.4    0.15     1.5   0.9    3.00    2.00    0.20     0.40"   );
	//thisShip -> rudder.fromString( " -1.1    0.0    1.77079632679 0.01   0.02     1.5   0.9    6.28    2.70    0.16     0.08"   );
	thisShip -> rudder.fromString( " -1.1    0.0    1.77079632679 0.01   0.05     0.5   0.9    5.00    2.70    0.25     0.40"   );
	thisShip -> mast  .fromString( " +0.05   0.0    0.0           3.0    0.2     -1.0   0.8    3.50    2.20    0.20     0.40"   );
	thisShip -> from_mass_points( 2, mass, (Vec2d*)poss );
	thisShip -> setDefaults( );

    evalComponentPolars( -M_PI/2, M_PI/2 );

	int shipShape = glGenLists(1);
	glNewList( shipShape , GL_COMPILE );
	glBegin   (GL_TRIANGLE_FAN);
		glNormal3f( 0.0f, 0.0f, 1.0f );
		glVertex3f( +1.5,  0.0, 0 );
 		glVertex3f( +0.5,  0.2, 0 );
		glVertex3f( -1.0,  0.2, 0 );
 		glVertex3f( -1.0, -0.2, 0 );
		glVertex3f( +0.5, -0.2, 0 );
		glVertex3f( +1.5,  0.0, 0 );
	glEnd();
	glEndList();
	thisShip->shape = shipShape;

    zoom = 50.0;
}

void TestAppSailPolar::draw(){

    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glPushMatrix();
    glTranslatef( 0.0f, -30.0f, 0.0f );



    const int nsteps = 200;
    Vec2d vels [ nsteps ];
    Vec2d rots [ nsteps ];
    Vec2d poss [ nsteps ];
    thisShip->pos.set( 0.0d );
    if( reseting ){
        thisShip->setAngle    ( ship_phi                       );
        //thisShip->vel.set( 0.0d );
        thisShip->vel.set_mul ( thisShip->rot, wind_speed*0.5d );
    }
    thisShip->mast  .setAngle( mast_phi   );
    thisShip->rudder.setAngle( rudder_phi );
    //printf( "mast %f rudder %f keel %f\n", thisShip->mast.phi, thisShip->rudder.phi, thisShip->keel.phi );
    //glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( {0.0,0.0}, 0.10 );

    //glColor3f( 0.1f, 0.1f, 0.1f );   Draw2D::drawBody2d_d( thisShip->rot, thisShip->pos, 1.0, 0.5 );
    glScalef ( 4.0f, 4.0f, 4.0f );
    glColor3f( 0.7f, 0.7f, 0.7f );   thisShip->draw_shape();
    glColor3f( 0.9f, 0.2f, 0.2f );   thisShip->rudder.draw( *thisShip );
    glColor3f( 0.2f, 0.2f, 0.9f );   thisShip->mast  .draw( *thisShip );
    glScalef ( 0.25f, 0.25f, 0.25f );

    glColor3f( 0.5f, 0.5f, 0.5f );

    int    nsub  = 10;
    double dt    = 0.2;
    thisShip->testSail( nsteps, nsub, dt/nsub, -wind_speed, poss, vels, rots );
    Draw2D::drawLines ( nsteps, poss );
    Draw2D::drawPoints( nsteps, poss, 0.05 );

/*
    Vec2d vel_conv, rot_conv;
    int nstepmax  = 1000;
    thisShip->pos.set( 0.0d );
    thisShip->vel.set( 0.0d );
    glColor3f( 0.2f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( {0.0,0.0}, 0.10 );
    int nstepconv = thisShip->convergeSail( nstepmax, 20, 0.001, -10.0, 1e-3, 1e-3, vel_conv, rot_conv );
    glColor3f( 0.9f, 0.2f, 0.2f );  Draw2D::drawPointCross_d( {0.0,0.0}, 0.10 );
    if ( nstepconv < nstepmax ){
        printf( " convergeSail : %i (%3.3f,%3.3f) (%3.3f,%3.3f) \n", nstepconv, vel_conv, rot_conv );
    }else{
        printf( " convergeSail not converged int %i ! \n", nstepmax );
    }
*/

    /*
    const int nsteps = 50;
    double phi_rudder [ nsteps ];
    double phi_mast   [ nsteps ];
    double wind_speed [ nsteps ];
    Vec2d  vels       [ nsteps ];
    Vec2d  rots       [ nsteps ];
    glColor3f( 0.2f, 0.2f, 0.2f );
    for( int i=0; i<nsteps; i++ ){
        //phi_rudder[ i ] = 1.57079632679 +  -0.3 + ( i*( 0.3-(-0.3) )/ nsteps );
        //phi_mast  [ i ] = M_PI * 0.25;
        //phi_mast  [ i ] = 0.0;

        phi_rudder[ i ] = 1.57079632679 + 0.09;
        phi_mast  [ i ] = M_PI * ( -0.5   +  i/(float)nsteps  );
        phi_mast  [ i ] = M_PI * (  0.0   +  i/(float)nsteps  );

        //phi_rudder[ i ] = 1.57079632679 + 0.1 * i/(float)nsteps ;
        //phi_mast  [ i ] = M_PI * +0.25;

        wind_speed[ i ] = -10.0d;
        vels[ i ].set( 0.0d, 0.0d );
        rots[ i ].set( 1.0d, 0.0d ); rots[ i ].normalize();
    }
    thisShip->evalPolar( nsteps, 0.01, 1e-300, 1e-300,  phi_rudder, phi_mast, wind_speed, vels, rots, true );
    glScalef( 10.0, 10.0, 10.0 );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::drawLines( nsteps, vels );
    glScalef( 1/10.0, 1/10.0, 1/10.0 );
    */
   //STOP = true;

    glPopMatrix();
};

void TestAppSailPolar::drawHUD(){
    glPushMatrix();

    /*
    glTranslatef( 100.0f, 100.0f, 0.0f );
    glScalef( 30.0f, 30.0f, 30.0f );
    glColor3f( 0.8f, 0.8f, 0.8f ); Draw2D::drawLine( {-3.0,   0.0 }, {+3.0,   0.0 } ); Draw2D::drawLine( {  0.0, -3.0 }, {  0.0, +3.0 } );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::plot( npolar, phis, hullCD );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::plot( npolar, phis, hullCL );
    glTranslatef( 0.0f, 6.0f, 0.0f );
    glColor3f( 0.8f, 0.8f, 0.8f ); Draw2D::drawLine( {-3.0,   0.0 }, {+3.0,   0.0 } ); Draw2D::drawLine( {  0.0, -3.0 }, {  0.0, +3.0 } );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::plot( npolar, phis, mastCD );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::plot( npolar, phis, mastCL );
    glTranslatef( 0.0f, 6.0f, 0.0f );
    glColor3f( 0.8f, 0.8f, 0.8f ); Draw2D::drawLine( {-3.0,   0.0 }, {+3.0,   0.0 } ); Draw2D::drawLine( {  0.0, -3.0 }, {  0.0, +3.0 } );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::plot( npolar, phis, rudderCD );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::plot( npolar, phis, rudderCL );
    */

    glTranslatef( 150.0f, 150.0f, 0.0f );
    glScalef( 50.0f, 50.0f, 50.0f );
    glColor3f( 0.8f, 0.8f, 0.8f ); Draw2D::drawLine( {-3.0,   0.0 }, {+3.0,   0.0 } ); Draw2D::drawLine( {  0.0, -3.0 }, {  0.0, +3.0 } );

    glColor3f( 1.0f, 0.0f, 0.0f ); Draw2D::plot( npolar, phis, hullCD );   Draw2D::drawLine( {-3.0,   +3.0  }, {-2.5,  +3.0 } );
    glColor3f( 0.0f, 0.0f, 1.0f ); Draw2D::plot( npolar, phis, hullCL );   Draw2D::drawLine( {-3.0,   +2.9 },  {-2.5,  +2.9 } );

    glColor3f( 1.0f, 0.3f, 0.6f ); Draw2D::plot( npolar, phis, rudderCD ); Draw2D::drawLine( {-3.0,   +2.7 },  {-2.5,  +2.7 } );
    glColor3f( 0.6f, 0.3f, 1.0f ); Draw2D::plot( npolar, phis, rudderCL ); Draw2D::drawLine( {-3.0,   +2.6 },  {-2.5,  +2.6 } );

    glColor3f( 1.0f, 1.0f, 0.0f ); Draw2D::plot( npolar, phis, mastCD );   Draw2D::drawLine( {-3.0,   +2.4 },  {-2.5,  +2.4 } );
    glColor3f( 0.0f, 1.0f, 1.0f ); Draw2D::plot( npolar, phis, mastCL );   Draw2D::drawLine( {-3.0,   +2.3 },  {-2.5,  +2.3 } );

    glPopMatrix();
}

void TestAppSailPolar::printSailSetup( ){
    printf( "ship %f rudder %f mast %f reseting %i \n", ship_phi, rudder_phi, mast_phi, reseting );
};

void TestAppSailPolar::keyStateHandling( const Uint8 *keys ){

    if( keys[ SDL_SCANCODE_LEFT  ] ){ camX0 -= camStep; }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ camX0 += camStep; }
	if( keys[ SDL_SCANCODE_UP    ] ){ camY0 += camStep; }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ camY0 -= camStep; }

    if( keys[ SDL_SCANCODE_KP_7 ] ){ ship_phi   += 0.01;  printSailSetup( ); }
	if( keys[ SDL_SCANCODE_KP_4 ] ){ ship_phi   -= 0.01;  printSailSetup( ); }
    if( keys[ SDL_SCANCODE_KP_8 ] ){ rudder_phi += 0.001; printSailSetup( ); }
	if( keys[ SDL_SCANCODE_KP_5 ] ){ rudder_phi -= 0.001; printSailSetup( ); }
	if( keys[ SDL_SCANCODE_KP_9 ] ){ mast_phi   += 0.005; printSailSetup( ); }
	if( keys[ SDL_SCANCODE_KP_6 ] ){ mast_phi   -= 0.005; printSailSetup( ); }
};

void TestAppSailPolar::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_KP_0:  reseting = !reseting; printSailSetup( ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
};

// ===================== MAIN

TestAppSailPolar * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSailPolar( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















