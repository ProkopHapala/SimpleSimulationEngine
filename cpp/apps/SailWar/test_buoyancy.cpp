
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
#include "VecN.h"

#include "PolyLinear1d.h"

#include "geom2D.h"
#include "Convex2d.h"
#include "drawMath2D.h"
#include "AppSDL2OGL.h"

#include "buoyancy.h"

#include "testUtils.h"

// ======================  TestApp

class TestAppTerrainCubic : public AppSDL2OGL {
	public:

    double angle =0.0d;
    double dAngle=0.01d;

    constexpr static int nPhis = 100;
    double phis[ nPhis ];
    double Ms  [ nPhis ];

    constexpr static int np = 4;
    Convex2d * hull;
    double yLs [ np ];
    double yRs [ np ];
    double xs  [ np ];

	// ---- function declarations
	virtual void draw   ();
	virtual void drawHUD();
	void evalPolar( double phi_min, double phi_max, double displacement );
	TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ );
};


TestAppTerrainCubic::TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    hull = new Convex2d( np );
	hull->corners[0].set( -1.0, -1.0 );
    hull->corners[1].set( +1.0, -1.0 );
    hull->corners[2].set( +1.5, +1.0 );
	hull->corners[3].set( -1.5, +1.0 );

    evalPolar( -1.5, 1.5, 2.5 );
    //exit(0);
}

void TestAppTerrainCubic::evalPolar( double phi_min, double phi_max, double displacement ){
    double dphi = ( phi_max - phi_min ) / (nPhis-1);
    for( int i=0; i<nPhis; i++ ){
        double phi = phi_min + dphi * i;
        Vec2d dir; dir.set( sin( phi  ), cos( phi ) );
        hull->projectToLine( dir, xs, yLs, yRs );
        double watterline;
        double moment = integrate_moment( np, xs, yLs, yRs, displacement, watterline );
        phis[i] = phi;
        Ms  [i] = moment;
        printf( " %i %f %f \n", i, phi, moment );
    }
}

void TestAppTerrainCubic::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    // ============== project to line Convex Polygons

    glColor3f( 0.7f, 0.7f, 0.7f ); Draw2D::drawConvexPolygon( np, hull->corners, false );

    //Vec2d dir; dir.set( +1.0, 0.1 ); dir.normalize( );
    //ouble angle = ( frameCount%100 - 50 ) * 0.007;

    angle += dAngle;
    if       ( ( angle >  1.5 ) && ( dAngle > 0 ) ){ dAngle = -dAngle; }
    else if  ( ( angle < -1.5 ) && ( dAngle < 0 ) ){ dAngle = -dAngle; };

    Vec2d dir; dir.set( sin( angle  ), cos( angle ) );
    hull->projectToLine( dir, xs, yLs, yRs );
    //printArray( base.n, xs );
    //printArray( base.n, yLs );
    //printArray( base.n, yRs );
    double ys  [ np ];
    VecN::sub( np, yLs, yRs, ys );

/*
    PolyLinear1d pline( np, xs, ys );
    //pline.n  = np;
    //pline.xs = xs;
    //pline.ys = ys;
    //double V = pline.integrate( -1000.0, 1000.0 );
    double V     = pline.integrate_ibounds( 0, pline.n );
    double Vtarget = V*0.1;
    double xhalf   = pline.x_of_integral( Vtarget );
    double Vhalf   = pline.integrate( -1000.0, xhalf );
    printf( " V %f xhalf %f Vhalf %f  Vtarget %f \n", V, xhalf, Vhalf, Vtarget );
    pline.detach();
*/

    double displacement = 2.5;
    double watterline;
    double moment = integrate_moment( np, xs, yLs, yRs, displacement, watterline );
    //printf( " V %f xtarget %f %f  Vtarget %f %f moment %f \n", V, xhalf, watterline, Vhalf, Vtarget, moment );
    printf( " %i %f %f \n", frameCount, watterline, moment );

    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine( {  -10.0, watterline }, {  +10.0, watterline } );
    glColor3f( 0.9f, 0.9f, 0.9f ); Draw2D::drawLine( {  0.0, -10.0 }, {  0.0, +10.0 } ); Draw2D::drawLine( {-10.0,   0.0 }, {+10.0,   0.0 } );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::plot( np, yLs, xs );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::plot( np, yRs, xs );
    glColor3f( 0.0f, 0.0f, 0.0f ); Draw2D::plot( np, ys,  xs );

    glColor3f( 0.9f, 0.2f, 0.9f ); Draw2D::plot( nPhis, phis, Ms );

	// ============== Mouse Raycast  Convex Polygons

   //STOP = true;

};

void TestAppTerrainCubic::drawHUD(){}


// ===================== MAIN

TestAppTerrainCubic * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppTerrainCubic( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















