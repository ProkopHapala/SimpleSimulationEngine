
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

	// ---- function declarations

    std::vector<Convex2d*>   isles;

	virtual void draw   ();
	virtual void drawHUD();

	TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppTerrainCubic::TestAppTerrainCubic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    //exit(0);

    int ncorners = 5;
    int nisles   = 30;
    std::vector<double> phi_buf(ncorners);
    for( int i=0; i<nisles; i++ ){
        for( double& phi : phi_buf ){ phi = randf() * M_PI * 2; }
        std::sort( phi_buf.begin(), phi_buf.end() );
        Convex2d * isle = new Convex2d( ncorners );
        double x0 = randf( -10.0, 10.0 );
        double y0 = randf( -8.0, 8.0 );
        for( int i=0; i<ncorners; i++ ){
            double phi = phi_buf[i];
            isle->corners[i].set( x0 + cos( phi ), y0 + sin(phi) );
        }
        isle->update_lines();
        isles.push_back( isle );
    }
}

void TestAppTerrainCubic::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );


    // ============== project to line Convex Polygons

    const int np = 7;
    double yLs [np];
    double yRs [np];
    double ys  [np];
    double xs  [np];
    Convex2d base( np );
	//base.corners[0].set( -1.0, +1.0 );
    //base.corners[1].set( +1.0, -1.5 );
    //base.corners[2].set( +2.0,  0.5 );
	base.corners[0].set( -1.0, -1.0 );
    base.corners[1].set( +1.0, -1.0 );
    base.corners[2].set( +2.0,  0.0 );
	base.corners[3].set( +1.0, +1.0 );
	base.corners[4].set( -1.0, +1.0 );
	base.corners[5].set( -1.5, +0.5 );
	base.corners[6].set( -1.5, -0.5 );

    glColor3f( 0.5f, 0.5f, 0.5f ); Draw2D::drawConvexPolygon( base.n, base.corners, false );

    //Vec2d dir; dir.set( +1.0, 0.1 ); dir.normalize( );
    double angle = M_PI * ( ( frameCount % 600 ) / 300.0 ) + 0.5;
    Vec2d dir; dir.set( cos( angle  ), sin( angle ) );
    base.projectToLine( dir, xs, yLs, yRs );
    //printArray( base.n, xs );
    //printArray( base.n, yLs );
    //printArray( base.n, yRs );
    VecN::sub( np, yLs, yRs, ys );

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

    double watterline;
    double moment = integrate_moment( np, xs, yLs, yRs, Vtarget, watterline );
    printf( " V %f xtarget %f %f  Vtarget %f %f moment %f \n", V, xhalf, watterline, Vhalf, Vtarget, moment );

    pline.detach();
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine( {  xhalf, -10.0 }, {  xhalf, +10.0 } );

    glColor3f( 0.5f, 0.5f, 0.5f );
    Draw2D::drawLine( {  0.0, -10.0 }, {  0.0, +10.0 } );
    Draw2D::drawLine( {-10.0,   0.0 }, {+10.0,   0.0 } );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::plot( base.n, xs, yLs );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::plot( base.n, xs, yRs );
    glColor3f( 0.0f, 0.0f, 0.0f ); Draw2D::plot( base.n, xs, ys );



	// ============== Mouse Raycast  Convex Polygons

/*
    //printf( " mouse x,y : %i %i %f %f \n", mouseX, mouseY, mouseRight(mouseX), mouseUp(mouseY) );
    Vec2d mvec;
    mvec.set( mouseRight(mouseX), mouseUp(mouseY) );
    //mvec.set( 0.0, 0.0 );
    Draw2D::drawPointCross_d( mvec, 1.0 );

	for( auto isle : isles ) {
        //printf( "draw isle: %i %f %f \n ", isle->n, isle->corners[0].x, isle->corners[0].y );
        if( isle->pointIn( mvec ) ){
            glColor3f( 1.0f, 0.2f, 0.2f );
            //printf( " mouse in isle \n " );
        }else{
            glColor3f( 0.2f, 0.2f, 0.2f );
        }
		Draw2D::drawConvexPolygon( isle->n, isle->corners, true );

		//exit(0);
	}
*/
    // ============== Cut Convex Polygon

/*
    Convex2d base( 5 );

	//base.corners[0].set( -1.0, +1.0 );
    //base.corners[1].set( +1.0, -1.5 );
    //base.corners[2].set( +2.0,  0.5 );

	base.corners[0].set( -1.0, -1.0 );
    base.corners[1].set( +1.0, -1.0 );
    base.corners[2].set( +2.0,  0.0 );
	base.corners[3].set( +1.0, +1.0 );
	base.corners[4].set( -1.0, +1.0 );

    glColor3f( 0.9f, 0.9f, 0.9f ); Draw2D::drawConvexPolygon( base.n, base.corners, false );


    Vec2d Acut,Bcut;
    Acut.set( -1.0, -5.0 );
    Bcut.set( +3.0, +5.0 );
    Line2d cutline; cutline.set( Acut, Bcut );
    glColor3f( 0.2f, 0.9f, 0.2f ); Draw2D::drawLine_d( Acut, Bcut );

    Convex2d left;
    Convex2d right;

    base.cut( cutline, left, right );
    glColor3f( 0.9f, 0.2f, 0.9f ); Draw2D::drawPointCross_d( left.corners[0], 0.3 );
    glColor3f( 0.2f, 0.9f, 0.9f ); Draw2D::drawPointCross_d( left.corners[1], 0.3 );
    glColor3f( 0.9f, 0.2f, 0.2f ); Draw2D::drawConvexPolygon( left.n,  left.corners,  true );
    glColor3f( 0.2f, 0.2f, 0.9f ); Draw2D::drawConvexPolygon( right.n, right.corners, true );

*/

   STOP = true;

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
















