
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "Multipoles.h"
#include "PotentialFlow.h"
#include "grids3D.h"
#include "MultipoleGrid.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#define R2SAFE  1.0e-8f

// ======= THE CLASS

typedef Vec3d (*VecFieldFunc)(Vec3d R);

Vec3d vInf      = {0.0,0.0,-1.0};
Vec3d flightDir = vInf*(-1.0/vInf.norm());

const int nCPs  = 3;
Vec3d CPs[nCPs] = { {0.0,0.0,0.0}, {1.0,0.0,0.0}, {2.0,0.0,-0.5} };
double Gamas[nCPs-1] = { 1.0, 1.0 };

Vec3d totVelField( Vec3d R ){
    Vec3d B = (Vec3d){0.0,0.0,0.0};
    //Vec3d B = vInf;
    B.add( vInf );
    for(int i=0; i<nCPs-1; i++){
        //B.add( horseshoe( R, CPs[i], CPs[i+1], flightDir, Gamas[nCPs] ) );
        horseshoe( B, R, CPs[i], CPs[i+1], flightDir, Gamas[i] );
        //horseshoeDecay( B, R, CPs[i], CPs[i+1], flightDir, Gamas[i], 16.0 );
    }
    //exit(0);
    //printf("=====\n");
    return B;

    //return horseshoe( R, {1.0,0.0,0.0}, {-1.0,0.0,0.0}, {0.0,0.0,1.0}, 1.0 );
};

void plotVortexFilaments( int n, Vec3d* CPs, Vec3d dDir ){
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,1.0,1.0);
    for(int i=0; i<nCPs;i++){
        glVertex3f( CPs[i].x, CPs[i].y, CPs[i].z );
    }
    glEnd();
    glColor3f(0.75,0.75,0.75);
    glBegin(GL_LINES);
    for(int i=0; i<nCPs;i++){
        Vec3d p = CPs[i];
        glVertex3f( p.x, p.y, p.z );
        p.add( dDir );
        glVertex3f( p.x, p.y, p.z );
    }
    glEnd();

};

void plotVecPlane( Vec2i n, Vec3d p0, Vec3d a, Vec3d b, double sz, double dt, VecFieldFunc func ){
    glBegin(GL_LINES);
    for(int ia=0; ia<n.a; ia++ ){
        for(int ib=0; ib<n.b; ib++ ){
            Vec3d p = p0 + a*ia + b*ib;
            if( sz>0 ){
                glVertex3f( p.x-sz, p.y   , p.z    ); glVertex3f( p.x+sz, p.y,    p.z    );
                glVertex3f( p.x   , p.y-sz, p.z    ); glVertex3f( p.x,    p.y+sz, p.z    );
                glVertex3f( p.x   , p.y   , p.z-sz ); glVertex3f( p.x,    p.y,    p.z+sz );
            }
            glVertex3f( p.x, p.y, p.z );
            Vec3d v = func(p); //printf( "(%f,%f,%f) (%f,%f,%f)\n", p.x, p.y, p.z,  v.x, v.y, v.z );
            p.add_mul( v, dt);
            glVertex3f( p.x, p.y, p.z );
        }
    }
    glEnd();
}

void plotStreamLine( int n, double dt, Vec3d p, VecFieldFunc func ){
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++ ){
        Vec3d v = func(p);
        p.add_mul(v,dt);
        glVertex3f( p.x, p.y, p.z );
        //printf( "(%f,%f,%f) (%f,%f,%f)\n", p.x, p.y, p.z,  v.x, v.y, v.z );
    }
    glEnd();
    //exit(0);
}

void plotStreamLinePlane( Vec2i n, int ns, Vec3d p0, Vec3d a, Vec3d b, double dt, VecFieldFunc func ){
    for(int ia=0; ia<n.a; ia++ ){
        for(int ib=0; ib<n.b; ib++ ){
            Vec3d p = p0 + a*ia + b*ib;
            //printf( "(%i,%i) (%f,%f,%f) \n", ia, ib, p.x, p.y, p.z );
            plotStreamLine( ns, dt, p, func );
        }
    }
}

void testILineFinite(){
    // test ILineFinite
    //Vec3d L  = (Vec3d){4.0,0.5,0.6};
    //Vec3d R0 = (Vec3d){-3.0,1.0,2.0};

    Vec3d L  = (Vec3d){0.0,4.0,0.0};
    Vec3d R0 = (Vec3d){1.0,1.0,0.0};

    double l = L.norm();
    Vec3d hL = L*(1.0/l);
    //Vec3d B_ana = ILineFinite( R0, hL, l );
    Vec3d B_ana = ISemiInfSheet( R0, {1.0,0.0,0.0}, hL, l );

    int nsamp=1000;
    Vec3d B_num = (Vec3d){0.0,0.0,0.0};
    double step = 1.0/nsamp;
    Vec3d dL    = L*step;
    Vec3d R     = R0 + dL*0.5;
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,1.0);
    for(int i=0; i<nsamp; i++){
        //B_num.add( dBiotSawart( R, dL ) );
        B_num.add_mul( ILineSemiInf( R, {1.0,0.0,0.0} ), step*l );
        //printf( "%f %f %f\n", R.y, ILineSemiInf( R, {1.0,0.0,0.0} ).norm()/VortexQuantum, B_num.norm()/VortexQuantum );
        R.add(dL);
        glVertex3f( l*step*i, B_num.z*20.0, 0 );
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,0.0,0.0);
    for(int i=0; i<nsamp; i++){
        //B_ana = ILineFinite( R0, hL, l*step*i );
        B_ana = ISemiInfSheet( R0, {1.0,0.0,0.0}, hL, l*step*i );
        glVertex3f( l*step*i, B_ana.z*20.0, 0 );
    }
    glEnd();
    //exit(0);
    glColor3f(0.0,0.0,0.0);
    Draw3D::drawLine( {0.0,0.0,0.0}, {10.0,0.0,0.0} );

    //printf(" B_ana (%f,%f,%f) B_num (%f,%f,%f) \n", B_ana.x,B_ana.y,B_ana.z, B_num.x,B_num.y,B_num.z );
    printf(" |B_ana| %f |B_num| %f \n", B_ana.norm()/VortexQuantum, B_num.norm()/VortexQuantum );
    //exit(0);
}

class TestAppMultipoles : public AppSDL2OGL_3D {
	public:

    MultipoleGrid grid;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMultipoles( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMultipoles::TestAppMultipoles( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //exit(0);


}

void TestAppMultipoles::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//LineSemiInf( {0.5,1.5,0.0}, {1.0,0.0,0.0} ); exit(0);
	//ISemiInfSheet( {1.0,1.0,0.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, 4.0 ); exit(0);

	testILineFinite(); //exit(0);
    //if(frameCount==0)testILineFinite();
    //glColor3f(1.0,1.0,1.0);
    //Draw3D::drawLine    ( (Vec3d){1.0,0.0, 0.0},  {-1.0,0.0,0.0} );
    //Draw3D::drawVecInPos( (Vec3d){0.0,0.0,10.0},  { 1.0,0.0,0.0} );
    //Draw3D::drawVecInPos( (Vec3d){0.0,0.0,10.0},  {-1.0,0.0,0.0} );

    /*
    plotVortexFilaments( nCPs, CPs, flightDir*-100.0 );
    glColor3f(0.0,0.0,0.0);
    //plotVecPlane( {21,21}, { -1.0,0.0,-2.0 }, {0.2,0.0,0.0}, {0.0,0.0,0.2}, 0.02, 1.0, totVelField );
    //plotVecPlane( {5,5}, { -2.0,0.0,-2.0 }, {1.0,0.0,0.0}, {0.0,0.0,1.0}, 0.02, 10.0, totVelField );
    //exit(0);
    plotStreamLinePlane( {20,1}, 200, {-1.1,0.0,3.0}, {0.2,0.0,0.0}, {0.0,0.2,0.0}, 0.1, totVelField );
    */

};


void TestAppMultipoles::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppMultipoles::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppMultipoles * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMultipoles( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















