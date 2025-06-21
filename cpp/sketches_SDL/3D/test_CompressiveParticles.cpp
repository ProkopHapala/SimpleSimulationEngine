
/// @file @brief  This program demonstrates a particle-based simulation of a compressible fluid, likely using a method similar to Smoothed Particle Hydrodynamics (SPH). It visualizes a high-velocity impact of the fluid particles against a static wedge-shaped boundary, showcasing effects like shockwaves, compression, and fluid dynamics without a mesh. The simulation is self-running, with camera controls for observation.
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "raytrace.h"
//#include "Body.h"
//#include "geom3D.h"

//#include "Solids.h"
//#include "SoftBody.h"
//#include "Truss.h"

#include "CompressiveParticles.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application

class TestAppCompressive: public AppSDL2OGL_3D { public:
    //Radiosity rad;
    //int ogl_complings;

    CompressiveParticles solver;

    bool bRun = false;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppCompressive( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppCompressive::TestAppCompressive( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    solver.realloc( 40, 2 );
    //solver.realloc( 1, 2 );
    solver.walls[0].fromPointAndNormal( {0.5,  1.0, 0.0}, {0.0,-1.0,0.0} );
    solver.walls[1].fromPointAndNormal( {0.5, -1.0, 0.0}, {0.0, 1.0,0.0} );
    solver.clearVelocity();
    double r = 0.1;
    for(int i=0; i<solver.n; i++){
        solver.pos [i].set( r*(i%4),r*(i/4 - 5), 0 );
        //solver.vpos[i] = Vec3dZero;
        solver.vpos[i] = Vec3dX * -1.0e+4;
        solver.setDensity( i, 534, 0.006, 10000, r*0.4 ); // lithium
    }

    /*
    solver.realloc( 6, 2 );
    solver.walls[0].fromPointAndNormal( {0.5,  1.0, 0.0}, {0.0,-1.0,0.0} );
    solver.walls[1].fromPointAndNormal( {0.5, -1.0, 0.0}, {0.0, 1.0,0.0} );
    solver.clearVelocity();
    double r = 0.1;
    for(int i=0; i<solver.n; i++){
        solver.pos [i].set( r*(i%2),r*(i/2), 0 );
        //solver.vpos[i] = Vec3dX * -1.0e+4;
        solver.vpos[i] = Vec3dX * 0;
        solver.setDensity( i, 534, 0.006, 10000, r*0.4 ); // lithium
    }
    */


    /*
    solver.realloc( 2, 0 );
    solver.pos [0]=Vec3dZero;
    solver.vpos[0]=Vec3dZero;
    solver.pos [1]={ 0.05,0.1,0.0};
    solver.vpos[1]={-1.0e+4,0.0,0.0};
    for(int i=0; i<solver.n; i++){
        solver.setDensity( i, 534, 0.006, 10000, 0.1 ); // lithium
    }
    */
}

void TestAppCompressive::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);

	double dt = 1e-7;
	if(bRun){
        printf( "running frame[%i] \n", frameCount );
        solver.clearForce();
        solver.eval();
        solver.moveMD( dt );

        for(int i=0; i<solver.n; i++){ solver.pos[i].z=0; } // fix
	}

	// --- drawParticles
	double vsc = 1.0e-4;
	double fsc = 1.0e-8;
    for(int i=0; i<solver.n; i++){
        //printf( "[%i] pos(%g,%g,%g) r %g \n", i, solver.pos[i].x, solver.pos[i].y, solver.pos[i].z, solver.Rs[i] );
        glColor3f(0.0,0.0,0.0); Draw3D::drawSphereOctLines( 16, solver.Rs[i], solver.pos[i] );
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( solver.fpos[i]*fsc, solver.pos[i] );
        //glColor3f(0.0,0.0,1.0); Draw3D::drawVecInPos( solver.vpos[i]*vsc, solver.pos[i] );
	}

	// --- draw walls
	glColor3f(1.0,1.0,1.0);
	for(int i=0; i<solver.nwall; i++){
        const Plane3D& w = solver.walls[i];
        Mat3d rot;
        rot.b = w.normal;
        rot.c = Vec3dZ;
        rot.a.set_cross ( rot.c, rot.b );
        Draw3D::drawPanel( w.normal*w.C, rot, {10,10} );
        Draw3D::drawVecInPos( w.normal*0.2, w.normal*w.C );
	}

	Draw3D::drawAxis( 5.0);
};


void TestAppCompressive::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE:  bRun = !bRun; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppCompressive::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppCompressive* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppCompressive( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
