
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

#include "Solids.h"
#include "SoftBody.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application

class TestAppElasticity : public AppSDL2OGL_3D { public:
    //Radiosity rad;
    //int ogl_complings;

    SoftBodyLinearized body;

    bool bRun = true;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppElasticity( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppElasticity::TestAppElasticity( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    body.init( Solids::Tetrahedron.nedge, Solids::Tetrahedron.nvert, Solids::Tetrahedron.verts, Solids::Tetrahedron.edges );
    printf("DEBUG 1 \n");
    body.prepare( 0 );
    printf("DEBUG 2 \n");

    body.Fextern[0].set(0.5,1.0,0.25);
    body.Fextern[3].set(-0.5,-1.0,-0.25);

    for(int i=0; i<body.npoints; i++){
        body.disps[i]=body.Fextern[i];
    };

    VecN::print_vector( body.n, body.x );
    VecN::print_vector( body.n, body.b );
    //print_vector( body.npoints, body.r  );

    /*
    VecN::print_vector( body.n, body.r  );
    VecN::print_vector( body.n, body.r2 );
    VecN::print_vector( body.n, body.p  );
    VecN::print_vector( body.n, body.Ap );
    */

}

void TestAppElasticity::draw(){
    //rintf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	//Draw3D::drawT

	if(bRun){
        double err2  = body.step_CG();
        printf( "istep %i err2 %f \n", body.istep, sqrt(err2)  );

        bRun = false;

        VecN::print_vector( body.n, body.x );
        VecN::print_vector( body.n, body.b );
        //print_vector( body.npoints, body.r  );

        VecN::print_vector( body.n, body.r  );
        VecN::print_vector( body.n, body.r2 );
        VecN::print_vector( body.n, body.p  );
        VecN::print_vector( body.n, body.Ap );
    }

	//printf("DEBUG 3 \n");
	glBegin(GL_LINES);
	float clrsc = 1.0;
    for(int il=0; il<body.nsticks; il++ ){
        Vec2i ij   = body.ijs [il];
        Vec3d hat  = body.dirs[il];
        Vec3d& pi  = body.poss[ij.a];
        Vec3d& pj  = body.poss[ij.b];
        double dfl = hat.dot( pi - pj );
        glColor3f (1-clrsc*dfl, 0.0, clrsc*dfl );
        glVertex3f(pi.x,pi.y,pi.z);
        glVertex3f(pj.x,pj.y,pj.z);
        //printf( "%i  (%i,%i)  (%f,%f,%f)  (%f,%f,%f) \n", il, ij.a, ij.b,  pi.x,pi.y,pi.z,  pj.x,pj.y,pj.z );
    }
    glEnd();

    //printf("DEBUG 4 \n");
    glBegin(GL_LINES);
    float sc = 1.0;
    for( int i=0; i<body.npoints; i++ ){
        const Vec3d& p  = body.poss[i];
        const Vec3d& p_ = p + body.disps[i]*sc;
        glColor3f (1-clrsc*sc, 0.0, clrsc*sc );
        glVertex3f(p.x,p.y,p.z);
        glVertex3f(p_.x,p_.y,p_.z);
    }
    glEnd();

    //printf("DEBUG 5 \n");

	/*
	glColor3f(0.8,0.8,0.8);  drawElements (rad);
	glColor3f(0.0,0.0,0.0);  drawObstacles(rad);

	glColor3f(0.0,0.0,1.0); glCallList(ogl_complings);

    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    double t;
    Vec3d hitpos, normal;
    */
};


void TestAppElasticity::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_s:  bRun = !bRun; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppElasticity::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppElasticity * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppElasticity( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















