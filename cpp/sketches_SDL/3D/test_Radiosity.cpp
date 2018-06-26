
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
#include "raytrace.h"
//#include "Body.h"

#include "geom3D.h"
#include "Radiosity.h"


#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application


Triangle3D trinagles[] = {
 (Triangle3D){ (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, (Vec3d){0.0,1.0,0.0}  },
 (Triangle3D){ (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, (Vec3d){0.0,0.0,1.0}  }
};


void drawObstacles( Radiosity& rad ){
    for( int i=0; i<rad.triangleObstacles.size(); i++ ){
        Triangle3D& tri = rad.triangleObstacles[i];
        Draw3D::drawTriangle( tri.a, tri.b, tri.c );
    }
}

void drawElements( Radiosity& rad ){
    for( SurfElement& el : rad.elements ){
        Draw3D::drawPointCross(el.pos, 0.1);
        Draw3D::drawVecInPos(el.normal, el.pos);
    }
}

class TestAppRadiosity : public AppSDL2OGL_3D { public:
    Radiosity rad;
    int ogl_complings;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppRadiosity( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRadiosity::TestAppRadiosity( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    ogl_complings = glGenLists(1);
    glNewList( ogl_complings, GL_COMPILE);
    //rad.processTriangles( 2, trinagles, 0.05 );

    rad.addTriangle( (Triangle3D){ (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, (Vec3d){0.0,1.0,0.0}  },   0.1 );
    rad.addTriangle( (Triangle3D){ (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, (Vec3d){0.0,0.0,1.0}  },   0.1 );
    rad.triangleObstacles.push_back( (Triangle3D){ (Vec3d){0.2,0.2,0.2}, (Vec3d){1.0,0.0,0.0}, (Vec3d){0.0,1.0,1.0}  } );

    rad.makeCouplingMatrix();

    glEndList();
}

void TestAppRadiosity::draw(){
    //rintf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glColor3f(0.8,0.8,0.8);  drawElements( rad );
	glColor3f(0.0,0.0,0.0);  drawObstacles(rad);

	glColor3f(0.0,0.0,1.0); glCallList(ogl_complings);

    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    double t;
    Vec3d hitpos, normal;

};


void TestAppRadiosity::eventHandling ( const SDL_Event& event  ){
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


void TestAppRadiosity::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppRadiosity * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRadiosity( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















