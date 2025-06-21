#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "geom3D.h"
#include "Radiosity.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

// ============= Application

void drawObstacles( const TriangleRayTracer& rad, bool filled ){
    for( int i=0; i<rad.triangleObstacles.size(); i++ ){
        const Triangle3D& tri = rad.triangleObstacles[i];
        Draw3D::drawTriangle( tri.a, tri.b, tri.c, filled );
    }
}

void drawElements( const TriangleRayTracer& rad ){
    for( const SurfElement& el : rad.elements ){
        Draw3D::drawPointCross(el.pos, 0.01);
        Draw3D::drawVecInPos(el.normal*0.05, el.pos);
    }
}

void drawElements2( const TriangleRayTracer& rt, double elementSize ){
    glBegin(GL_LINES);
    for( const SurfElement& s: rt.elements ){
        Draw3D::vertex( s.pos );
        Draw3D::vertex( s.pos+s.normal );
    }
    glEnd();
}

void drawPanels( int n, const SurfElement* elems,  const double* vals, float sz, float sc ){
    glBegin(GL_QUADS);
    for(int i=0; i<n; i++){
        Vec3d p = elems[i].pos;
        Vec3d nr, up,fw;
        nr = elems[i].normal;
        elems[i].normal.getSomeOrtho(up,fw);
        up.mul(sz);
        fw.mul(sz);
        float c = 0.0f;
        if(vals) c = vals[i]*sc;
        glColor3f( c,c,c );
        glVertex3f( p.x+up.x, p.y+up.y, p.z+up.z );
        glVertex3f( p.x-fw.x, p.y-fw.y, p.z-fw.z );
        glVertex3f( p.x-up.x, p.y-up.y, p.z-up.z );
        glVertex3f( p.x+fw.x, p.y+fw.y, p.z+fw.z );
    }
    glEnd();
}

class TestAppRadiosity : public AppSDL2OGL_3D { public:
    Radiosity solver;
    int ogl_complings;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppRadiosity( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRadiosity::TestAppRadiosity( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    ogl_complings = glGenLists(1);
    glNewList( ogl_complings, GL_COMPILE );

    solver.addTriangle( (Triangle3D){Vec3d{-1.0,-1.0,5.0},   Vec3d{1.0,-1.0,5.0},    Vec3d{0.0,1.0,5.0}},  0.5, true );
    solver.addTriangle( (Triangle3D){Vec3d{-1.0,-1.0,10.0},  Vec3d{1.0,-1.0,10.0},   Vec3d{0.0,1.0,10.0}}, 0.5, true );

    printf( " nElements %i \n", solver.elements.size() );

    solver.makeCouplingMatrix();
    solver.prepare();
    for( int i=0; i<solver.elements.size(); i++ ){ solver.sources[i]=0.0; solver.vals[i]=0.0; }
    for( int i=0; i<solver.elements.size(); i++ ){ solver.sources[i]=1.0; }

    glEndList();
}

void TestAppRadiosity::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glDisable( GL_LIGHTING );
    glColor3f(0.0,0.0,0.0);
	drawElements2(solver, 0.02 );

    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    glColor3f(0.6,0.6,0.6);
    drawObstacles(solver, false);

    solver.step_Direct();
	drawPanels( solver.elements.size(), &solver.elements[0],  solver.vals, 0.2, 10.0 );
};


void TestAppRadiosity::eventHandling ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
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
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRadiosity( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}