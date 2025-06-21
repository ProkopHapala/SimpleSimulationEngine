/// @file @brief  This demo implements a particle or light scattering simulation using the direct iterative approach found in `Scatterer2.h`. Instead of solving a large matrix like in radiosity, this method models the transport of flux through a network of pre-defined channels connecting different surface elements. It's an alternative approach to global illumination that can be more efficient for certain scene types.
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
#include "Scatterer2.h"
// #include "Scatterer.h" // Not needed if only using Scatterer2

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
        glVertex3f( p.x+up.x, p.y+up.y, p.z+up.z ); // Typo: fz.z should be fw.z
        glVertex3f( p.x-fw.x, p.y-fw.y, p.z-fw.z );
        glVertex3f( p.x-up.x, p.y-up.y, p.z-up.z );
        glVertex3f( p.x+fw.x, p.y+fw.y, p.z+fw.z );
    }
    glEnd();
}

// New function to draw ScatterElem2s based on their flux
void drawScatterElems(const Scattering2& solver, float sz, float sc) {
    glBegin(GL_QUADS);
    for (const auto& scatterElem : solver.scattelems) {
        Vec3d p = scatterElem.pos;
        Vec3d nr = Vec3dZ; // Simplified normal for drawing
        Vec3d up, fw;
        nr.getSomeOrtho(up, fw);
        up.mul(sz);
        fw.mul(sz);
        float c = scatterElem.fluxIn * sc; // Color based on fluxIn
        glColor3f(c, c, c);
        glVertex3f(p.x + up.x, p.y + up.y, p.z + up.z);
        glVertex3f(p.x - fw.x, p.y - fw.y, p.z - fw.z);
        glVertex3f(p.x - up.x, p.y - up.y, p.z - up.z);
        glVertex3f(p.x + fw.x, p.y + fw.y, p.z + fw.z);
    }
    glEnd();
}

// New function to draw channels based on their flux
void drawChannels(const Scattering2& solver, float line_width, float flux_scale) {
    glLineWidth(line_width);
    glBegin(GL_LINES);
    for (const auto& channel : solver.channels) {
        // Get positions of connected elements
        const Vec3d& p0 = solver.scattelems[channel.ends[0].elemIdx].pos;
        const Vec3d& p1 = solver.scattelems[channel.ends[1].elemIdx].pos;
        
        // Calculate average flux for coloring
        float flux_val = (channel.ends[0].fluxOut + channel.ends[1].fluxOut) * 0.5f;
        float c = flux_val * flux_scale;
        glColor3f(c, c, c); // Grayscale based on flux
        glVertex3f(p0.x, p0.y, p0.z);
        glVertex3f(p1.x, p1.y, p1.z);
    }
    glEnd();
    glLineWidth(1.0f); // Reset line width to default
}

class TestAppScattering : public AppSDL2OGL_3D { public:
    Scattering2 solver;
    int ogl_complings;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppScattering( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppScattering::TestAppScattering( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    ogl_complings = glGenLists(1);
    glNewList( ogl_complings, GL_COMPILE );

    solver.addTriangle( (Triangle3D){Vec3d{-1.0,-1.0,5.0},   Vec3d{1.0,-1.0,5.0},    Vec3d{0.0,1.0,5.0}},  0.5, true );
    solver.addTriangle( (Triangle3D){Vec3d{-1.0,-1.0,10.0},  Vec3d{1.0,-1.0,10.0},   Vec3d{0.0,1.0,10.0}}, 0.5, true );

    printf( "solver.elements.size() %i \n", solver.elements.size() );
    solver.setupScatterElems(Vec3d{0.1, 0.1, 0.1}, Vec3d{1.0, 1.0, 1.0}, 0.1);
    printf( "solver.scattelems.size() %i \n", solver.scattelems.size() );

    // Initialize some flux sources
    if (!solver.scattelems.empty()) {
        solver.scattelems[0].fluxIn = 10.0; // Source on the first element
        if (solver.scattelems.size() > 1) {
            solver.scattelems[1].fluxIn = 5.0; // Another source
        }
    }

    solver.makeChannles();
    printf( "solver.channels.size() %i \n", solver.channels.size() );

    glEndList();
}

void TestAppScattering::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glDisable( GL_LIGHTING );
    glColor3f(0.0,0.0,0.0);
	// drawElements2(solver, 0.02 ); // This draws normals of SurfElements, keep it if desired for debugging geometry

    glEnable( GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    glColor3f(0.6,0.6,0.6);
    drawObstacles(solver, false);

    solver.step_Direct();
	drawScatterElems(solver, 0.2, 1.0); // Draw scattering elements with flux

    drawChannels(solver, 2.0, 0.1); // Draw channels with flux
};


void TestAppScattering::eventHandling ( const SDL_Event& event  ){
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


void TestAppScattering::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppScattering * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppScattering( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}