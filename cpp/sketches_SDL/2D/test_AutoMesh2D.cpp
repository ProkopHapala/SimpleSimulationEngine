
/**
 * @file test_AutoMesh2D.cpp
 * @brief This is example of 2D finite element simulation, and comparison of different methods implemnented in MechGrid2D.h, MechPIC2D.h, and MechMesh2D.h
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "MechGrid2D.h"
#include "MechPIC2D.h"
#include "MechMesh2D.h"

CompressibleMaterial materials[] = {
    {1.,1.}
};

class TestAppMech2D : public AppSDL2OGL { public:

    MechMesh2D mmesh;


    Mesh2D::Builder builder;

    //bool bRun = false;
    bool bRun = true;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppMech2D::TestAppMech2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    //                0            1           2        3
    double vs[]{ -1.0,-1.0,   1.0,-1.0,   -1.0,1.0,  1.0,1.0  };
    int    es[]{  0,1,   0,2,  1,3,  2,3,   0,3    };

    builder.fromArrayVerts( 4, (Vec2d*)vs);
    builder.fromArrayEdges( 5, (Vec2i*)es);

    //builder.sortEdges();
    for( Mesh2D::Edge e: builder.edges ){ printf( " %i %i \n", e.verts.a, e.verts.b ); };
    builder.makeTris();

}

void TestAppMech2D::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glDisable(GL_LIGHTING);

    glColor3f(1.0,0.0,0.0);

    double dt = 1e-7;
    if(bRun){


    }

    glBegin( GL_LINES );
    for( Mesh2D::Edge e : builder.edges ){
        Vec2d va = builder.verts[e.verts.a];
        Vec2d vb = builder.verts[e.verts.b];
        glVertex3f(va.x,va.y,0);
        glVertex3f(vb.x,vb.y,0);
    };
    glEnd();


    //glColor3f(0.,0.,0.);
    //double vsc   = 1e-3;
    //double molsc = 5.1;
    //for(int i=0; i<mpic.np; i++){
    //    Draw2D::drawPointCross_d( mpic.pos[i]*mpic.invStep, mpic.pmoles[i]*molsc );
    //    Draw2D::drawVecInPos_d  ( mpic.vel[i]*vsc, mpic.pos[i]*mpic.invStep );
    //}

};

void TestAppMech2D::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
}

void TestAppMech2D::eventHandling ( const SDL_Event& event  ){
    //printf( "TestAppMech2D::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //pickParticle( world.picked );
                    //ipick = pline1.nearestPoint( {mouse_begin_x,mouse_begin_y} );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    //world.picked = NULL;
                    //ipick = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppMech2D::drawHUD(){

}

// ===================== MAIN

TestAppMech2D * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new TestAppMech2D( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















