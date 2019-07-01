
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

#include "Voronoi.h"

using namespace VoronoiNamespace;

class TestAppVoronoi : public AppSDL2OGL { public:

    Voronoi voronoi;
    Edges    edges;
    Vertices vertices;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppVoronoi( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppVoronoi::TestAppVoronoi( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    voronoi.edges    = &edges;
    //voronoi.vertices = &vertices;

    /*
    // warrning :  FOR THIS DEGENERATE POINTS SEEMS TO FAIL
	vertices.push_back( new (Vec2d){-0.1, -0.1} );
	vertices.push_back( new (Vec2d){ 0.1,  0.1} );
	vertices.push_back( new (Vec2d){-0.5,  1  } );
	vertices.push_back( new (Vec2d){-0.5, -1  } );
	vertices.push_back( new (Vec2d){ 1  , -0.5} );
	vertices.push_back( new (Vec2d){-1  , -0.5} );
	vertices.push_back( new (Vec2d){ 0  ,  1  } );
	vertices.push_back( new (Vec2d){ 0  , -1  } );
	vertices.push_back( new (Vec2d){ 1  ,  0  } );
	vertices.push_back( new (Vec2d){-1  ,  0  } );
	vertices.push_back( new (Vec2d){ 0.5,  1  } );
	vertices.push_back( new (Vec2d){ 0.5, -1  } );
	vertices.push_back( new (Vec2d){ 1  ,  0.5} );
	vertices.push_back( new (Vec2d){-1  ,  0.5} );
	vertices.push_back( new (Vec2d){-1  ,  1  } );
	vertices.push_back( new (Vec2d){-1  , -1  } );
	vertices.push_back( new (Vec2d){ 1  ,  1  } );
	vertices.push_back( new (Vec2d){ 1  , -1  } );
    */


    Vec2d span = (Vec2d){3.0,5.0};
    for(int i=0; i<100; i++){
        //ver.push_back( new (Vec2d){ randf(-span.x,span.x)  , randf(-span.y,span.y)  } );
        vertices.push_back( new (Vec2d){ randf(-span.x,span.x)  , randf(-span.y,span.y)  } );
    }


    voronoi.GetEdges(&vertices, 20, 20);
    printf( "edges.size() %i \n", edges.size()  );
    int n=0;
    for ( VEdge* e : edges ){
        Vec2d& p1 = *e->start;
        Vec2d& p2 = *e->end;
        printf( "edge[%i] %g,%g %g,%g \n", n, p1.x, p1.y, p1.x, p1.y );
        n++;
    };

}

void TestAppVoronoi::draw(){

    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    glDisable(GL_LIGHTING);

    // draw points
    glColor3f(1.0,1.0,1.0);
	for ( Vec2d* p : vertices ){
        Draw2D::drawPointCross_d( *p, 0.1 );
    }

    // draw edges
    glColor3f(0.0,0.0,1.0);
	glBegin(GL_LINES);
	for ( VEdge* e : edges ){
       Vec2d& p1 = *e->start;
       Vec2d& p2 = *e->end;
       glVertex2d( p1.x, p1.y );
       glVertex2d( p2.x, p2.y );
    }
	glEnd();

    // draw edges
    glColor3f(1.0,0.0,0.0);
	glBegin(GL_LINES);
	for ( VEdge* e : edges ){
       Vec2d& p1 = *e->left;
       Vec2d& p2 = *e->right;
       glVertex2d( p1.x, p1.y );
       glVertex2d( p2.x, p2.y );
    }
	glEnd();


    glFlush();
};

void TestAppVoronoi::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
}

void TestAppVoronoi::eventHandling ( const SDL_Event& event  ){
    //printf( "TestAppVoronoi::eventHandling() \n" );
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

void TestAppVoronoi::drawHUD(){

}

// ===================== MAIN

TestAppVoronoi * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new TestAppVoronoi( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















