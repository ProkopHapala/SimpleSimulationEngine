
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

#include "voronoi.h"


class TestAppVoronoi : public AppSDL2OGL { public:

    Voronoi* vdg;
    std::vector<Vec2d*> ver;
    std::vector<VEdge> edges;


	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppVoronoi( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppVoronoi::TestAppVoronoi( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

for (std::vector<Vec2d*>::iterator i = ver.begin(); i != ver.end(); i++)
		delete((*i));
	ver.clear();
	edges.clear();
    /*
	ver.push_back( new (Vec2d){-0.1, -0.1} );
	ver.push_back( new (Vec2d){ 0.1,  0.1} );
	ver.push_back( new (Vec2d){-0.5,  1  } );
	ver.push_back( new (Vec2d){-0.5, -1  } );
	ver.push_back( new (Vec2d){ 1  , -0.5} );
	ver.push_back( new (Vec2d){-1  , -0.5} );
	ver.push_back( new (Vec2d){ 0  ,  1  } );
	ver.push_back( new (Vec2d){ 0  , -1  } );
	ver.push_back( new (Vec2d){ 1  ,  0  } );
	ver.push_back( new (Vec2d){-1  ,  0  } );
	ver.push_back( new (Vec2d){ 0.5,  1  } );
	ver.push_back( new (Vec2d){ 0.5, -1  } );
	ver.push_back( new (Vec2d){ 1  ,  0.5} );
	ver.push_back( new (Vec2d){-1  ,  0.5} );
	ver.push_back( new (Vec2d){-1  ,  1  } );
	ver.push_back( new (Vec2d){-1  , -1  } );
	ver.push_back( new (Vec2d){ 1  ,  1  } );
	ver.push_back( new (Vec2d){ 1  , -1  } );
    */

    Vec2d span = (Vec2d){3.0,5.0};
    for(int i=0; i<100; i++){
        ver.push_back( new (Vec2d){ randf(-span.x,span.x)  , randf(-span.y,span.y)  } );
    }

	vdg = new Voronoi();
	double minY = -10;
	double maxY = 10;
	edges = vdg->ComputeVoronoiGraph(ver, minY, maxY);
	//delete vdg;
}

void TestAppVoronoi::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

    // draw points
    glColor3f(1,1,1);
	//glBegin(GL_POINTS);
	glColor3f(0.0f, 1.0f, 1.0f);
	for (std::vector<Vec2d*>::iterator i = ver.begin(); i != ver.end(); i++){
		//glVertex2d((*i)->x, (*i)->y);
        Draw2D::drawPointCross_d( **i, 0.1 );
    }
    //glEnd();

    // Draw Voronoi Edges
    glColor3f(0,0,1);
	glBegin(GL_LINES);
	glColor3f(0.0f, .8f, .5f);
	for (std::vector<VEdge>::iterator j = edges.begin(); j != edges.end(); j++){
		glVertex2d(j->VertexA.x, j->VertexA.y);
		glVertex2d(j->VertexB.x, j->VertexB.y);
	}
	glEnd();
    
    glColor3f(1,0,0);
    glBegin(GL_LINES);
    int ne = 0; 
    for( GraphEdge& e : vdg->graph_edges ){
        //printf( " %i (%g,%g)  (%g,%g) \n", e.x1, e.y1, e.x2, e.y2 );
        glVertex2d( e.x1, e.y1 );
        glVertex2d( e.x2, e.y2 );
        ne++;
    }
    //printf( "n edges %i \n", ne );
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
















