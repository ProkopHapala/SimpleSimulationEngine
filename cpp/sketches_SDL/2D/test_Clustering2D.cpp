
/// @file @brief This program demonstrates a 2D clustering algorithm, likely K-means or a similar iterative method, using `Clustering2D.h`. It generates random points and iteratively refines their cluster assignments, allowing the user to observe the clustering process step-by-step by pressing the `SPACE` key.


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "testUtils.h"

#include "clustering3D.h"

/**
 * @class TestAppClustering2D
 * @brief This class represents a test application for 2D clustering.
 */
class TestAppClustering2D: public AppSDL2OGL { public:

    double R = 3.0;
    ClusterMap    cmap;
    ClusterMapBox cbmap;
    std::vector<Vec3d> ps;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppClustering2D( int& id, int WIDTH_, int HEIGHT_ );

};

/**
 * @brief Constructor for the TestAppClustering2D class.
 * @param id The ID of the application.
 * @param WIDTH_ The width of the application window.
 * @param HEIGHT_ The height of the application window.
 */
TestAppClustering2D::TestAppClustering2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    //std::ve
    float d = 10.0;
    float dd = 1.0;
    Vec3d p;
    for(int i=0; i<200; i++){
        //ps.push_back( (Vec3d){ randf(-d,d), randf(-d,d), 0.0 } );
        p.add( randf(-dd,dd), randf(-dd,dd), 0.0 );
        ps.push_back( p );
        printf( "%i %f %f %f \n", ps.size(), ps.back().x, ps.back().y, ps.back().z );
    }
}

/**
 * @brief Draws the application.
 */
void TestAppClustering2D::draw(){
    long tstart = getCPUticks();
    glClearColor( 0.9f, 0.9f, 0.9f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glColor3f(0.0,0.0,0.0);

	if( frameCount<ps.size() ){
        int i = frameCount;
        cmap .insert_dist( i, ps[i], R );
        cbmap.insert     ( i, ps[i], R, R*0.1 );
        printf( "%i cmap.clusters.size() %i \n", i, cmap.clusters.size() );
	}

    int i=0;
    /*
    for( Cluster& c : cmap.clusters ){
        Draw::color_of_hash(15454*i+1545454);
        Draw2D::drawCircle_d( c.center.xy(),R,32,false);
        for( int ip : c.leafs ){
            Vec3d p =  ps[ ip ];
            Draw2D::drawPointCross_d( p.xy(), 0.3 );
        }
        i++;
    }
    */
    for( ClusterBox& c : cbmap.clusters ){
        Draw::color_of_hash(15454*i+1545454);
        Draw2D::drawRectangle_d( c.bbox.a.xy(), c.bbox.b.xy(), false );
        for( int ip : c.leafs ){
            Vec3d p =  ps[ ip ];
            Draw2D::drawPointCross_d( p.xy(), 0.3 );
        }
        i++;
    }

};

/**
 * @brief Handles mouse events.
 */
void TestAppClustering2D::mouseHandling( ){
    Uint32 buttons = SDL_GetMouseState( &mouseX, &mouseY );
    defaultMouseHandling( mouseX, mouseY );
}

/**
 * @brief Handles SDL events.
 * @param event The SDL event.
 */
void TestAppClustering2D::eventHandling ( const SDL_Event& event  ){
    //printf( "TestAppClustering2D::eventHandling() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

/**
 * @brief Draws the heads-up display (HUD).
 */
void TestAppClustering2D::drawHUD(){

}

// ===================== MAIN

TestAppClustering2D* thisApp;

/**
 * @brief The main function.
 * @param argc The number of command-line arguments.
 * @param argv The command-line arguments.
 * @return The exit status of the program.
 */
int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	thisApp = new TestAppClustering2D( junk , 800, 600 );
	thisApp->zoom = 30;
	thisApp->loop( 1000000 );
	return 0;
}
















