
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
#include "quaternion.h"
#include "raytrace.h"
//#include "Body.h"

#include "arrayAlgs.h"

#include "grids3D.h"
#include "HashMapT64.h"
#include "HashMap3D.h"
//#include "GridMap3D.h"   //   objects2cells instead of open-indexing

#include "Shooter.h"
#include "broadPhaseCollision.h" // Move to Sooter later
//#include "sweep.h"
//#include "kBoxes.h"
//#include "arrayAlgs.h"
//#include "AABBTree3D.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

/*

HashMap vs GridMap
 - Depending on filling density there are several variations of Grid Storage structures
    - Open-Indexing HashMap
    -
    - Array Based Map
    - Array Based SkipMap

*/

void collideSelfBruteForce(int n, Box* bodies, std::vector<Vec2i>& colPairs ){
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if( bodies[i].overlap( bodies[j] ) ){
                colPairs.push_back( {i,j} );
            }
        }
    }
}

class TestAppGridHash : public AppSDL2OGL_3D { public:
    Shooter world;



	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppGridHash( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppGridHash::TestAppGridHash( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {


    long t1; double T;


    printf( "SETUP DONE \n"  );
}

void TestAppGridHash::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


};

void TestAppGridHash::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );

}

void TestAppGridHash::eventHandling ( const SDL_Event& event  ){
    /*
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
            }
            break;
    };
    */
    AppSDL2OGL::eventHandling( event );
}


// ===================== MAIN

TestAppGridHash * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppGridHash( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















