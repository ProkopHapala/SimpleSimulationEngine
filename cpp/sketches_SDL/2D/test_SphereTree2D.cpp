
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"


#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"
#include "SDL_utils.h"

//#include "SphereTreeND.h"
#include "SphereTreeND_2.h"

// ======================  TestApp

constexpr int nLevel = 4;
double RIs[4]        = { 1.0, 0.5, 0.25, 0.125 };

class TestAppPlotting : public AppSDL2OGL{
	public:

    SphereTreeND store;
    int fontTex;

	virtual void draw   ();
    //virtual void eventHandling( const SDL_Event& event );

	TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ );
};

TestAppPlotting::TestAppPlotting( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    printf("sizeof(std::vector<int>) %i\n", sizeof(std::vector<int>) );

    store.init( nLevel, 2, 256, RIs );
    for(int i=0;i<store.nLevel; i++){
        printf( "%i %f %f \n", i, store.RIs[i], store.ROs[i] );
    }

    double * xs = new double[2];
    xs[0]=1.0; xs[1]=2.0;
    store.insert( 2, xs, 0, 0 );

    for(int ilevel=0;ilevel<store.nLevel; ilevel++){
        printf(" ---- level %i ", ilevel );
        for(int i=0;i<store.nodes[ilevel].size(); i++){
            SphereNodeND& nd = store.nodes[ilevel][i];
            printf( "%i (%3.3f,%3.3f) ", i, nd.center[i], nd.center[0] );
            for(int j=0;j<nd.branches.size(); j++){
                printf( "%i ", nd.branches[j] );
            }
            printf( "\n" );
        }
    }

    //store.init();

}

void TestAppPlotting::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

};

// ===================== MAIN

TestAppPlotting * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppPlotting( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}







