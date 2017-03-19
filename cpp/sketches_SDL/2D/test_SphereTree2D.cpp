
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


void drawSphereTreeND( SphereTreeND& store ){
    for(int ilevel=0;ilevel<store.nLevel; ilevel++){
        //printf(" ---- level %i \n", ilevel );
        double R = store.RIs[ilevel];
        for(int i=0;i<store.nodes[ilevel].size(); i++){
            SphereNodeND& nd = store.nodes[ilevel][i];
            Draw2D::drawCircle_d( *((Vec2d*)nd.center), R, 16, false );

            for(int j=0; j<nd.branches.size(); j++){
                int isub = nd.branches[j];
                if(ilevel<(store.nLevel-1)){
                    Draw2D::drawLine_d( *((Vec2d*) nd.center), *((Vec2d*) store.nodes[ilevel+1][isub].center) );
                }else{
                    //Draw2D::drawLine_d( *((Vec2d*) nd.center), *((Vec2d*) store.leafs[isub] ) );
                }
            }

        }
    }
}

void printSphereTreeND( SphereTreeND& store ){
    for(int ilevel=0;ilevel<store.nLevel; ilevel++){
        printf(" ---- level %i \n", ilevel );
        for(int i=0;i<store.nodes[ilevel].size(); i++){
            SphereNodeND& nd = store.nodes[ilevel][i];
            printf( "i %i c=(%3.3f,%3.3f) branches : ", i, nd.center[i], nd.center[0] );
            for(int j=0;j<nd.branches.size(); j++){
                printf( "%i ", nd.branches[j] );
            }
            printf( "\n" );
        }
    }
}

class TestAppSphereTree2D : public AppSDL2OGL{
	public:

    SphereTreeND store;
    int fontTex;

	virtual void draw   ();
    virtual void eventHandling( const SDL_Event& event );

    void testInsert( Vec2d p );

	TestAppSphereTree2D( int& id, int WIDTH_, int HEIGHT_ );
};

TestAppSphereTree2D::TestAppSphereTree2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    printf("sizeof(std::vector<int>) %i\n", sizeof(std::vector<int>) );

    store.init( nLevel, 2, 256, RIs );
    for(int i=0;i<store.nLevel; i++){
        printf( "%i %f %f   %f %f \n", i, store.RIs[i], store.ROs[i], store.R2Is[i], store.R2Os[i] );
    }

    double * xs;
    xs = new double[2]; xs[0]=1.0; xs[1]=1.0;
    store.insert( 2, xs, 0, 0 );

    /*
    xs = new double[2]; xs[0]=1.0; xs[1]=1.15;
    //store.insert( 2, xs, 1, 0 );
    store.insert( 2, xs );

    xs = new double[2]; xs[0]=1.0; xs[1]=1.0-0.35;
    store.insert( 2, xs );

    xs = new double[2]; xs[0]=1.75; xs[1]=1.00;
    store.insert( 2, xs );
    */


    printf("==== print store done ====\n");

    //store.init();

}

void TestAppSphereTree2D::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable( GL_DEPTH_TEST );

	glColor3f(0.5,0.5,0.5);
	drawSphereTreeND(  store );

};

void TestAppSphereTree2D::testInsert( Vec2d p ){
    printf( " ==== testInsert (%3.3f,%3.3f) \n", p.x, p.y );
    double * p_ = new double[2]; // to make sure p is permananetly sored on heap
    ((Vec2d*)p_)->set(p);
    store.insert( 2, p_ );
    printSphereTreeND( store );
}

void TestAppSphereTree2D::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    testInsert( {mouse_begin_x,mouse_begin_y} );
                break;
            }
            /*
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "left button pressed !!!! " );
                    world.picked = NULL;
                    break;
            }
            break;
            */
    };
    AppSDL2OGL::eventHandling( event );
};

// ===================== MAIN

TestAppSphereTree2D * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppSphereTree2D( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}







