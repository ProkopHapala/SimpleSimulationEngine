
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "SDL_utils.h"

#include "fastmath.h"
#include "Vec2.h"
#include "Solids.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

#include "AABBTree3D.h"

// ======================  TestApp

class TestAppAABBTree : public AppSDL2OGL_3D { public:

    // ---- function declarations
    virtual void draw   ();
    virtual void drawHUD();
    virtual void eventHandling ( const SDL_Event& event  );

    TestAppAABBTree( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppAABBTree::TestAppAABBTree( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

}

void TestAppAABBTree::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

};

void TestAppAABBTree::drawHUD(){
}

void TestAppAABBTree::eventHandling ( const SDL_Event& event  ){
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppAABBTree * testApp;

int main(int argc, char *argv[]){
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    int junk;
    testApp = new TestAppAABBTree( junk , 800, 600 );
    testApp->loop( 1000000 );
    return 0;
}
















