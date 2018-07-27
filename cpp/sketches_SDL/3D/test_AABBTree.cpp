
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
#include "kBoxes.h"


// ======================  TestApp

class TestAppAABBTree : public AppSDL2OGL_3D { public:

    std::vector<Box> bodies;
    KBoxes kboxes;

    // ---- function declarations
    virtual void draw   ();
    virtual void drawHUD();
    virtual void eventHandling ( const SDL_Event& event  );

    TestAppAABBTree( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppAABBTree::TestAppAABBTree( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    //double d    = 0.1
    //double xmax = 5.0;
    //Box mapSpan; mapSpan.a.set(-xmax,-xmax,-xmax); mapSpan.b.set(xmax,xmax,xmax);
    //Box span;    span.a.set(-d,-d,-d); span.b.set(d,d,d);

    Box mapSpan; mapSpan.setSymetric(5.0);    mapSpan.scale( (Vec3d){1.0,1.0,0.1} );
    Box span; span.setSymetric(0.5); span.a=Vec3dZero;
    printf( "span (%g,%g,%g)  (%g,%g,%g) \n", span.a.x,span.a.y,span.a.z,   span.b.x,span.b.y,span.b.z );

    for(int i=0;i<1600;i++){
        Vec3d c = mapSpan.genRandomSample();
        bodies.push_back( {c, c+span.genRandomSample()} );
    }

    int n = bodies.size();
    int K = (int)sqrt(n);
    kboxes.build( K, n, false, bodies.data() );

    kboxes.collideSelf();

    for(int i=0; i<kboxes.branches.size(); i++  ){
        const KBox& b = kboxes.branches[i];
        printf( "KBox %i : n %i V %g R %g \n", i, b.n, b.span.volume(), b.span.dimensions().norm2() );
    }

}

void TestAppAABBTree::draw   (){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    glColor3f(0.0,0.0,0.0);
    for( const Box& b: bodies ){
        Draw3D::drawBBox(b.a,b.b);
    }

    glColor3f(1.0,0.0,0.0);
    srand(121);
    for( const KBox& b: kboxes.branches ){
        Draw::color_of_hash( rand() );
        Draw3D::drawBBox(b.span.a,b.span.b);
    }

    for( const Vec2i& ab: kboxes.collisionPairs ){
        Vec3d ca = kboxes.bodies[ab.a].center();
        Vec3d cb = kboxes.bodies[ab.b].center();
        Draw3D::drawLine( ca, cb );
    }

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
















