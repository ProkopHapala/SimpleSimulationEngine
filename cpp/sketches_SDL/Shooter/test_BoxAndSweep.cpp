
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

#include "Shooter.h"
#include "broadPhaseCollision.h" // Move to Sooter later
#include "kBoxes.h"

#include "arrayAlgs.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

const int ndirs = 1;
Vec3d dirs[ndirs] = { {1.0,0.0,0.0} };

void drawSweepPermut( int n, int* permut, sweep::Span* data , float y0, float dy, float WIDTH ){
    float xmin  = data[permut[0  ]].xmax;
    float xmax  = data[permut[n-1]].xmin;
    float sc    = WIDTH/(xmax-xmin);
    glBegin(GL_LINES);
    for( int i=0; i<n; i++ ){
        float x0 = ( data[permut[i]].xmin-xmin )*sc;
        float x1 = ( data[permut[i]].xmax-xmin )*sc;
        glVertex3f(x0,y0+i*dy,0);
        glVertex3f(x1,y0+i*dy,0);
    }
    glEnd();
}


void BoxSweepCollideSelf( KBoxes& kboxes ){


}













class TestAppSweepAndPrune : public AppSDL2OGL_3D { public:
    Shooter world;

    int ipicked;
    int glo_burst = 0;
    //ObjectType     objType1;
    //VehicleType    vehType1;

    //SAPbuff sapObj,sapVeh;

    std::vector<Box> bodies;
    KBoxes kboxes;
    int*         Kpermut = 0;
    sweep::Span* Kintervals = 0;
    sweep::Span* Bintervals = 0;

    int oglSphereHitShape=0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppSweepAndPrune( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSweepAndPrune::TestAppSweepAndPrune( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

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

    int nk = kboxes.branches.size();
    Kpermut     = new int        [ nk ];
    Kintervals  = new sweep::Span[ nk ];
    Bintervals  = new sweep::Span[ n ];
    for(int i=0; i<nk; i++){ Kpermut[i]=i; };

    for(int i=0; i<nk; i++){
        Box& span =  kboxes.branches[i].span;
        Kintervals[i] = (sweep::Span){(float)span.a.x,(float)span.b.x};
    };

    int niter = insertSort( nk, Kpermut, Kintervals ); printf( "insertSort N: %i niters: %i \n",  nk, niter );
    niter     = insertSort( nk, Kpermut, Kintervals ); printf( "insertSort N: %i niters: %i \n",  nk, niter );

    printf("DEBUG 1 \n");
    for(int i=0; i<n; i++){
        Box& span =  kboxes.bodies[i];
        Bintervals[i] = (sweep::Span){(float)span.a.y,(float)span.b.y};
    };

    printf("DEBUG 2 \n");
    for(int i=0; i<nk; i++){
        KBox& B =  kboxes.branches[i];
        //printf( " \n", B.i0, B.n );
        insertSort( B.n, kboxes.permut+B.i0, Bintervals );
    };

    printf("DEBUG 3 \n");



    printf( "SETUP DONE \n"  );
}

void TestAppSweepAndPrune::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


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

};

void TestAppSweepAndPrune::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );

    glColor3f(0.0,0.5,0.0); drawSweepPermut( kboxes.branches.size(),Kpermut, Kintervals, 200, 3, WIDTH );

    srand(121);
    for(int i=0; i<kboxes.branches.size(); i++){
        Draw::color_of_hash( rand() );
        KBox& B =  kboxes.branches[i];
        drawSweepPermut( B.n, kboxes.permut+B.i0, Bintervals, 300, 3, WIDTH );
    };


}

void TestAppSweepAndPrune::eventHandling ( const SDL_Event& event  ){
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

TestAppSweepAndPrune * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppSweepAndPrune( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















