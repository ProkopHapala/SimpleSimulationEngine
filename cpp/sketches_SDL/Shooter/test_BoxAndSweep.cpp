
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
//#include "broadPhaseCollision.h" // Move to Sooter later
#include "sweep.h"
#include "kBoxes.h"

#include "arrayAlgs.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

const int ndirs = 1;
Vec3d dirs[ndirs] = { {1.0,0.0,0.0} };

void drawSweepPermut( int n, int* permut, sweep::Span* data , float y0, float dy, float WIDTH ){
    float xmin  = +1e+100;
    float xmax  = -1e+100;
    for(int i=0; i<n; i++){ xmin=fmin(xmin,data[permut[i]].xmin);  xmax=fmax(xmax,data[permut[i]].xmax); }
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

void collideSelfBruteForce(int n, Box* bodies, std::vector<Vec2i>& colPairs ){
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if( bodies[i].overlap( bodies[j] ) ){
                colPairs.push_back( {i,j} );
            }
        }
    }
}

class TestAppSweepAndPrune : public AppSDL2OGL_3D { public:
    Shooter world;

    int ipicked;
    int glo_burst = 0;
    //ObjectType     objType1;
    //VehicleType    vehType1;

    //SAPbuff sapObj,sapVeh;

    int ncol;
    Int2* collisionPairs;

    std::vector<Box> bodies;
    KBoxes kboxes;
    //int*         Kpermut = 0;
    //sweep::Span* Kintervals = 0;
    //sweep::Span* Bintervals = 0;


    std::vector<Vec2i> colPairsRef;


    int oglSphereHitShape=0;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppSweepAndPrune( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSweepAndPrune::TestAppSweepAndPrune( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {


    long t1; double T;

    Box mapSpan; mapSpan.setSymetric(5.0);    mapSpan.scale( (Vec3d){1.0,1.0,0.1} );
    Box span; span.setSymetric(0.2); span.a=Vec3dZero;
    //Box span; span.setSymetric(0.5); span.a=Vec3dZero;
    //Box span; span.setSymetric(1.0); span.a=Vec3dZero;
    printf( "span (%g,%g,%g)  (%g,%g,%g) \n", span.a.x,span.a.y,span.a.z,   span.b.x,span.b.y,span.b.z );

    //for(int i=0;i<160;i++){
    for(int i=0;i<16000;i++){
        Vec3d c = mapSpan.genRandomSample();
        bodies.push_back( {c, c+span.genRandomSample()} );
    }


    t1 = getCPUticks();
    kboxes.build(  (int)sqrt(bodies.size()), bodies.size(), false, bodies.data(), true );  //printf(" DEBUG 1 \n");
    T = getCPUticks()-t1; printf( "kboxes.build   %g [Mticks] \n", T*1.0e-6 );

    //kboxes.printBranchInfo();

    //return;

    //kboxes.updateSweep(); printf(" DEBUG 2 \n");

    //kboxes.collisionPairs.resize( kboxes.nbodies * kboxes.nbodies );
    //kboxes.collisionPairs.clear();

    //collisionPairs = new Int2[ kboxes.nbodies * kboxes.nbodies ];

    //printf("DEBUG 4 \n");

    colPairsRef.clear();
    t1 = getCPUticks();
    collideSelfBruteForce( bodies.size(), bodies.data(), colPairsRef );
    T = getCPUticks()-t1; printf( "collideSelfBruteForce   %g [Mticks] | npairs %i \n", T*1.0e-6, colPairsRef.size() );

    kboxes.collisionPairs.clear();
    t1 = getCPUticks();
    ncol = kboxes.collideSelf();
    T = getCPUticks()-t1; printf( "kboxes.collideSelf      %g [Mticks] | npairs %i \n", T*1.0e-6, kboxes.collisionPairs.size() );

    kboxes.collisionPairs.clear();
    t1 = getCPUticks();
    ncol = kboxes.collideSelf();
    T = getCPUticks()-t1; printf( "kboxes.collideSelfNew   %g [Mticks] | npairs %i \n", T*1.0e-6, kboxes.collisionPairs.size() );

    kboxes.collisionPairs.clear();
    t1 = getCPUticks();
    ncol = kboxes.collideSelfSweep();
    T = getCPUticks()-t1; printf( "kboxes.collideSelfSweep %g | npairs %i \n", T*1.0e-6, kboxes.collisionPairs.size() );

    t1 = getCPUticks();
    kboxes.updateSweep();
    T = getCPUticks()-t1; printf( "kboxes.updateSweep %g  \n", T*1.0e-6 );

    kboxes.collisionPairs.clear();
    t1 = getCPUticks();
    ncol = kboxes.collideSelfSweep();
    T = getCPUticks()-t1; printf( "kboxes.collideSelfSweep %g | npairs %i \n", T*1.0e-6, kboxes.collisionPairs.size() );


    //kboxes.collideSelf();
    //printf( "found %i collision pairs\n", ncol );



    //printf( " collision pairs  %i brute %i \n", kboxes.collisionPairs.size(), colPairsRef.size()  );

    /*
    for( Vec2i& p: kboxes.collisionPairs ){ p.order(); printf( "colPair    (%i,%i)\n",  p.x, p.y ); }
    for( Vec2i& p: colPairsRef           ){ p.order(); printf( "colPairRef (%i,%i)\n",  p.x, p.y ); }

    int i1 = 12;
    int ib1;
    int ip1 = kboxes.findBody( 12, ib1 );
    printf( " %i -> (%i,%i) %i %i \n", i1, ib1, ip1, kboxes.branches[ib1].i0, kboxes.branches[ib1].i0+kboxes.branches[ib1].n  );

    int i2 = 58;
    int ib2;
    int ip2 = kboxes.findBody( 58, ib2 );
    printf( " %i -> (%i,%i) %i %i \n", i2, ib2, ip2, kboxes.branches[ib2].i0, kboxes.branches[ib2].i0+kboxes.branches[ib2].n  );
    */

    /*
    bDEBUG = true;
    KBox& Bi = kboxes.branches[0];
    KBox& Bj = kboxes.branches[3];
    //sweep::collideCross( int ni, int* permuti, Span* boundsi,  int nj, int* permutj, Span* boundsj, Int2* cols,  bool inv=false );
    int nc = sweep::collideCross(  Bi.n, kboxes.permut+Bi.i0, kboxes.Bintervals, Bj.n, kboxes.permut+Bj.i0, kboxes.Bintervals, collisionPairs );
    bDEBUG = false;
    */

    printf( "SETUP DONE \n"  );
}

void TestAppSweepAndPrune::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


	//if(false){
	if(true){
        float dt = 0.01;
        float fvel  = 0.3;
        float fdiff = 0.0;
        for( Box& b : bodies){
            Vec3d vel,diff;
            vel.set( randf(-fvel,fvel), randf(-fvel,fvel), randf(-fvel,fvel) );
            diff.set( randf(-fdiff,fdiff), randf(-fdiff,fdiff), randf(-fdiff,fdiff) );
            b.a.add_lincomb( dt, vel, dt, diff );
            diff.set( randf(-fdiff,fdiff), randf(-fdiff,fdiff), randf(-fdiff,fdiff) );
            b.b.add_lincomb( dt, vel, dt, diff );
        }

        //kboxes.update();

        long   t1up,t1col;
        double Tup,Tcol;

        bool bFullUpdate = (frameCount%5)==0;

        t1up = getCPUticks();
        //kboxes.updateStable();
        if(bFullUpdate){ kboxes.updateSweep();           }
        else           { kboxes.updateSweepStable(true); }
        Tup = getCPUticks()-t1up;

        t1col = getCPUticks();
        kboxes.collisionPairs.clear();
        kboxes.collideSelfSweep();
        Tcol = getCPUticks()-t1col;

        printf( "Tup %g  Tcol %g  ncol %i full %i \n", Tup*1.0e-6, Tcol*1.0e-6, kboxes.collisionPairs.size(), bFullUpdate );
    }
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

    glColor3f(0.0,1.0,0.0);
    //printf( "=========\n" );
    /*
    for( int i=0; i<ncol; i++){
        const Int2& p = collisionPairs[i];
        printf( "col pair (%i,%i)\n",  p.i, p.j );
        Draw3D::drawLine( kboxes.bodies[p.i].center(), kboxes.bodies[p.j].center() );
    }
    */
    for( const Vec2i& p: kboxes.collisionPairs ){
        //printf( "col pair (%i,%i)\n",  p.x, p.y );
        Draw3D::drawLine( kboxes.bodies[p.x].center(), kboxes.bodies[p.y].center() );
    }

    /*
    glColor3f(1.0,0.0,0.0);
    for( const Vec2i& p: colPairsRef ){
        //printf( "col pair (%i,%i)\n",  p.x, p.y );
        Draw3D::drawLine( kboxes.bodies[p.x].center(), kboxes.bodies[p.y].center() );
    }
    */

};

void TestAppSweepAndPrune::drawHUD(){
    glDisable ( GL_LIGHTING );
    glColor3f( 0.0f, 1.0f, 0.0f );

    glColor3f(0.0,0.5,0.0); drawSweepPermut( kboxes.branches.size(), kboxes.Kpermut, kboxes.Kintervals, 200, 3, WIDTH );

    srand(121);
    for(int i=0; i<kboxes.branches.size(); i++){
        Draw::color_of_hash( rand() );
        KBox& B =  kboxes.branches[i];
        drawSweepPermut( B.n, kboxes.permut+B.i0, kboxes.Bintervals, 300, 3, WIDTH );
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
















