
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
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
#include "HashMap64.h"
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

[Mticks]: build 1.73258 col 5.37771 brute 205.428 | ncol 3438 ref 3438
[Mticks]: build 1.69094 col 5.41739 brute 205.648 | ncol 3446 ref 3446
[Mticks]: build 1.90722 col 5.61058 brute 208.294 | ncol 3444 ref 3444
[Mticks]: build 1.73773 col 5.41648 brute 210.531 | ncol 3446 ref 3446
[Mticks]: build 1.73743 col 5.41527 brute 207.719 | ncol 3460 ref 3460
[Mticks]: build 1.69063 col 5.52472 brute 213.158 | ncol 3404 ref 3404
[Mticks]: build 1.69083 col 5.43928 brute 207.623 | ncol 3368 ref 3368
[Mticks]: build 1.78166 col 5.47265 brute 207.541 | ncol 3394 ref 3394
[Mticks]: build 1.81032 col 5.772 brute 206.514 | ncol 3402 ref 3402
[Mticks]: build 1.83466 col 5.42177 brute 207.632 | ncol 3376 ref 3376
[Mticks]: build 2.22474 col 6.59574 brute 213.105 | ncol 3430 ref 3430



*/


void countBucket(const HashMap64& hmap, int h){
    int DEBUG_counter=0;
    int n = hmap.hits[h];
    int i = h;
    //printf( "countBucket hits[%i] == %i \n", h, n );
    while( n>0 ){
        if( hmap.hashs[i]==h ){ n--; }
        i=(i+1)&hmap.mask;
        DEBUG_counter++;
        if(DEBUG_counter>hmap.capacity){
            printf("hits[%i] == %i but just %i found \n", h, hmap.hits[h], h-n );
            return;
        }
    }
}

void countObjects(const HashMap64& hmap){
    int nI = 0;
    int nP = 0;
    int nH = 0;
    int nhitsum = 0;
    for(int i=0; i<hmap.capacity; i++){
        if( hmap.iboxs[i] != hmap.EMPTY_I ){ nI++; }
        if( hmap.store[i] != hmap.EMPTY_P ){ nP++; }
        if( hmap.hashs[i] != hmap.EMPTY_H ){ nH++; }
        nhitsum += hmap.hits[i];
    }
    printf( "countObjects: filled %i nI %i nP %i nH %i nhitsum %i \n", hmap.filled, nI, nP, nH, nhitsum );
    for(int i=0; i<hmap.capacity; i++){
        countBucket(hmap, i);
    }
}

void collideBruteForce(int n, int o, const Box& b, Box* boxes, std::vector<Vec2i>& colPairs ){
    for(int j=0; j<n; j++){
        if( j==o ) continue;
        if( b.overlap( boxes[j] ) ){
            colPairs.push_back( {o,j} );
            //printf( "Ref[%i] (%g,%g,%g) (%g,%g,%g) j %i\n", colPairs.size(), boxes[j].a.x, boxes[j].a.y, boxes[j].a.z,     boxes[j].b.x, boxes[j].b.y, boxes[j].b.z, j );
        }
    }
}

void collideSelfBruteForce(int n, Box* bodies, std::vector<Vec2i>& colPairs ){
    for(int i=0; i<n; i++){
        //collideSelfBruteForce(n, bodies[i], bodies, colPairs );
        for(int j=i+1; j<n; j++){
            if( bodies[i].overlap( bodies[j] ) ){
                colPairs.push_back( {i,j} );
            }
        }
    }
}

class TestAppGridHash : public AppSDL2OGL_3D { public:
    Shooter world;

    HashMap3D hmap;
    std::vector<Box> bodies;
    //std::vector<Box> bodies2;

    Box cursor;
    Vec2i * colPairs = 0;
    int ncol = 0;

    std::vector<Vec2i> colPairsRef;

    bool bMove = true;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

	TestAppGridHash( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppGridHash::TestAppGridHash( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    long t1; double T;

    srand(15454);

    Vec3i ip = {-54,-69,-15}, ip_;
    uint64_t il;

    il = hmap.ruler.ixyz2long(  ip  );
    hmap.ruler.long2ixyz( il, ip_ );
    printf( " (%i,%i,%i) %li (%i,%i,%i) \n", ip.x, ip.y, ip.z,  il,  ip_.x, ip_.y, ip_.z );

    Vec3d p = {-54.05,-69.69,-15.66}, p_;
    il = hmap.ruler.pos2long(  p  );
    hmap.ruler.long2pos(  il, p_  );
    printf( " (%g,%g,%g) %li (%g,%g,%g) \n", p.x, p.y, p.z,  il,  p_.x, p_.y, p_.z );
    //exit(0);

    //Box mapSpan; mapSpan.setSymetric(2.0);    mapSpan.scale( (Vec3d){1.0,1.0,0.1} );
    Box span; span.setSymetric(1.0); //span.a=Vec3dZero;
    //Box span; span.setSymetric(0.5); span.a=Vec3dZero;
    //Box span; span.setSymetric(1.0); span.a=Vec3dZero;
    printf( "span (%g,%g,%g)  (%g,%g,%g) \n", span.a.x,span.a.y,span.a.z,   span.b.x,span.b.y,span.b.z );

    int nbodies = 8000;
    //for(int i=0;i<160;i++){
    Vec3d c = Vec3dZero;
    for(int i=0;i<nbodies;i++){
        c.add_mul( span.genRandomSample(), 1.3 );
        bodies.push_back( {c, c+span.genRandomSample()} );
        bodies.back().order();
    }

    //bodies.push_back( { (Vec3d){-3.5,0.5,-6.5}, (Vec3d){-2.6,0.6,6.6}} );

    int ncap = hmap.realloc( bodies.size() );
    hmap.ruler.setStep( 2.5 );
    hmap.setBodyBuff( bodies.size() );
    t1 = getCPUticks();
    for(int i=0;i<bodies.size();i++){
        hmap.insert(i,bodies[i]);
    }
    T = getCPUticks()-t1; printf( "build hmap    %g [Mticks]  %g [/body] \n", T*1.0e-6, T/bodies.size() );

    printf( "hmap power %i capacity %i filled %i ratio %g \n", hmap.hmap.power, hmap.hmap.capacity, hmap.hmap.filled, hmap.hmap.filled/(float)hmap.hmap.capacity );

    countObjects(hmap.hmap);
    //printf( "hmap.DEBUG_ninserts %i \n", hmap.DEBUG_ninserts );

    cursor.a.set(-0.5,0.23,1.2);
    cursor.b.set(0.5,0.8,0.2);
    cursor.order();

    printf( "boides:    n %i nbytes %i \n", bodies.size(),      bodies.size()*sizeof(Box) );
    printf( "hmap.hmap: n %i nbytes %i \n", hmap.hmap.capacity, hmap.hmap.nbytes()        );

    colPairs = new Vec2i[10000];
    ncol     = hmap.collide( -1, cursor, bodies.data(), colPairs );

    collideBruteForce( bodies.size(), -1, cursor, bodies.data(), colPairsRef );

    printf( "ncol %i ref %i \n", ncol, colPairsRef.size() );

    // ==== colide NxN (self)

    ncol = 0;
    t1 = getCPUticks();
    for(int i=0; i<bodies.size(); i++){
        ncol += hmap.collide( i, bodies[i], bodies.data(), colPairs );
    }
    T = getCPUticks()-t1; printf( "self+NxN hmap    %g [Mticks] \n", T*1.0e-6 );

    colPairsRef.reserve( bodies.size()*bodies.size() );
    colPairsRef.clear();
    t1 = getCPUticks();
    for(int i=0; i<bodies.size(); i++){
        collideBruteForce( bodies.size(), i, bodies[i], bodies.data(), colPairsRef );
    }
    T = getCPUticks()-t1; printf( "self+NxN brute   %g [Mticks] \n", T*1.0e-6 );

    printf( "self-NxN: ncol %i ref %i \n", ncol, colPairsRef.size() );

    //exit(0);

    printf( "SETUP DONE \n"  );
}

void TestAppGridHash::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    if( bMove ){
        float dt = 0.1;
        float fvel  = 0.3;
        float fdiff = 0.1;
        for( Box& b : bodies){
            Vec3d vel,diff;
            vel.set( randf(-fvel,fvel), randf(-fvel,fvel), randf(-fvel,fvel) );
            diff.set( randf(-fdiff,fdiff), randf(-fdiff,fdiff), randf(-fdiff,fdiff) );
            b.a.add_lincomb( dt, vel, dt, diff );
            diff.set( randf(-fdiff,fdiff), randf(-fdiff,fdiff), randf(-fdiff,fdiff) );
            b.b.add_lincomb( dt, vel, dt, diff );
            b.order();
        }
    }

    {

        long   t1build, t1col, t1brute;
        double Tbuild,  Tcol,  Tbrute;

        // TODO : we can speed up update ~2x if we do not rebuild map each frame
        //        but instead we reisert bboxes which entered new cells
        //    - conservative update - boxes which leaved cell but did not entered new does not need to be updated

        t1build = getCPUticks();
        hmap.hmap.clear();
        for(int i=0;i<bodies.size();i++){
            hmap.insert(i,bodies[i]);
        }
        Tbuild = getCPUticks()-t1build;
        //printf( "hmap power %i capacity %i filled %i ratio %g \n", hmap.hmap.power, hmap.hmap.capacity, hmap.hmap.filled, hmap.hmap.filled/(float)hmap.hmap.capacity );

        ncol = 0;
        t1col = getCPUticks();
        for(int i=0; i<bodies.size(); i++){
            //printf( "col %i %i \n", i, ncol );
            ncol += hmap.collide( i, bodies[i], bodies.data(), colPairs );
        }
        Tcol = getCPUticks()-t1col;

        colPairsRef.clear();
        t1brute = getCPUticks();
        for(int i=0; i<bodies.size(); i++){
            collideBruteForce( bodies.size(), i, bodies[i], bodies.data(), colPairsRef );
        }
        Tbrute = getCPUticks()-t1brute;

        printf( "[Mticks]: build %g col %g brute %g | ncol %i ref %i \n", Tbuild*1.0e-6, Tcol*1.0e-6, Tbrute*1.0e-6, ncol, colPairsRef.size() );

        if( ncol != colPairsRef.size() ){
            printf( "ERROR ncol(%i) != ref(%i) \n", ncol, colPairsRef.size() );
            countObjects( hmap.hmap );

            for(int i=0; i<bodies.size(); i++){
                ncol = hmap.collide( i, bodies[i], bodies.data(), colPairs );
                colPairsRef.clear();
                collideBruteForce( bodies.size(), i, bodies[i], bodies.data(), colPairsRef );
                if( ncol != colPairsRef.size() ){
                    printf( "body[%i]   ncol(%i) != ref(%i) \n", i, ncol, colPairsRef.size() );
                }
            }
            exit(0);
            bMove = false;
        }
    }

	glColor3f(1.0,0.0,0.0);
	Draw3D::drawBBox(cursor.a,cursor.b);

    glColor3f(0.0,0.0,0.0);
    for( const Box& b: bodies ){
        Draw3D::drawBBox(b.a,b.b);
    }

    glColor3f(0.0,0.4,0.0);
    HashMap64& hmp = hmap.hmap;
    for( int i=0; i<hmp.capacity; i++ ){
        if( hmp.store[i] != hmp.EMPTY_P ){
            Vec3i ip;
            Vec3d p0,p1;
            //printf( "%i not empty \n", i );
            //hmap.ruler.long2pos( hmp.iboxs[i], pos );
            hmap.ruler.long2ixyz( hmp.iboxs[i], ip );
            hmap.ruler.ixyz2pos(ip,p0);
            hmap.ruler.ixyz2pos({ip.x+1,ip.y+1,ip.z+1},p1);
            //printf( "%i (%g,%g,%g) (%g,%g,%g)\n", i, p0.x, p0.y, p0.z,    p1.x, p1.y, p1.z   );
            Draw3D::drawBBox(p0,p1);
        }
    }

    ncol = hmap.collide( -1, cursor, bodies.data(), colPairs );

    {
        colPairsRef.clear();
        collideBruteForce( bodies.size(), -1, cursor, bodies.data(), colPairsRef );
        if( ncol != colPairsRef.size() ){
            printf( " ncol(%i) != ref(%i) \n", ncol, colPairsRef.size() );
        }
    }
    glColor3f(1.0,1.0,0.0);
    for( int i=0; i<ncol; i++ ){
        //printf( "col pair (%i,%i)\n",  p.x, p.y );
        Draw3D::drawLine( cursor.center(), bodies[colPairs[i].b].center() );
    }

    //exit(0);
};

void TestAppGridHash::drawHUD(){
    glDisable ( GL_LIGHTING );

    float y0    = 200;
    float hstep = 3;
    HashMap64& hmp = hmap.hmap;
    for( int i=0; i<hmp.capacity; i++ ){
        int hn = hmp.hits[i];
        glColor3f( 0.0f, 1.0f, 0.0f );
        Draw2D::drawLine( {i, y0}, {i, y0+hn*hstep} );

        int hh = hmp.hashs[i];
        //glColor3f( 0.0f, 1.0f, 0.0f );
        Draw::color_of_hash(hh);
        Draw2D::drawLine( {i, y0-hstep}, {i, y0-hstep*2} );

    }
    glColor3f( 1.0f, 1.0f, 1.0f );
    Draw2D::drawLine( {0, y0}, {hmp.capacity, y0} );

}

void TestAppGridHash::keyStateHandling( const Uint8 *keys ){
    float step = 0.1;
	if( keys[ SDL_SCANCODE_KP_6 ] ){ cursor.shift({+step,0.0,0.0}); }
	if( keys[ SDL_SCANCODE_KP_4 ] ){ cursor.shift({-step,0.0,0.0}); }
    if( keys[ SDL_SCANCODE_KP_8 ] ){ cursor.shift({0.0,+step,0.0}); }
	if( keys[ SDL_SCANCODE_KP_5 ] ){ cursor.shift({0.0,-step,0.0}); }
    if( keys[ SDL_SCANCODE_KP_7 ] ){ cursor.shift({0.0,0.0,+step}); }
	if( keys[ SDL_SCANCODE_KP_9 ] ){ cursor.shift({0.0,0.0,-step}); }
	AppSDL2OGL::keyStateHandling( keys );
};

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
















