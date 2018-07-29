
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

void collideBruteForce(int n, int o, const Box& b, Box* bodies, std::vector<Vec2i>& colPairs ){
    for(int j=0; j<n; j++){
        if( j==o ) continue;
        if( b.overlap( bodies[j] ) ){
            colPairs.push_back( {o,j} );
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

    Box cursor;
    Vec2i * colPairs = 0;
    int ncol = 0;

    std::vector<Vec2i> colPairsRef;

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

    int nbodies = 50;
    //for(int i=0;i<160;i++){
    Vec3d c = Vec3dZero;
    for(int i=0;i<nbodies;i++){
        c.add_mul( span.genRandomSample(), 0.3 );
        bodies.push_back( {c, c+span.genRandomSample()} );
        bodies.back().order();
    }

    //bodies.push_back( { (Vec3d){-3.5,0.5,-6.5}, (Vec3d){-2.6,0.6,6.6}} );

    int ncap = hmap.realloc( bodies.size() );
    hmap.ruler.setStep( 1.0 );
    hmap.setBodyBuff( bodies.size() );
    for(int i=0;i<bodies.size();i++){
        hmap.insert(i,bodies[i]);
    }

    cursor.a.set(-0.5,0.23,1.2);
    cursor.b.set(0.5,0.8,0.2);
    cursor.order();

    printf( "boides:    n %i nbytes %i \n", bodies.size(),      bodies.size()*sizeof(Box) );
    printf( "hmap.hmap: n %i nbytes %i \n", hmap.hmap.capacity, hmap.hmap.nbytes()        );

    colPairs = new Vec2i[1000];
    ncol     = hmap.collide( -1, cursor, bodies.data(), colPairs );

    collideBruteForce( bodies.size(), -1, cursor, bodies.data(), colPairsRef );

    printf( "ncol %i ref %i \n", ncol, colPairsRef.size() );

    printf( "SETUP DONE \n"  );
}

void TestAppGridHash::draw(){
   // printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

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
















