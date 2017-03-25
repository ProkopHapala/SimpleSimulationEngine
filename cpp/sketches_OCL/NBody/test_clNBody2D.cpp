
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "Vec2.h"

#include "testUtils.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"



#define R2SAFE  1.0e-2
#define F2MAX  10.0

//int n = 1024;
//constexpr int n = 32;
//constexpr int n = 256;
constexpr int n = 1024;
float damp = 0.8;
Vec2f pos  [n];
Vec2f vel  [n];
Vec2f force[n];

inline void acum_force( const Vec2f& p1, const Vec2f& p2, Vec2f& f ){
    Vec2f dp; dp.set_sub( p2, p1 );
    float ir2 = 1/( dp.norm2() + R2SAFE );
    //float ir  = sqrt(ir2);
    float ir6 = ir2*ir2*ir2;
    float fr  = ir6 - ir6*ir6;
    f.add_mul( dp, fr );
}

void add_ineraction_forces(){
    for(int i=0; i<n; i++){
        Vec2f f; f.set(0.0);
        for(int j=0; j<n; j++){
            acum_force( pos[i], pos[j], f );
        }
        force[i] = f;
    }
}

void clean_force( ){ for(int i=0; i<n; i++){ force[i].set(0.0); } }

void add_confine_force( const Vec2f& center, float strength ){
    for(int i=0; i<n; i++){ force[i].add_mul( pos[i], strength ); }
}

void add_external_force( const Vec2f& center, float strength ){
    for(int i=0; i<n; i++){
        Vec2f dp; dp.set_sub( pos[i], center );
        float ir2 = 1/(dp.norm2() + R2SAFE);
        float ir  = sqrt(ir2);
        force[i].add_mul( dp, strength*ir2*ir );
    }
}

float move_leap_frog( float dt ){
    double f2max = 0;
    for(int i=0; i<n; i++){
        float f2 = force[i].norm2();
        if(f2>f2max) f2max = f2;
    }
    if(f2max > F2MAX ) dt *= sqrt( sqrt(F2MAX/f2max) );
    for(int i=0; i<n; i++){
        vel[i].mul(damp);
        vel[i].add_mul( force[i], dt );
        pos[i].add_mul( vel[i]  , dt );
    }
    return f2max;
}

// ==================== Declaration

class TestApp_clNBody2D : public AppSDL2OGL {
	public:
	int per_frame = 10;

	//virtual void loop( int nframes );
	virtual void draw   ();
	//virtual void drawHUD();
	TestApp_clNBody2D( int& id, int WIDTH_, int HEIGHT_ );

};

// ==================== Implementation

TestApp_clNBody2D::TestApp_clNBody2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    float span = 5.0;
    for(int i=0; i<n; i++){
        pos[i].set(randf(-span,span),randf(-span,span));
        vel[i].set(0.0);
    }

}

void TestApp_clNBody2D::draw(){
    long t1,t2,t3;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //delay = 100;
    t1=getCPUticks();
    double f2err = 0;
	for(int itr=0; itr<per_frame; itr++){
        clean_force();
        add_ineraction_forces( );
        add_confine_force( {0.0,0.0}, -2.0 );
        //if( LMB ) add_external_force( {mouse_begin_x, mouse_begin_y}, -10.0 );
        f2err = move_leap_frog( 0.01 );
	}
	double fticks = getCPUticks() - t1;
    printf( "%i Ferr= %g T= %g [Mtick] %g [t/op]\n", frameCount*per_frame, sqrt(f2err), fticks*1e-6, fticks/(n*n) );


	glBegin(GL_POINTS);
    for(int i=0; i<n; i++){  glVertex3f( pos[i].x, pos[i].y, 0.0 ); }
    glEnd();


};

// ===================== MAIN

TestApp_clNBody2D * app;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	app = new TestApp_clNBody2D( junk , 800, 600 );
	app->loop( 1000000 );
	return 0;
}
















