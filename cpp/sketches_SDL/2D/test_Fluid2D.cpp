
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw2D.h"
#include "AppSDL2OGL.h"

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"

#include "Fluid2D.h"
//#include "Fluid2D.cpp"
#include "testUtils.h"
/*
#include "TerrainHydraulics.h"
#include "testUtils.h"


#include "GridIndex2D.h"
#include "Grid2DAlgs.h"
#include "Grid2DAlgs.cpp" // FIXME
#include "SquareRuler.h"
*/

//#include "Grid.h"

// ======================  TestApp


const int nVSource = 4;
//const Vec2i pSources[nVSource] = { 16,32,    48,32,      32,16,     32,48   };
const Vec2i pSources[] = { 16,16+1,    48,16,      16,48-2,     48,48   };
//const Vec2d vSources[nVSource] = { 0.0,-1.0,  0.0,+1.0,  -1.0,0.0,  1.0,0.0 };
//const Vec2d vSources[nVSource] = { 0.0,+1.0,  0.0,-1.0,  -1.0,0.0,  1.0,0.0 };
const Vec2d vSources[] = { 1.0,0.0,  -1.0,0.0,  1.0,0.0,  -1.0,0.0 };


class TestAppFluid2D : public AppSDL2OGL{
	public:
    int job_type = 1;
    int perframe = 3;
    //TerrainHydraulics terrain;
    int shape;
    bool running = true;

    Fluid2D fluid;

    int nParticles;
    Vec2d * particles = 0;

	// ---- function declarations

	virtual void draw   ();
	virtual void drawHUD();
    virtual void eventHandling( const SDL_Event& event );
	//virtual int tileToList( float x0, float y0, float x1, float y1 );

    //void renderMapContent ( float x0, float y0, float scale, float csc, float hsc );
    //double terrain_color( int i );
	//void drawSimplexGrid( int n, float step );

	TestAppFluid2D( int& id, int WIDTH_, int HEIGHT_ );

};

/*
double TestAppFluid2D::terrain_color( int i ){
    float g   = terrain.ground[i];
    float w   = terrain.water [i];
    if( w > g ){
        float c = (1-20*(w-g)); if(c<0) c=0;
        glColor3f(g*g*c,0.2+0.8*g*(1-g)*c,0.5);
    }else{
        glColor3f(g*g,0.2+0.8*(1-g)*g,0);
    }
    return g;
}

void TestAppFluid2D::renderMapContent( float x0, float y0, float scale, float csc, float hsc ){
    //glColor3f( 0.1f,0.1f,0.1f );
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); b.mul(scale);
    //b.set( 0.0d, 1.0 ); b.mul(scale);
    //glDisable(GL_SMOOTH);
    int ii = 0;
    for (int iy=0; iy<terrain.n.y-1; iy+=1){
    //for (int iy=0; iy<terrain.ny-1; iy+=1){
        glBegin( GL_TRIANGLE_STRIP );
        for (int ix=0; ix<terrain.n.x; ix+=1){
        //for (int ix=0; ix<terrain.nx; ix+=1){
            p.set( ix*a.x+iy*b.x + x0, ix*a.y+iy*b.y + y0 );
            terrain_color( ii               ); glVertex3f( p.x    , p.y    , 0 );
            terrain_color( ii + terrain.n.x  ); glVertex3f( p.x+b.x, p.y+b.y, 0 );
            ii++;
        }
        glEnd();
    }
}
*/

TestAppFluid2D::TestAppFluid2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {
    srand(154978);
    //fluid.allocate( {64,64} );
    fluid.allocate( {128,128} );
    //fluid.allocate( {256,256} );
    //fluid.allocate( {128,256} );
    fluid.ndiffuse = 2;

    nParticles = fluid.ntot;
    particles = new Vec2d[nParticles];


    //fluid.source[ fluid.ip2i({16,32}) ] =  1.0;
    //fluid.source[ fluid.ip2i({48,32}) ] = -1.0;

    shape=glGenLists(1);
}

void TestAppFluid2D::draw(){

    //delay = 100;
    glClearColor( 0.5f, 0.5f, 0.5f, 0.0f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glColor4f( 0.5f, 0.5f, 0.5f, 0.05f );
	glColor4f( 1.0f, 1.0f, 1.0f, 0.05f );
	Draw2D::drawRectangle( {-100.0,-100.0},{100.0,100.0},true);
    glDisable(GL_BLEND);

	glDisable( GL_DEPTH_TEST );
    long t0,t1;

    double vsc = 5;
    glBegin(GL_LINES);
    glColor3f(1.0,0.0,0.0);
    for(int ii=0; ii<nVSource; ii++){
        const Vec2i& ip = pSources[ii];
        const Vec2d& v  = vSources[ii];
        int i = fluid.ip2i( ip );
        //printf( "%i (%i,%i) (%f,%f) \n", ii, ip.x, ip.y, v.x, v.y  );
        fluid.vx[i] = v.x;
        fluid.vy[i] = v.y;
        Vec2d p = (Vec2d){ ip.x*0.1,ip.y*0.1 };
        glVertex3f(p.x,p.y,0);
        glVertex3f(p.x+fluid.vx[i]*vsc,p.y+fluid.vy[i]*vsc,0);
    }
    glEnd();

    double dt = 0.1;

    glColor3f(1.0,0.0,1.0);
    t0 = getCPUticks();
    fluid.fluidStep( dt );
    t0 = getCPUticks()-t0;

    //glBegin(GL_LINES);
    glBegin(GL_POINTS);
    glColor3f(1.0,0.0,0.0);
    float psc = 10000.0;
    for (int iy=1; iy<fluid.n.y-1; iy++){ for (int ix=1; ix<fluid.n.x-1; ix++){
        int i = fluid.ip2i( {ix,iy} );
        Vec2d p = (Vec2d){ ix*0.1,iy*0.1 };
        float c = fluid.p[i]*psc;
        //if( c<0 ){ glColor3f(1.0,1.0+c,1.0+c); }else{ glColor3f(1-c,1-c,1.0); }
        //glVertex3f(p.x,p.y,0);
        //glVertex3f(p.x+fluid.vx[i]*vsc,p.y+fluid.vy[i]*vsc,0);
    } }
    glEnd();

    glBegin(GL_POINTS);
    glColor3f(0.0,0.0,0.0);
    for(int i=0; i<nParticles; i++){
        Vec2d& p = particles[i];
        if( randf()<0.01 ){ p.set(randf(0.0,fluid.n.x),randf(0.0,fluid.n.y)); };
        double vx = fluid.interpBilinear( p, fluid.vx );
        double vy = fluid.interpBilinear( p, fluid.vy );
        p.add_mul( (Vec2d){vx,vy}, dt*50 );
        glVertex3f( p.x*0.1, p.y*0.1,0.0 );
    }
    glEnd();

    printf( " %f Mticks %f op/pix \n", t0*1.0e-6 ,((double)t0)/fluid.ntot );
};

void TestAppFluid2D::drawHUD(){}


void TestAppFluid2D::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
               // case SDLK_r:  terrain.initErrosion( 0.8 ); running=true;  break;
            }
            break;
         /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //paintSimplex( mouse_begin_x, mouse_begin_y );
                    mouse_left  = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = true;
                    //eraseSimplex( mouse_begin_x, mouse_begin_y );
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    mouse_left = false;
                    break;
                case SDL_BUTTON_RIGHT:
                    mouse_right = false;
                    break;
            }
            break;
        */
    };
    AppSDL2OGL::eventHandling( event );
};

// ===================== MAIN

TestAppFluid2D * testApp;

int main(int argc, char *argv[]){

	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestAppFluid2D( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















