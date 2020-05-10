
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

#include "SDL_utils.h"
#include "GUI.h"

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


    GUI gui;
    GUIPanel panel_dt;
    GUIPanel panel_diff;
    GUIPanel panel_visc;
    GUIPanel panel_ndiff;
    GUIPanel panel_npress;

    GUIPanel panel_CPUtime;

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


/*
void testClamp(double a, double amin, double amax){
    //printf(" %g in (%g ... %g) -> %g \n", a, amin, amax, _clamp(a, amin, amax) );
}
*/


TestAppFluid2D::TestAppFluid2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL( id, WIDTH_, HEIGHT_ ) {

    GUI_fontTex = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    //((GUIPanel*)gui.addPanel( new GUIPanel( "dt [Rad]", 5,5,105,35, true, true ) )) -> setRange(0.0, 0.2);
    //    ->command = &command_example;
    panel_dt.initPanel( "dt", 5,5,90,35 );
    panel_dt.setRange(0.0, 0.2);
    panel_dt.value = 0.05;
    gui.addPanel( &panel_dt );

    panel_visc.initPanel( "viscosity", 100,5,190,35 );
    panel_visc.setRange(0.0, 0.001);
    panel_visc.value = 0.0001*4;
    gui.addPanel( &panel_visc );

    panel_diff.initPanel( "diffuse", 200,5,290,35 );
    panel_diff.setRange(0.0, 5.0);
    panel_diff.value = 3.0;
    gui.addPanel( &panel_diff );

    panel_ndiff.initPanel( "ndiffuse", 300,5,390,35 );
    panel_ndiff.setRange(0, 10);
    panel_ndiff.value = 5;
    panel_ndiff.isInt=true;
    gui.addPanel( &panel_ndiff );

    panel_npress.initPanel( "npressure", 400,5,490,35 );
    panel_npress.setRange(0, 10);
    panel_npress.value = 5;
    panel_npress.isInt=true;
    gui.addPanel( &panel_npress );

    panel_CPUtime.initPanel( "CPUtime [tick/pix]", 5,40,90,70 );
    panel_CPUtime.setRange(0, 1000);
    gui.addPanel( &panel_CPUtime );

    //int ndiffuse  = 5;
    //int npressure = 5;
    //double	visc = 0.0001*4;
	//double	diff = 3.0;

    srand(154978);
    //fluid.allocate( {64,64} );
    fluid.allocate( {128,128} );
    //fluid.allocate( {256,256} );
    //fluid.allocate( {128,256} );
    fluid.ndiffuse = 2;

    nParticles = fluid.ntot;
    particles = new Vec2d[nParticles];

    /*
    testClamp(  0.6,  0.5, 1.5 );
    testClamp(  0.6,  0.5, 1.5 );
    testClamp(  2.3,  0.5, 1.5 );
    testClamp( -0.2, -0.7,  1.5 );
    testClamp( -0.7, -0.5,  0.7 );
    testClamp( -0.1, -0.5, -0.2 );
    */
    //fluid.interpBilinear( {16.3,13.6}, fluid.vx );
    //exit(0);

    //fluid.source[ fluid.ip2i({16,32}) ] =  1.0;
    //fluid.source[ fluid.ip2i({48,32}) ] = -1.0;

    shape=glGenLists(1);
}

void TestAppFluid2D::draw(){

    //delay = 100;
    glClearColor( 1.0f, 1.0f, 1.0f, 0.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glDisable(GL_LIGHTING);
	/*
	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glColor4f( 0.5f, 0.5f, 0.5f, 0.05f );
	//glColor4f( 1.0f, 1.0f, 1.0f, 0.05f );
	glColor4f( 1.0f, 1.0f, 1.0f, 0.2f );
	Draw2D::drawRectangle( {-100.0,-100.0},{100.0,100.0},true);
    glDisable(GL_BLEND);
    */


    if(frameCount==1){  camX0=fluid.n.x*0.05; camY0=fluid.n.y*0.05; };

	glDisable( GL_DEPTH_TEST );

    double vsc = 5;

    glBegin(GL_LINES);
    glColor3f(0.0,1.0,0.0);
    for(int ii=0; ii<nVSource; ii++){
        const Vec2i& ip = pSources[ii];
        const Vec2d& v  = vSources[ii];
        int i = fluid.ip2i( ip );
        //printf( "%i (%i,%i) (%f,%f) \n", ii, ip.x, ip.y, v.x, v.y  );
        fluid.vx[i] = v.x;
        fluid.vy[i] = v.y;
        fluid.vx_[i] = v.x;
        fluid.vy_[i] = v.y;
        Vec2d p = (Vec2d){ ip.x*0.1,ip.y*0.1 };
        glVertex3f(p.x,p.y,0);
        glVertex3f(p.x+fluid.vx[i]*vsc,p.y+fluid.vy[i]*vsc,0);
    }
    glEnd();

    double dt  = panel_dt.value;
    fluid.diff = panel_diff.value;
    fluid.visc = panel_visc.value;
    fluid.npressure=(int)panel_npress.value;
    fluid.ndiffuse =(int)panel_ndiff .value;
    //double dt = 0.05;

    long t0;
    glColor3f(1.0,0.0,1.0);
    t0 = getCPUticks();
    fluid.fluidStep( dt );
    t0 = getCPUticks()-t0;

    panel_CPUtime.value = ((double)t0)/fluid.ntot;
    panel_CPUtime.redraw=true;

    //int i0 = 128*50+50;
    //printf( " vx,vy %g %g  \n", fluid.vx[i0], fluid.vx[i0] );

    //SDL_Delay(500);

    glBegin(GL_LINES);
    //glBegin(GL_POINTS);
    glColor3f(1.0,0.0,0.0);
    float psc = 10000.0;
    vsc = 1.0;
    for (int iy=1; iy<fluid.n.y-1; iy++){ for (int ix=1; ix<fluid.n.x-1; ix++){
        int i = fluid.ip2i( {ix,iy} );
        Vec2d p = (Vec2d){ ix*0.1,iy*0.1 };
        float c = fluid.p[i]*psc;
        //if( c<0 ){ glColor3f(1.0,1.0+c,1.0+c); }else{ glColor3f(1-c,1-c,1.0); }

        glColor3f(1.0,0.0,0.0);
        glVertex3f(p.x,p.y,0);
        glVertex3f( p.x+fluid.vx[i]*vsc, p.y+fluid.vy[i]*vsc, 0);

        /*
        double vx = fluid.interpBilinear( (Vec2d){ix,iy}, fluid.vx );
        double vy = fluid.interpBilinear( (Vec2d){ix,iy}, fluid.vy );
        glColor3f(0.0,0.0,1.0);
        glVertex3f(p.x,p.y,0);
        glVertex3f( p.x+vx*vsc, p.y+vy*vsc, 0);
        */

    } }
    glEnd();

    /*
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
    */


    /*
    glBegin(GL_LINES);
    glColor3f(0.0,0.0,0.0);
    double L = 1.5;
    for(int i=0; i<nParticles/2; i++){
        int i2 = i<<1;
        Vec2d& p1 = particles[i2+0];
        Vec2d& p2 = particles[i2+1];
        if( (randf()<0.0001)||(frameCount==1) ){
            double ang = randf(0.0,M_PI*2);
            double dx = cos(ang)*L*0.5;
            double dy = sin(ang)*L*0.5;
            double x = randf(0.0,fluid.n.x);
            double y = randf(0.0,fluid.n.y);
            p1.set(x+dx,y+dy);
            p2.set(x-dx,y-dy);
        }else{
            double vx,vy;
            vx = fluid.interpBilinear( p1, fluid.vx );
            vy = fluid.interpBilinear( p1, fluid.vy );
            p1.add_mul( (Vec2d){vx,vy}, dt*50 );
            vx = fluid.interpBilinear( p2, fluid.vx );
            vy = fluid.interpBilinear( p2, fluid.vy );
            p2.add_mul( (Vec2d){vx,vy}, dt*50 );
            Vec2d d = p2-p1;
            double l = d.norm();
            double c = (l-L*0.5)/l;
            p1.add_mul(d, c);
            p2.add_mul(d,-c);
            //if(i==0) printf("%f %f \n",l,c);
        };
        glVertex3f( p1.x*0.1, p1.y*0.1,0.0 );
        glVertex3f( p2.x*0.1, p2.y*0.1,0.0 );
    }
    glEnd();
    */

    glColor3f(0.0,0.0,0.0);
    double L = 1.5;
    int nl   = 10;
    for(int i=0; i<nParticles/nl; i++){
        if( (randf()<0.0005)||(frameCount==1) ){
            double ang = randf(0.0,M_PI*2);
            double dx  = cos(ang)*L*0.5;
            double dy  = sin(ang)*L*0.5;
            double x   = randf(0.0,fluid.n.x);
            double y   = randf(0.0,fluid.n.y);
            for(int j=0; j<nl; j++){
                particles[i*nl+j].set(x+dx*j,y+dy*j);
            }
        }else{
            double vx,vy;
            glBegin(GL_LINE_STRIP);
            Vec2d* op=0;
            for(int j=0; j<nl; j++){
                Vec2d& p = particles[i*nl+j];

                vx = fluid.interpBilinear( p, fluid.vx );
                vy = fluid.interpBilinear( p, fluid.vy );
                p.add_mul( (Vec2d){vx,vy}, dt*50 );
                if( j>0 ){
                    Vec2d  d = p-(*op);
                    double l = d.norm();
                    double c = (l-L*0.5)/l;
                    p    .add_mul(d,-c*0.5);
                    (*op).add_mul(d,+c*0.5);
                }

                glVertex3f( p.x*0.1, p.y*0.1,0.0 );
                op=&p;
            }
            glEnd();
            //if(i==0) printf("%f %f \n",l,c);
        };
    }

    //printf( " %f Mticks %f op/pix \n", t0*1.0e-6 ,((double)t0)/fluid.ntot );
};

void TestAppFluid2D::drawHUD(){
    glDisable( GL_LIGHTING );
    glDisable(GL_DEPTH_TEST);


    Draw2D::drawText( "AHOJ !!!!", 0, {10, 100}, 0.0, GUI_fontTex, fontSizeDef );

    gui.draw();
}


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
    gui.onEvent( mouseX, HEIGHT-mouseY, event );
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
















