
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
//#include "Body.h"

#include "geom3D.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "appliedPhysics.h"
#include "PotentialFlow.h"
#include "DrawField.h"


/*

Magnetic Bottle

https://www.youtube.com/watch?v=Sf1MGTD9xGY
https://physics.stackexchange.com/questions/30261/analogy-between-magnetic-bottle-and-van-allens-radiation-belt



*/

const int nCoillSteps = 8;
//const int nCoils = 5;
//double coil_R[nCoils] =  {  1.0,  2.0, 3.0,  2.0,  1.0 };
//double coil_z[nCoils] =  { -5.0, -3.0, 0.0,  3.0,  5.0 };
//double coil_I[nCoils] =  { 10.0,  0.0, 1.0,  0.0,  10.0 };


//const int nCoils = 4;
//double coil_R[nCoils] =  {  0.5,  3.0,  3.0,  0.5 };
//double coil_z[nCoils] =  { -4.0, -1.5,  1.5,  4.0 };
//double coil_I[nCoils] =  { 10.0,  1.0,  1.0,  10.0 };

//const int nCoils = 6;
//double coil_R[nCoils] =  {  0.5, 2.0,  3.0,  3.0, 2.0, 0.5 };
//double coil_z[nCoils] =  { -4.0,-3.2, -1.5,  1.5, 3.2, 4.0 };
//double coil_I[nCoils] =  { 10.0, 1.0,  1.0,  1.0, 1.0, 10.0*0.0 };

float sc = 1.0;
const int nCoils = 5;
double coil_R[nCoils] =  {  0.5, 2.0,  3.0,  3.0, 2.0 };
double coil_z[nCoils] =  { -4.0,-3.2, -1.5,  1.5, 4.0 };
double coil_I[nCoils] =  { 1.5*sc, 0.125*sc,  0.07*sc,  0.05*sc, 0.25*sc  };

struct Particle{
    double invMass;
    Vec3d pos;
    Vec3d vel;
    void move( double dt, Vec3d force ){
        vel.add_mul( force, dt*invMass );
        pos.add_mul( vel,   dt );
    }
    void makeThermal( double T, double M ){
        double m = M*const_massProton;               //  [kg] proton mass
        invMass=1/m;
        double v = sqrt( T * 3*const_Bonltzman/m );
        printf( "v %g [km/s] E %g [keV]\n", v*1e-3, 0.5*m*sq(v)/(1e+3*const_eV) );
        vel.fromRandomSphereSample();
        vel.mul(v);
    }
};

void drawParticleTrj( int n, double dt, Particle p, VecFieldFunc func ){
    //Vec3d opos=p.pos;
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++ ){
        Vec3d B = func( p.pos );
        Vec3d f;
        f.set_cross( B, p.vel );
        f.add_cross( B, p.vel+f*(const_ElectronCharge*p.invMass*dt) );
        f.mul( const_ElectronCharge*0.5 );
        //printf( "trj[%i] p(%g,%g,%g) v(%g,%g,%g) f(%g,%g,%g) B(%g,%g,%g)\n", i, p.pos.x,p.pos.y,p.pos.z,  p.vel.x,p.vel.y,p.vel.z,  f.x,f.y,f.z, B.x, B.y, B.z );
        //printf( "trj[%i] v %g [km/s]\n", i, p.vel.norm()*1e-3   );
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( f*1e+11, p.pos );
        //glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( B, p.pos );
        //glColor3f(0.0,0.0,0.0); Draw3D::drawVecInPos( p.vel*1e-8, p.pos );
        p.move(dt,f);
        Draw3D::vertex(p.pos);

    }
    glEnd();
}

Vec3d coilField( Vec3d pos, Vec3d h, double R, int n, bool bDraw=false ){
    double dang = M_PI*2/n;
    double dl   = R*dang;
    Vec2d rot=Vec2dX;
    Vec2d drot; drot.fromAngle(dang*0.5);
    if(n<16){
        Vec2d r=rot;
        r.mul_cmplx( drot );
        r.mul_cmplx( drot );
        r.sub(rot);
        dl=r.norm()*R;
    }
    Vec3d B = Vec3dZero;
    Vec3d a,b; h.getSomeOrtho(a,b);

    for(int i=0; i<n; i++){
        Vec3d p;  p .set_lincomb( rot.a, a,  rot.b, b ); p.mul(R); p.add(pos);
        rot.mul_cmplx( drot );
        Vec3d dh; dh.set_lincomb( rot.a, b, -rot.b, a );
        rot.mul_cmplx( drot );
        B.add( ILineFinite( p, dh, dl ) );
        if(bDraw)Draw3D::drawLine( p, p+dh*dl );
    }
    return B;
}

Vec3d coilField( Vec3d pos ){
    Vec3d B=Vec3dZero;
    for(int i=0; i<nCoils; i++){
        B.add_mul( coilField( pos-(Vec3d){0.0,0.0,coil_z[i]}, Vec3dZ, coil_R[i], nCoillSteps, false ), coil_I[i]  );
    }
    //return coilField( (Vec3d){0.0,0.0,0.0} - pos, {0.0,0.0,1.0}, 1.0, 16, false )
    //+      coilField( (Vec3d){0.0,0.0,2.0} - pos, {0.0,0.0,1.0}, 1.0, 16, false )
    return B;
}

void drawCoils(){
    for(int i=0; i<nCoils; i++){
        coilField( (Vec3d){0.0,0.0,coil_z[i]}, Vec3dZ, coil_R[i], nCoillSteps, true );
    }
}


// ============= Application

class TestAppElectromagnetic : public AppSDL2OGL_3D { public:

    int ogl;
    bool bIntegrate=false;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppElectromagnetic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppElectromagnetic::TestAppElectromagnetic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    ogl = Draw::list( );
    glColor3f(0.0,0.0,1.0);
    //plotVecPlane( {30,30}, {0.0,0.0,0.0}, {0.0,0.0,0.25}, {0.25,0.0,0.0},  -1.0, 1.0, coilField );
    plotStreamLinePlane( {20,1}, 500, {0.0,0.0,0.0}, Vec3dX*0.25, Vec3dY*0.25, 0.25, coilField );
    glColor3f(0.0,0.0,0.0);
    drawCoils();
    glColor3f(1.0,0.0,0.0);

    Particle p;
    p.pos.fromRandomBox({-0.5,-0.5,-0.5},{0.5,0.5,0.5});
    p.makeThermal(100e+6,1);
    drawParticleTrj( 5000, 0.5e-8, p, coilField );
    glEndList();
}

void TestAppElectromagnetic::draw(){
    //rintf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    if(bIntegrate){
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0,1.0,1.0,0.1);

        Particle p;
        //p.pos.fromRandomBox({-1.0,-1.0,-1.0},{1.0,1.0,1.0});
        p.pos.fromRandomBox({-0.5,-0.5,-0.5},{0.5,0.5,0.5});
        p.makeThermal(100e+6,1);
        drawParticleTrj( 5000, 0.5e-8, p, coilField );
    }else{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glCallList(ogl);

    }

};


void TestAppElectromagnetic::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_m:  bIntegrate=!bIntegrate; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}


void TestAppElectromagnetic::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppElectromagnetic * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppElectromagnetic( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















