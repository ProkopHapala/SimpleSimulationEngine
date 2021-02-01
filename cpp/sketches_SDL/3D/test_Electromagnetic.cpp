
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

#include "Ruler2DFast.h"


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


//float sc = 0.8;
//const int nCoils = 5;
//double coil_R[nCoils] =  { 0.25, 2.0,  3.0,  3.0, 2.0 };
//double coil_z[nCoils] =  { -4.0,-3.2, -1.3,  1.3, 4.0 };
//double coil_I[nCoils] =  { 2.5*sc, 0.125*sc,  0.07*sc,  0.05*sc, 0.40*sc  };

float sc = 0.8;
const int nCoils = 6;
double coil_R[nCoils] =  { 0.25, 2.0,  3.0,  3.0, 2.0, 1.0 };
double coil_z[nCoils] =  { -4.0,-3.2, -1.3,  1.3, 4.0, 5.0 };
double coil_I[nCoils] =  { 2.5*sc, 0.125*sc,  0.07*sc,  0.05*sc, 0.20*sc, 0.60*sc  };

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

class GridFF2d : public Ruler2DFast{ public:
    Vec2d*  data;
    double* hits;

    GridFF2d( Vec2i n_, Vec2d pmin, Vec2d pmax ){
        Vec2d span = pmax-pmin;
        setup( pmin, {span.x/n_.x,span.y/n_.y} );
        n_.x++; n_.y++;
        setN( n_ );
        _realloc(data, ntot);
        _realloc(hits, ntot);
        for(int i=0; i<ntot; i++){ hits[i]=0; };
    }

    Vec2d interpolate( Vec2d p )const{
        Vec2d d;
        Vec2i ipos;
        pos2index(p, d, ipos);
        if( (ipos.x>0)&&(ipos.x<n.x-1) && (ipos.y>0)&&(ipos.y<n.y-1) ){
            int i = ip2i(ipos);
            double mx=1-d.x;
            double my=1-d.y;
            Vec2d f;
            f.set_mul( data[i      ], my*mx  );
            f.add_mul( data[i    +1], my*d.x );
            f.add_mul( data[i+n.x  ],d.y*mx  );
            f.add_mul( data[i+n.x+1],d.y*d.x );
            return f;
        }
        return Vec2dZero;
    }

    void acumHit( Vec2d p, double w){
        Vec2d d;
        Vec2i ipos;
        pos2index(p, d, ipos);
        if( (ipos.x>0)&&(ipos.x<n.x-1) && (ipos.y>0)&&(ipos.y<n.y-1) ){
            int i = ip2i(ipos);
            double mx=1-d.x;
            double my=1-d.y;
            hits[i      ]+=w* my*mx  ;
            hits[i    +1]+=w* my*d.x ;
            hits[i+n.x  ]+=w*d.y*mx  ;
            hits[i+n.x+1]+=w*d.y*d.x ;
        }
    }

};

GridFF2d gridFF( {121,41}, {-5.,0.}, {7.,4.} );

void drawParticleTrj( int n, double dt, Particle p, VecFieldFunc func ){
    //Vec3d opos=p.pos;
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++ ){
        Vec3d B = func( p.pos );
        Vec3d f;
        f.set_cross( B, p.vel );
        f.add_cross( B, p.vel+f*(const_ElectronCharge*p.invMass*dt) );
        f.mul( const_ElectronCharge*0.5 );

        gridFF.acumHit( { p.pos.z, sqrt(sq(p.pos.x)+sq(p.pos.y)) }, 1);
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

Vec3d coilFieldGrid( Vec3d pos ){
    Vec2d p2d;
    p2d.x     = pos.z;
    p2d.y     = sqrt( pos.x*pos.x + pos.y*pos.y );
    double invr = 1/p2d.y;
    Vec2d f2d = gridFF.interpolate( p2d );
    return (Vec3d){ f2d.y*pos.x*invr, f2d.y*pos.y*invr, f2d.x };
}

void drawCoils(){
    for(int i=0; i<nCoils; i++){
        coilField( (Vec3d){0.0,0.0,coil_z[i]}, Vec3dZ, coil_R[i], nCoillSteps, true );
    }
}

void prepareFFgrid( GridFF2d& grid ){
    for(int iy=0; iy<grid.n.y; iy++){
        for(int ix=0; ix<grid.n.x; ix++){
            Vec2d p;
            grid.index2pos( {ix,iy},{0,0}, p );
            Vec3d B = coilField( {0,p.y,p.x} );
            int   i = grid.ip2i( {ix,iy} );
            grid.data[i].set( B.z, B.y );
        }
    }
};

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

    prepareFFgrid( gridFF );

    ogl = Draw::list( );
    glColor3f(0.0,0.0,1.0);
    //plotVecPlane( {30,30}, {0.0,0.0,0.0}, {0.0,0.0,0.25}, {0.25,0.0,0.0},  -1.0, 1.0, coilFieldGrid );
    //plotStreamLinePlane( {20,1}, 500, {0.0,0.0,0.0}, Vec3dX*0.25, Vec3dY*0.25, 0.25, coilField );
    plotStreamLinePlane( {20,1}, 500, {0.0,0.0,0.0}, Vec3dX*0.25, Vec3dY*0.25, 0.25, coilFieldGrid );
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
        long t = getCPUticks();
        int nstep=5000;
        //drawParticleTrj( 5000, 0.5e-8, p, coilField );
        drawParticleTrj( nstep, 0.5e-8, p, coilFieldGrid );
        t= getCPUticks()-t;
        printf( " %g CPUticks/step \n", t/(double)nstep );
    }else{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
        glCallList(ogl);
        double logHist[gridFF.ntot];
        for(int i=0; i<gridFF.ntot; i++){ logHist[i]=log(gridFF.hits[i]+1); };
        Draw3D::drawScalarGrid( gridFF.n, {0,gridFF.pos0.y,gridFF.pos0.x}, {0.0,0.0,gridFF.step.x}, {-gridFF.step.y,0.0,0.0}, logHist,  0.0, 12.0 );

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
















