
/// @file @brief  This is a physics simulation of charged particles in a magnetic field. It models a thin, hot plasma where ions are trapped within a "magnetic bottle". The simulation calculates the aggregate electric and magnetic fields using a Poisson solver (`Fourier.h`) and `PotentialFlow.h`, then integrates the motion of the ions. The description suggests pressing the 'm' key to run or control the simulation.
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

#include "VecN.h"

#include "appliedPhysics.h"
#include "PotentialFlow.h"
#include "DrawField.h"

#include "Ruler2DFast.h"

#include "Fourier.h"


/*

Magnetic Bottle

https://www.youtube.com/watch?v=Sf1MGTD9xGY
https://physics.stackexchange.com/questions/30261/analogy-between-magnetic-bottle-and-van-allens-radiation-belt



*/

bool   bElFieldFromGrid = false;
double scaleFgrid       = 1;
double scaleFwire       = 1;

double Rwire = 1.0;

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

const int nWire = 2;
double wire_y[nCoils] =  { -0.0,-0.5 };
double wire_z[nCoils] =  { -2.0, 2.0 };
double wire_Q[nCoils] =  {  1000000.0, -1000000.0 };

struct Particle{
    double Q       = 1.0; // charge
    double invMass = 1.0;
    Vec3d pos;
    Vec3d vel;
    void move( double dt, Vec3d force ){
        vel.add_mul( force, dt*invMass*Q );
        pos.add_mul( vel,   dt );
    }
    void makeThermal( double T, double M ){
        double m = M*const_massProton;               //  [kg] proton mass
        invMass=1/m;
        //Q = 1.0;
        double v = sqrt( T * 3*const_Bonltzman/m );
        printf( "v %g [km/s] E %g [keV]\n", v*1e-3, 0.5*m*sq(v)/(1e+3*const_eV) );
        vel.fromRandomSphereSample();
        vel.mul(v);
    }
};

class GridFF2d : public Ruler2DFast{ public:
    Vec2d*  force=0;
    double* dens=0;
    double* Vpot=0;
    double* Vpot_=0;

    GridFF2d( Vec2i n_, Vec2d pmin, Vec2d pmax ){
        Vec2d span = pmax-pmin;
        setup( pmin, {span.x/n_.x,span.y/n_.y} );
        n_.x++; n_.y++;
        setN( n_ );
        _realloc(force, ntot);
        _realloc(dens, ntot);
        for(int i=0; i<ntot; i++){ dens[i]=0; };
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
            f.set_mul( force[i      ], my*mx  );
            f.add_mul( force[i    +1], my*d.x );
            f.add_mul( force[i+n.x  ],d.y*mx  );
            f.add_mul( force[i+n.x+1],d.y*d.x );
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
            dens[i      ]+=w* my*mx  ;
            dens[i    +1]+=w* my*d.x ;
            dens[i+n.x  ]+=w*d.y*mx  ;
            dens[i+n.x+1]+=w*d.y*d.x ;
        }
    }

    void poissonStep( double cRho , double cR ){
        for(int iy=1; iy<n.y-1; iy++ ){
            for(int ix=1; ix<n.x-1; ix++ ){
                int i = ip2i({ix,iy});
                // V <= V + a * R
                // R  = 0.25 * (  rho[0,0]*h^2/eps  V[+,0] + V[-,0] + V[0,+] + V[0,-] )
                // see file:///home/prokop/Dropbox/KnowDev/NumMath/PDE/Solving_Generalized_Poisson_Equation_using_FDM.pdf
                double R = Vpot[i-1]+Vpot[i+1] + Vpot[i-n.x] + Vpot[i+n.x];
                R       += cRho*dens[i];
                R*=0.25;
                //R-=Vpot[i];
                Vpot_[i] = Vpot[i] + R*cR;
            }
        }
        _swap(Vpot,Vpot_);
    }

    void potential2force(){
        for(int iy=0; iy<n.y; iy++ ){
            for(int ix=0; ix<n.x; ix++ ){
                int i = ip2i({ix,iy});
                // Calculate force.y
                if(iy < n.y - 1){ // Not the last row
                    force[i].y = (Vpot[i+n.x]-Vpot[i])*invStep.y;
                } else { // Last row, no "below" neighbor
                    force[i].y = 0;
                }
                // Calculate force.x
                if(ix < n.x - 1){ // Not the last column
                    force[i].x = (Vpot[i+1]-Vpot[i])*invStep.x;
                } else { // Last column, no "right" neighbor
                    force[i].x = 0;
                }
            }
        }
    }

    void solvePoisson( int nstep, double eps, bool bMakeFF=true){
        _realloc(Vpot ,ntot);
        _realloc(Vpot_,ntot);
        for(int i=0;i<ntot;i++){ Vpot[i]=0; Vpot_[i]=0; }
        double h    = sqrt(step.x*step.y);
        double cRho =  1e-6 *  h*h/eps;
        double t    = cos(M_PI/n.x) + cos(M_PI/n.y);
        //double cR   = ( 8 - sqrt( 64 - 16*t*t ) )/(t*t);
        double cR   = 1.0;
        // see file:///home/prokop/Dropbox/KnowDev/NumMath/PDE/Solving_Generalized_Poisson_Equation_using_FDM.pdf
        printf( "solvePoisson t %g cR %g cRho %g \n", t, cR, cRho );
        for(int i=0; i<nstep;  i++){
            poissonStep(cRho,cR);
        }
        double vmin=+1e+300,vmax=-1e+300;
        for(int i=0;i<ntot;i++){ vmin=_min(vmin,Vpot[i]); vmax=_max(vmax,Vpot[i]); }
        printf( "solvePoisson vmin %g vmax %g \n", vmin, vmax );
        if(bMakeFF)potential2force();
    }

    void solvePoissonFFT(bool bMakeFF=true){
        _realloc(Vpot ,ntot);
        _realloc(Vpot_,ntot);
        int nnx = log2(n.x);
        int nny = log2(n.y);
        printf( "nn(%i,%i) \n", nnx, nny );
        FFT_2D(nnx, nny, dens,  1, Vpot_ );
        /*
        double C = 1;
        int nhx=n.x/2;
        int nhy=n.y/2;
        for(int iy=0; iy<n.y; iy++ ){
            for(int ix=0; ix<n.x; ix++ ){
                int i = ip2i({ix,iy});
                int ix_=ix-nhx;
                int iy_=iy-nhy;
                Vpot_[i] *= C/( ix_*ix_ + iy_*iy_ );
            }
        }
        FFT_2D(nnx, nny, Vpot, -1, Vpot  );
        */
        if(bMakeFF)potential2force();
    }

    void solvePoisson(){
        solvePoisson(20,1.0);
        //solvePoissonFFT();
    }

    void dampDens( double bmix ) { // ToDo : we may use VecN::mult()
        double damp = 1-bmix;
        double vmin =  1e+300;
        double vmax = -1e+300;

        for(int i=0;i<ntot;i++){ dens[i]*=damp;
            double v = dens[i];  _setmin(vmin,v); _setmax(vmax,v);
        }

        /*
        // DEBUG - replace by simple fixed test density distribution
        for(int iy=0; iy<n.y; iy++ ){
            for(int ix=0; ix<n.x; ix++ ){
                int i = ip2i({ix,iy});
                if( (ix<260)&&(ix>250)  &&  (iy<260)&&(iy>250) ){ dens[i]=1000; }else{ dens[i]=0; }
            }
        }
        */

        printf( "dampDens() vmin %g vmax %g \n", vmin, vmax );
    }

};

//GridFF2d gridFF( {121,41}, {-5.,0.}, {7.,4.} );

//GridFF2d gridFF( {128,64}, {-5.,-3.0}, {7.7,3.3} );
GridFF2d gridFF( {512,512}, {-25.0,-25.0}, {25.0,25.0} );

void traceParticleTrj( int n, double dt, Particle p, VecFieldFunc func, bool bDraw=true ){
    //Vec3d opos=p.pos;
    //const bool bCylindrical = true;
    //const bool bMagnetic    = true;
    const bool bCylindrical = false;
    const bool bMagnetic    = false;
    if(bDraw)glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++ ){
        Vec3d f;
        if(bMagnetic){ // Magnetic Coils
            Vec3d B = func( p.pos );
            f.set_cross( B, p.vel );
            f.add_cross( B, p.vel+f*(const_ElectronCharge*p.invMass*dt) ); // Some drift force - have to look on plasmahydrodynamics
            f.mul( const_ElectronCharge*0.5 );
        }else{  // Electro-Static
            //f.set_mul( func( p.pos ), const_ElectronCharge );
            f.set_mul( func( p.pos ), 1e+9 );
            //f.mul( q );
        }
        Vec2d p2d;
        if(bCylindrical){
            p2d = { p.pos.z,  sqrt(sq(p.pos.x)+sq(p.pos.y)) };
        }else{
            p2d = { p.pos.z, -p.pos.x };
        }
        gridFF.acumHit( p2d, p.Q );
        //printf( "trj[%i] p(%g,%g,%g) v(%g,%g,%g) f(%g,%g,%g) B(%g,%g,%g)\n", i, p.pos.x,p.pos.y,p.pos.z,  p.vel.x,p.vel.y,p.vel.z,  f.x,f.y,f.z, B.x, B.y, B.z );
        //printf( "trj[%i] v %g [km/s]\n", i, p.vel.norm()*1e-3   );
        //glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( f*1e+11, p.pos );
        //glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( B, p.pos );
        //glColor3f(0.0,0.0,0.0); Draw3D::drawVecInPos( p.vel*1e-8, p.pos );


        float dtfac = 1. + 0.5*randf(-1,1);
        p.move(dt*dtfac,f);
        if(bDraw){
            Draw3D::vertex(p.pos);
            //Draw3D::vertex(p.pos + f*dt*10000.0);
            //Draw3D::vertex(p.pos);
        }
        //if(i==0){ printf( "f(%g,%g,%g) dt %g m %g \n", f.x,f.y,f.z,  dt, 1/p.invMass ); }

    }
    if(bDraw)glEnd();
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


Vec3d elecField( Vec3d pos ){
    Vec3d E=Vec3dZero;
    double R2wire = Rwire*Rwire;
    for(int i=0; i<nWire; i++){
        Vec3d dp; dp.set_sub(pos, (Vec3d){wire_y[i], pos.y, wire_z[i] } );
        E.add_mul( dp, scaleFwire*wire_Q[i]/( dp.norm2() + R2wire ) );

        if(bElFieldFromGrid){
            Vec2d p2d;
            p2d.x     = pos.z;
            p2d.y     = pos.x;
            Vec2d f2d = gridFF.interpolate( p2d ) * -scaleFgrid;
            E.add( f2d.y, 0, f2d.x );
            //E.add( f2d.y, 0, f2d.x );
        }

    }
    return E;
}


void drawCoils(){
    for(int i=0; i<nCoils; i++){
        coilField( (Vec3d){0.0,0.0,coil_z[i]}, Vec3dZ, coil_R[i], nCoillSteps, true );
    }
}

void drawWires(double L){
    for(int i=0; i<nWire; i++){
        Draw3D::drawLine((Vec3d){wire_y[i],-L,wire_z[i]},(Vec3d){wire_y[i],L,wire_z[i]});
    }
}

void prepareFFgrid( GridFF2d& grid ){
    for(int iy=0; iy<grid.n.y; iy++){
        for(int ix=0; ix<grid.n.x; ix++){
            Vec2d p;
            grid.index2pos( {ix,iy},{0,0}, p );
            Vec3d B = coilField( {0,p.y,p.x} );
            int   i = grid.ip2i( {ix,iy} );
            grid.force[i].set( B.z, B.y );
        }
    }
};

// ============= Application

class TestAppElectromagnetic : public AppSDL2OGL_3D { public:

    int ogl;
    bool bIntegrate=false;
    bool bDrawTrj  = true;
    bool bPlotGridForceField = false;
    int perFrame   = 100;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppElectromagnetic( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppElectromagnetic::TestAppElectromagnetic( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    prepareFFgrid( gridFF );

    VecN::set( gridFF.ntot, 0.0, gridFF.dens );



    // test FFT
    int nn = 10;
    int n  = 1 << nn;
    Vec2d data[n];
    double f =0;
    double df=0.2;
    double kf = (2*M_PI/n);
    for(int i=0; i<n; i++){
        double w = ( 1-cos(i*2*kf) );
        //data[i].set( cos( i*kf*16 )*w, sin( i*kf*8 )*w*0 );
        //data[i].set( cos( i*(0.04 +0.08*sin(i*0.04) )  ), 0.0 );

        //f+=randf(-df,df);
        //data[i].set(f,0);
        //data[i].set( cos(i*0.04) + cos(i*0.3) + cos(i*1.0)  , 0 );

        data[i].set( 0, 0 );

        //data[i].set( cos(i*0.04) + cos(i*0.3) + cos(i*1.0)  , sin(i*0.04) + sin(i*0.3) + sin(i*1.0) );
    }
    data[n/4].set(100.0,0);

    FFT( (double*)data, n/2,  1);
    for(int i=0; i<(n/4); i++){
        //double k = i;
        double k = i+1;
        double f=0;
        if(k!=0){ f=5.0/(k*k); }
        data[    i  ].mul(f);
        data[n/2-i-1].mul(f); // it is symmetric
    } // poisson kernell
    FFT( (double*)data, n/2, -1);

    ogl = Draw::list( );
    for(int i=0; i<n; i++){
        double x=i*0.1;
        glColor3f(1,0,0); Draw3D::drawLine( (Vec3d){x,0,0}, (Vec3d){x,data[i].x,0} );
        glColor3f(0,0,1); Draw3D::drawLine( (Vec3d){x,0,0}, (Vec3d){x,data[i].y,0} );

        glColor3f(0,0,0); Draw3D::drawLine( (Vec3d){x,0,0}, (Vec3d){x,0.1,0} );
    }
    glEndList();




    /*
    ogl = Draw::list( );
    glColor3f(0.0,0.0,1.0);
    //plotVecPlane( {30,30}, {0.0,0.0,0.0}, {0.0,0.0,0.25}, {0.25,0.0,0.0},  -1.0, 1.0, coilFieldGrid );
    //plotStreamLinePlane( {20,1}, 500, {0.0,0.0,0.0}, Vec3dX*0.25, Vec3dY*0.25, 0.25, coilField );
    //plotStreamLinePlane( {20,1}, 500, {0.0,0.0,0.0}, Vec3dX*0.25, Vec3dY*0.25, 0.25, coilFieldGrid );
    glColor3f(0.0,0.5,0.0);
    //drawCoils();
    drawWires(3.);
    glColor3f(1.0,0.0,0.0);
    //Particle p;
    //p.pos.fromRandomBox({-0.5,-0.5,-0.5},{0.5,0.5,0.5});
    //p.makeThermal(100e+6,1);
    //drawParticleTrj( 5000, 0.5e-8, p, coilField );
    glEndList();
    */

    //qCamera = Quat4fFront;
    //qCamera = Quat4fLeft;
    qCamera = Quat4fTop;


}

void TestAppElectromagnetic::draw(){
    //rintf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	//glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );




    if(bIntegrate){
        glDisable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);



        glCallList(ogl);
        //return;


        bElFieldFromGrid = false;
        scaleFgrid=1e+6 * 0;
        scaleFwire=1;
        //perFrame = 0;
        perFrame = 1000;
        int nstep=100;
        long t = getCPUticks();
        bDrawTrj = true;
        for(int i=0; i<perFrame; i++){
            Particle p;
            //p.pos.fromRandomBox({-1.0,-1.0,-1.0},{1.0,1.0,1.0});
            //p.pos.fromRandomBox({-0.5,-0.5,-0.5},{0.5,0.5,0.5});
            p.pos.fromRandomBox({-2.5,-0.5,-4.5},{+2.5,0.5,-3.5});
            //p.makeThermal(100e+6,1);
            p.vel.fromRandomBox({-0.1,-0.1,-0.1},{0.1,0.1,0.1});
            p.vel.z = 1.0;
            p.vel.mul( 1e+8 );
            //p.vel.mul( 5e+7 );
            //p.vel.mul( 2e+7 );
            p.Q = 1.0;
            if((frameCount%2)==0){ p.Q*=-1; p.invMass*=10.0;  };


            //int nstep=5000;
            //drawParticleTrj( 5000, 0.5e-8, p, coilField );
            //drawParticleTrj( nstep, 0.5e-8, p, coilFieldGrid );
            //glColor4f(1.0,1.0,1.0,0.1);
            //glColor4f(0.0,0.0,0.0,0.1);
            if(p.Q<0){ glColor4f(0.0,0.0,1.0,0.1); }else{ glColor4f(1.0,0.0,0.0,0.1); };

            if(i>20)bDrawTrj=false;
            traceParticleTrj( nstep, 0.5e-8, p, elecField, bDrawTrj );
        }

        gridFF.dampDens(0.01);
        if((frameCount%5)==0){
            //gridFF.solvePoisson(20,1.0);
            gridFF.solvePoisson();
        }




        t= getCPUticks()-t;
        printf( " %g CPUticks/step \n", t/(double)nstep );
    }else{
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );


        glColor3f(0.,0.,0.);
        Vec3d p{mouse_begin_x, 0.0, mouse_begin_y};
        Vec3d f = elecField( p );
        Draw3D::drawVecInPos (f*0.00001,p);
        Draw3D::drawPointCross(p, 0.5 );
        //printf( "p(%g,%g,%g) f(%g,%g,%g)\n", p.x,p.y,p.z,   f.x,f.y,f.z );


        glCallList(ogl);
        double logHist[gridFF.ntot];
        //for(int i=0; i<gridFF.ntot; i++){ logHist[i]=log(gridFF.dens[i]+1); };
        //Draw3D::drawScalarGrid( gridFF.n, {0,gridFF.pos0.y,gridFF.pos0.x}, {0.0,0.0,gridFF.step.x}, {-gridFF.step.y,0.0,0.0}, logHist,  0.0, 12.0 );
        //DEBUG
        //for(int i=0; i<gridFF.ntot; i++){ logHist[i]=sinh(gridFF.dens[i]); };
        //Draw3D::drawScalarGrid( gridFF.n, {-gridFF.pos0.y,0.0,gridFF.pos0.x}, {0.0,0.0,gridFF.step.x}, {-gridFF.step.y,0.0,0.0}, gridFF.dens,  -5.0, 5.0, Draw::colors_RWB );
        //DEBUG
        if(gridFF.Vpot){
            printf( "Vpot min %g max %g \n", VecN::min(gridFF.ntot,gridFF.Vpot), VecN::max(gridFF.ntot,gridFF.Vpot) );
            Draw3D::drawScalarGrid( gridFF.n, {-gridFF.pos0.y,0.0,gridFF.pos0.x}, {0.0,0.0,gridFF.step.x}, {-gridFF.step.y,0.0,0.0}, gridFF.Vpot, -1.0, 1.0, Draw::colors_RWB );
            //Draw3D::drawScalarGrid( gridFF.n, {0,gridFF.pos0.y,gridFF.pos0.x}, {0.0,0.0,gridFF.step.x}, {-gridFF.step.y,0.0,0.0}, logHist,  -5.0, 5.0 );
        }
        //DEBUG

        if(bPlotGridForceField){
        glColor3f(0.,0.,0.);
        Draw3D::drawVectorGrid( gridFF.n, {-gridFF.pos0.y,0.0,gridFF.pos0.x}, {0.0,0.0,gridFF.step.x}, {-gridFF.step.y,0.0,0.0}, gridFF.force, 10.0 );
        }


    }

};


void TestAppElectromagnetic::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_m:  bIntegrate=!bIntegrate;
                    if(bIntegrate==false){ gridFF.solvePoisson(); }
                    break;
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
