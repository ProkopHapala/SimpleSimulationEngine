
/// @file @brief  This is a computational fluid dynamics (CFD) demo for aerodynamics. It uses the vortex lattice method (`PotentialFlow.h`) to calculate the airflow and lift generated by a wing. The wing is represented as a "lift-line," and the demo visualizes the resulting velocity field around it. It also includes functions for numerical integration to verify the analytical results.
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

#include "SDL_utils.h"
#include "Plot2D.h"

#include "PotentialFlow.h"
//#include "grids3D.h"
//#include "Multipoles.h"
//#include "MultipoleGrid.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"

#include "DrawField.h"

#define R2SAFE  1.0e-8f

// ======= THE CLASS

typedef Vec3d (*VecFieldFunc)(Vec3d R);

Vec3d vInf      = {0.0,0.0,-1.0};
Vec3d flightDir = vInf*(-1.0/vInf.norm());

const int nCPs  = 3;
Vec3d CPs[nCPs] = { {0.0,0.0,0.0}, {1.0,0.0,0.0}, {2.0,0.0,-0.5} };
double Gamas[nCPs-1] = { 1.0, 1.0 };

Vec3d totVelField( Vec3d R ){
    Vec3d B = (Vec3d){0.0,0.0,0.0};
    //Vec3d B = vInf;
    B.add( vInf );
    for(int i=0; i<nCPs-1; i++){
        //B.add( horseshoe( R, CPs[i], CPs[i+1], flightDir, Gamas[nCPs] ) );
        //horseshoe( B, R, CPs[i], CPs[i+1], flightDir, Gamas[i] );
        horseshoeDecay( B, R, CPs[i], CPs[i+1], flightDir, Gamas[i], 0.01 );
        //B.add( sourceDipol( R, {0.0,0.0,0.0,1.0} ) );
        //B.add( sourceDipol( R+(Vec3d){0.0,0.2,0.0}, {0.0,0.0,1.0,0.0} ) );
        //B.add( pointSource( R ) * 1.0 );
    }
    //exit(0);
    //printf("=====\n");
    return B;

    //return horseshoe( R, {1.0,0.0,0.0}, {-1.0,0.0,0.0}, {0.0,0.0,1.0}, 1.0 );
};




/*
template<typename Type_f, typename Type_F>
void TestIntegralScalar1D( int n, Vec3d p, Vec3d dp, Type_f f, Type_F f){

    for(){

    }
}
*/


template<typename Type_f>
double numIntegralScalarLine1D( int n, Vec3d p0, Vec3d p1, Type_f f, double* ys=0 ){
    Vec3d dp = p1-p0;
    dp.mul(1./n);
    double dl = dp.norm();
    //double dl = 1;
    Vec3d p  = p0 + dp*0.5;
    double sum=0;
    for(int i=0; i<n; i++){
        double y =  f( p ) * dl;
        sum  += y;
        if(ys)ys[i] = sum;
        p.add(dp);
    }
    return sum;
}



template<typename Type_f>
Vec3d numIntegralVec3Line1D( int n, Vec3d p0, Vec3d p1, Type_f f, Vec3d* ys=0 ){
    Vec3d dp = p1-p0;
    dp.mul(1./n);
    double dl = dp.norm();
    //double dl = 1;
    Vec3d p  = p0 + dp*0.5;
    Vec3d sum=Vec3dZero;
    for(int i=0; i<n; i++){
        Vec3d y =  f( p ) * dl;
        sum.add(y);
        if(ys)ys[i] = sum;
        p.add(dp);
    }
    return sum;
}

template<typename Type_F>
double definiteIntegral2Array( int n, Vec3d p0, Vec3d p1, Type_F F, double * ys){
    Vec3d dp = p1-p0;
    dp.mul(1./n);
    Vec3d p  = p0;
    double Y0 = F( p0 );
    double y;
    for(int i=0; i<n; i++){
        p.add(dp);
        y     = F( p ) - Y0;
        ys[i] = y;
    }
    return y;
}

template<typename Type_F>
Vec3d definiteIntegral2ArrayVec3( int n, Vec3d p0, Vec3d p1, Type_F F, Vec3d * ys){
    Vec3d dp = p1-p0;
    dp.mul(1./n);
    Vec3d p  = p0;
    Vec3d Y0 = F( p0 );
    Vec3d y;
    for(int i=0; i<n; i++){
        p.add(dp);
        y     = F( p ) - Y0;
        ys[i] = y;
    }
    return y;
}

template<typename Type_f, typename Type_F>
double testIntegralScalarLine1D( int n, Vec3d p0, Vec3d p1, Type_f f, Type_F F){
    double y_num = numIntegralScalarLine1D( n, p0, p1, f);
    double y_ana = F( p0, p1 );
    double err   = y_ana - y_num;
    printf( "Integral error %g | analytic %g numerical %g \n", y_ana, y_num, err );
    return err;
}





void testILineFinite(){
    // test ILineFinite
    //Vec3d L  = (Vec3d){4.0,0.5,0.6};
    //Vec3d R0 = (Vec3d){-3.0,1.0,2.0};

    Vec3d L  = (Vec3d){0.0,4.0,0.0};
    Vec3d R0 = (Vec3d){1.0,1.0,0.0};

    double l = L.norm();
    Vec3d hL = L*(1.0/l);
    //Vec3d B_ana = ILineFinite( R0, hL, l );
    Vec3d B_ana = ISemiInfSheet( R0, {1.0,0.0,0.0}, hL, l );

    int nsamp=1000;
    Vec3d B_num = (Vec3d){0.0,0.0,0.0};
    double step = 1.0/nsamp;
    Vec3d dL    = L*step;
    Vec3d R     = R0 + dL*0.5;
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0,0.0,1.0);
    for(int i=0; i<nsamp; i++){
        //B_num.add( dBiotSawart( R, dL ) );
        B_num.add_mul( ILineSemiInf( R, {1.0,0.0,0.0} ), step*l );
        //printf( "%f %f %f\n", R.y, ILineSemiInf( R, {1.0,0.0,0.0} ).norm()/VortexQuantum, B_num.norm()/VortexQuantum );
        R.add(dL);
        glVertex3f( l*step*i, B_num.z*20.0, 0 );
    }
    glEnd();

    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,0.0,0.0);
    for(int i=0; i<nsamp; i++){
        //B_ana = ILineFinite( R0, hL, l*step*i );
        B_ana = ISemiInfSheet( R0, {1.0,0.0,0.0}, hL, l*step*i );
        glVertex3f( l*step*i, B_ana.z*20.0, 0 );
    }
    glEnd();
    //exit(0);
    glColor3f(0.0,0.0,0.0);
    Draw3D::drawLine( (Vec3f){0.0,0.0,0.0}, (Vec3f){10.0,0.0,0.0} );

    //printf(" B_ana (%f,%f,%f) B_num (%f,%f,%f) \n", B_ana.x,B_ana.y,B_ana.z, B_num.x,B_num.y,B_num.z );
    printf(" |B_ana| %f |B_num| %f \n", B_ana.norm()/VortexQuantum, B_num.norm()/VortexQuantum );
    //exit(0);
}

class TestAppVortexLattice : public AppSDL2OGL_3D {
	public:

    //MultipoleGrid grid;

    Plot2D plot1;
    //int fontTex;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppVortexLattice( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppVortexLattice::TestAppVortexLattice( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    /*
    //exit(0);
    double x0 = 6;
    double y0 = 10;
    testIntegralScalarLine1D( 10000, {-0.1,0,0}, {0.6,0,0},
        [x0,y0](Vec3d p          )->double{ double dx  = p.x-x0;                        return     dx/(dx*dx + y0*y0);                                }, // f
        [x0,y0](Vec3d p0,Vec3d p1)->double{ double dx0 = p0.x-x0; double dx1 = p1.x-x0; return 0.5*log(dx1*dx1 + y0*y0 ) - 0.5*log(dx0*dx0 + y0*y0 ); }  // F
    );
    exit(0);
    */

    plot1.init();
    //plot1.fontTex = fontTex;
    plot1.bGrid=false;
    //plot1.bAxes=false;
    plot1.bTicks=false;
    plot1.scaling.y = (5.0);

    int nsamp = 100;
    double x0 = 0.2;
    double y0 = 0.3;
    Vec3d p0 = (Vec3d){-0.3,   0,0};
    Vec3d p1 = (Vec3d){+0.6,   0,0};
    double dx = (p1.x-p0.x)/nsamp;
    DataLine2D* line_FnunX = new DataLine2D(nsamp,p0.x,dx*5, 0xFF0000 ); plot1.add(line_FnunX);
    DataLine2D* line_FnunY = new DataLine2D(nsamp,p0.x,dx*5, 0x0000FF ); plot1.add(line_FnunY);
    DataLine2D* line_FanaX = new DataLine2D(nsamp,p0.x,dx*5, 0xFFFF00 ); plot1.add(line_FanaX); line_FanaX->lineStyle=' '; line_FanaX->pointStyle='.'; //line_Fana->pointSize=0.0003;
    DataLine2D* line_FanaY = new DataLine2D(nsamp,p0.x,dx*5, 0x00FFFF ); plot1.add(line_FanaY); line_FanaY->lineStyle=' '; line_FanaY->pointStyle='.'; //line_Fana->pointSize=0.0003;


    // #### Coulomb Law   vs   Biot-Savart Law
    //   Coulomb       E = kQ* r/|r|^3     V=kQ/|r|                 Lapalce{ V } = rho   (charge  density)
    //   Biot-Savart   B = J x r/|r|^3     A= J/|r|                 Laplace{ A } = J     (current density)
    //                            (Ax,Ay,Az)=(Jx,Jy,Jz)/|r|


//    // ## source line: rho(x) = 1
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return Fx; }, line_FnunX->ys );
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return Fy; }, line_FnunY->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow0(x,y0,Fx,Fy); return Fx; }, line_FanaX->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow0(x,y0,Fx,Fy); return Fy; }, line_FanaY->ys );

//    // ##  source line rho(x) = x
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return x*Fx; }, line_FnunX->ys );
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return x*Fy; }, line_FnunY->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow1(x,y0,Fx,Fy); return Fx;    }, line_FanaX->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow1(x,y0,Fx,Fy); return Fy;    }, line_FanaY->ys );

//    // ##  source line rho(x) = x^2
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return x*x*Fx; }, line_FnunX->ys );
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return x*x*Fy; }, line_FnunY->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow2(x,y0,Fx,Fy); return Fx;     }, line_FanaX->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow2(x,y0,Fx,Fy); return Fy;     }, line_FanaY->ys );

    // ##  source line rho(x) = x^3
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return x*x*x*Fx; }, line_FnunX->ys );
//    numIntegralScalarLine1D( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY        (x,y0,Fx,Fy); return x*x*x*Fy; }, line_FnunY->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow3(x,y0,Fx,Fy); return Fx;       }, line_FanaX->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [x0,y0](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_pow3(x,y0,Fx,Fy); return Fy;       }, line_FanaY->ys );

    double C0=0.5,C1=-0.7,C2=0.8,C3=-5.6;

    //    // ##  source line rho(x) = C0 + C1*x
//    numIntegralScalarLine1D( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY          (x,y0,Fx,Fy);       return (C0+C1*x)*Fx; }, line_FnunX->ys );
//    numIntegralScalarLine1D( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY          (x,y0,Fx,Fy);       return (C0+C1*x)*Fy; }, line_FnunY->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_linear(x,y0,Fx,Fy,C0,C1); return Fx;           }, line_FanaX->ys );
//    definiteIntegral2Array ( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_linear(x,y0,Fx,Fy,C0,C1); return Fy;           }, line_FanaY->ys );

    numIntegralScalarLine1D( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY          (x,y0,Fx,Fy);            return (C0+x*(C1+x*(C2+C3*x)))*Fx; }, line_FnunX->ys );
    numIntegralScalarLine1D( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceXY          (x,y0,Fx,Fy);            return (C0+x*(C1+x*(C2+C3*x)))*Fy; }, line_FnunY->ys );
    definiteIntegral2Array ( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_poly3(x,y0,Fx,Fy,C0,C1,C2,C3); return Fx;                         }, line_FanaX->ys );
    definiteIntegral2Array ( nsamp , p0, p1, [=](Vec3d p){ double x = p.x-x0; double Fx,Fy; sourceLineF_poly3(x,y0,Fx,Fy,C0,C1,C2,C3); return Fy;                         }, line_FanaY->ys );

    plot1.render();

}

void TestAppVortexLattice::draw(){
    printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    plot1.view();
    //return;


	//LineSemiInf( {0.5,1.5,0.0}, {1.0,0.0,0.0} ); exit(0);
	//ISemiInfSheet( {1.0,1.0,0.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, 4.0 ); exit(0);

	//testILineFinite(); //exit(0);
    //if(frameCount==0)testILineFinite();
    //glColor3f(1.0,1.0,1.0);
    //Draw3D::drawLine    ( (Vec3d){1.0,0.0, 0.0},  {-1.0,0.0,0.0} );
    //Draw3D::drawVecInPos( (Vec3d){0.0,0.0,10.0},  { 1.0,0.0,0.0} );
    //Draw3D::drawVecInPos( (Vec3d){0.0,0.0,10.0},  {-1.0,0.0,0.0} );


    plotVortexFilaments( nCPs, CPs, flightDir*-100.0 );
    glColor3f(0.0,0.0,0.0);
    //plotVecPlane( {21,21}, { -1.0,0.0,-2.0 }, {0.2,0.0,0.0}, {0.0,0.0,0.2}, 0.02, 1.0, totVelField );
    //plotVecPlane( {81,11}, { -1.0,0.0,-2.0 }, {0.05,0.0,0.0}, {0.0,0.0,0.2}, 0.02, 1.0, totVelField );
    //plotVecPlane( {5,5}, { -2.0,0.0,-2.0 }, {1.0,0.0,0.0}, {0.0,0.0,1.0}, 0.02, 10.0, totVelField );
    //exit(0);
    plotStreamLinePlane( {20,1}, 200, {-1.1,0.0,3.0}, {0.2,0.0,0.0}, {0.0,0.2,0.0}, 0.1, totVelField );


};


void TestAppVortexLattice::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppVortexLattice::drawHUD(){
    glDisable ( GL_LIGHTING );
}

// ===================== MAIN

TestAppVortexLattice * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppVortexLattice( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
