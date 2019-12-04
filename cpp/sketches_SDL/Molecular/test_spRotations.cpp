
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "DrawIso.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Mat3.h"
#include "Mat4.h"
//#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"

#include "Plot2D.h"

#include "integration.h"
#include "AOIntegrals.h"
#include "AOrotations.h"

#include "testUtils.h"

#include "Lingebra.h"
#include "approximation.h"

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"

class TestAppSp3Space: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //RigidAtomType type1,type2;

    bool bRun = false;

    Mat3d rot;         // rotation matrix between atoms
    Vec3d p1,p2;       // position of atoms 1,2
    Quat4d psi1,psi2;  // spxyz orbital expansion on atoms 1,2

    std::vector<Vec3d> ps{ 1,Vec3dZero };


    Mat4d orbs;

    int ogl=0;

    int ipicked  = -1, ibpicked = -1;

    Plot2D plot1;

    //double Emin,Emax;
    //int     npoints;
    //Vec3d*  points  =0;
    //double* Energies=0;
    //Vec3d * Forces  =0;

    int      fontTex;

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppSp3Space( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppSp3Space::TestAppSp3Space( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );

    double beta1=1.0,beta2=2.0;

    p1 = (Vec3d ){-1.00,0.50,-0.25};
    p2 = (Vec3d ){ 1.00,-0.25,0.52};  // position of atoms 1,2
    psi1=(Quat4d){ 4.00,-2.00,-3.00, 0.};
    psi2=(Quat4d){ 4.00,3.10,2.01, 0.};

    psi1.normalize();
    psi2.normalize();

    Vec3d H = (Vec3d){1,1,1};

    Vec3d  h = p1-p2;
    double r = h.normalize();
    rot.fromDirUp(h,psi1.p);

    double h12 = hrot_sp( h, psi1, psi2, H );
    printf( "h12 %g \n", h12 );

    printf( "norms %g %g \n", psi1.norm2(), psi2.norm2() );

    auto wf = [&](const Vec3d& p)->double{
        //double w1 = 1-p.norm();
        double w1 = slater( p-p1, psi1, beta1 );
        double w2 = slater( p-p2, psi2, beta2 );
        //printf( "%g | %g %g \n", w1*w2, w1, w2 );
        return  w1 * w2;
        //return  w2;
    };

    auto wf2D_1 = [&](double rxy, double z)->double{
        //printf( "wf2D_1 xy %g z %g \n", rxy, z );
        return slater( {0,rxy,z}, psi1, beta1 );
    };

    auto wf2D_2 = [&](double rxy, double z)->double{
        return slater( {0,rxy,z}, psi2, beta2 );
    };

    // ===== Check Numerical Integral

    double xmax = 8.0;
    int    n    = 40;
    double dx   = xmax/n;

    //void projectFr( Func func, int nr, int nz, double dz, const double* rs, double* f );

    // ===== Check Numerical Integral
    beta1=1.8; beta2=1.8;
    p1 = (Vec3d ){ 0.0,0.0,0.0};
    p2 = (Vec3d ){ 0.0,0.0,0.0};  // position of atoms 1,2
    //psi1=(Quat4d){ 0.0, 0.0, 0.0, 1.0 };
    //psi2=(Quat4d){ 0.0, 0.0, 0.0, 1.0 };
    //psi1=(Quat4d){ 0.0, 0.0, 1.0, 0.0 };
    //psi2=(Quat4d){ 0.0, 0.0, 1.0, 0.0 };
    psi1=(Quat4d){ 0.0, 1.0, 0.0, 0.0 };
    psi2=(Quat4d){ 0.0, 1.0, 0.0, 0.0 };

    double* IsRef = new double[n];
    double* Is2D  = new double[n];


    plot1.init();
    plot1.clrGrid = 0xFF858585;
    DataLine2D* lref = plot1.add( new DataLine2D(n,0,dx,0xFF000000) );
    //DataLine2D* l2D  = plot1.add( new DataLine2D(n,0,dx,0xFF008000) );

    DataLine2D* lss  = plot1.add( new DataLine2D(n,0,dx,0xFF00FF00) );
    DataLine2D* lsz  = plot1.add( new DataLine2D(n,0,dx,0xFF0080FF) );
    DataLine2D* lzs  = plot1.add( new DataLine2D(n,0,dx,0xFF00FFFF) );
    DataLine2D* lzz  = plot1.add( new DataLine2D(n,0,dx,0xFF0000FF) );
    DataLine2D* lyy  = plot1.add( new DataLine2D(n,0,dx,0xFFFF0000) );

    DataLine2D* lpoly = plot1.add( new DataLine2D(n,0,dx,0xFFFFFFFF) );


    Vec3d pmin={-4,-4,-4};
    Vec3d pmax={4,4,4+xmax};
    double Rmax = pmax.y;

    //integrateCylFunc( wf2D_1, wf2D_2, 0, n, dx, Rmax, Is2D, 1., 1. );

    double dr   = 0.1;
    int    nr   = Rmax/dr + 3;
    double* frs = new double[nr];

    for(int ir=0; ir<nr; ir++){
        double r = ((ir-1)*dr); // WARRNING : don't forget [-1] !!!!
        frs[ir] = exp(-beta1*r);
    }

    double *Iss = new double[n];
    double *Isz = new double[n];
    double *Izs = new double[n];
    double *Izz = new double[n];
    double *Iyy = new double[n];

    const double* fr1[]={ frs, frs, frs };
    const double* fr2[]={ frs, frs, frs };
    double*       Is []={ Iss, Isz, Izs, Izz, Iyy };

    long t1 = getCPUticks();
    integrateSP( nr, 0, n, dx, Rmax, dr, fr1, fr2, Is );
    long t12 = getCPUticks() - t1;
    printf( "time{integrateSP} %g [Mtick] | n(%i,%i,%i) %g tick/op \n", t12*1e-6, 14, n, n, t12/(double)(14*n*n) );

    //double rs[1] = {0.2};
    //projectFr(         1, n, nr, dx, dr, 1.0, rs, frs, Isz, 0 );
    //projectFr( wf2D_1, 1, n,     dx, 1.0, rs,      Isz    );

    //exit(0);

    /*
    int nint    = n;
    double dz   = dx;
    double Rmax = pmax.y;

    constexpr const int nr = 14;
    constexpr const double *ws_ = GaussQuadrature::ws_14;
    constexpr const double *rs  = GaussQuadrature::xs_14;
    double *ws=new double[nr];
    double cw = Rmax*Rmax*(M_PI*2)*dz;
    for(int i=0; i<nr; i++){ ws[i] = ws_[i]*rs[i]*cw; };
    const int nz  = nint+1;
    const int nrz = nr*nz;
    double * f1s  = new double[nrz];
    double * f2s  = new double[nrz];
    projectFr( wf2D_1, nr, nz, dz, Rmax, rs, f1s);
    projectFr( wf2D_2, nr, nz, dz, Rmax, rs, f2s);
    for(int i=0; i<nint; i++){ Is2D[i]=0; }
    intCyl_shift( nr, nz, nint, f1s, f2s, ws, Is2D,  1, 1 );

    DataLine2D* lw1  = plot1.add( new DataLine2D(n,0,dx,0xFFFF4000) );
    DataLine2D* lw2  = plot1.add( new DataLine2D(n,0,dx,0xFF0040FF) );
    int iroff = (nr-1)*nz;
    for(int i=0; i<n; i++){
        lw1->ys[i] = f1s[iroff+i];
        lw2->ys[i] = f2s[iroff+i];
        printf( "lw1[%i] x,y %g,%g \n", i, lw1->xs[i], lw1->ys[i] );
    }
    */

    /*
    {
        int nn=10;
        double fs1[nn];
        double fs2[nn];
        double fh1[nn];
        double fh2[nn];
        for( int i=0; i<nn; i++ ){ fs1[i]=0; fs2[i]=0; fh1[i]=0; fh2[i]=0; }
        int i0=nn/2;
        fs1[i0-1]=0.5; fs1[i0]=1; fs1[i0+1]=0.5;
        fs2[i0-1]=0.5; fs2[i0]=1; fs2[i0+1]=0.5;
        fh1[0]=1; fh1[1]=0.5;
        fh2[0]=1; fh2[1]=0.5;
        for( int ishift=0; ishift<5; ishift++ ){
            double I1 = dot_rolled     ( nn, ishift, fs1, fs2 );
            double I2 = dot_shifted_sym( nn, ishift, fh1, fh2, 1, 1 );
            printf( "ishift %i I1 %g  I2 %g \n", ishift, I1, I2 );
        }
    }
    exit(0);
    */




    for(int i=0; i<n; i++){
        //double x = i*dx;
        double x = lref->xs[i];
        p2.z = x;
        double Iref =0;
        //double Iref = integrateMidpoint3D( wf, 0.2, pmin, pmax );
        IsRef[i] = Iref;
        lref->ys[i] = Iref   ;
        //l2D ->ys[i] = Is2D[i];
        lss ->ys[i] = Iss [i];
        lsz ->ys[i] = Isz [i];
        lzs ->ys[i] = Izs [i];
        lzz ->ys[i] = Izz [i];
        lyy ->ys[i] = Iyy [i];
        //printf( "Iref[%02i] x,f(x):  %g    %g %g   %g %g \n", i, x, Iref, Is2D[i], Iss[i], Isz[i] );
    }

    int npoly = 6;
    double coefs[npoly];
    Approx ::polyFit ( lsz->n,   npoly, lsz->xs,   lsz->ys,   coefs );
    Approx ::polyeval( lpoly->n, npoly, lpoly->xs, lpoly->ys, coefs );
    //polyfit( n, 5, xs, double* ys, double* BB, double* By );
    for(int i=0; i<lpoly->n; i++){ printf( "lpoly[%i] %g -> %g | %g \n", i, lpoly->xs[i], lpoly->ys[i], lsz->ys[i] ); }
    for(int i=0; i<npoly;    i++){ printf( "poly[%i] %g \n", i, coefs[i] ); }


    /*
    int i0 = 10;
    int n_ = n-i0;
    int    npows = 4;
    double ypows[npows]{ 1., 1./2, 1./4, 1./8 };
    double errs [2*npows];
    double coefss[npoly*npows];

    Approx::ypowsApprox( npows, n_, npoly, lsz->xs+i0, lzs->ys+i0, coefss, ypows, errs );
    for(int i=0; i<npows; i++){
        //printf( "[%i] ypow %g err %g rmse %g \n", i, ypows[i], errs[i*2], errs[i*2+1] );
        printf( "[%i] ypow %g err %g rmse %g \n", i, ypows[i], errs[i*2], sqrt(errs[i*2+1]/n_) );
    }
    */


    // ======= test AutoAprox

    Approx::AutoApprox aaprox;

    int npoints = 100;
    int npows   = 4;
    int npolys  = 15;
    //double pows [npows] {1,2,4,8};
    double pows [npows] {-1,-2,-4,-8};
    int    polys[npolys]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    aaprox.bindOrRealloc( npolys, npows, npoints, polys, pows );
    aaprox.ws.alloc(aaprox.np);
    //aaprox.ws = aaprox.ys_ref;
     DEBUG
    for(int i=0; i<npoints; i++){
        double x         = i*0.1;
        aaprox.xs    [i] = x;
        aaprox.ys_ref[i] = exp(-x);
        //aaprox.ws    [i] = exp(-x*(7./8));
        aaprox.ws    [i] = exp(-x*(9./8.));
        //aaprox.ync   [i] = 1+x;
        //aaprox.ws    [i] = exp(-x);
        //aaprox.ys_ref[i] = x*exp(-x);
        //aaprox.ys_ref[i] = pow(x,3);
    }
    DEBUG
    aaprox.preparePowers();
    DEBUG
    for(int i=0; i<aaprox.npows; i++){
        int order = aaprox.tryVariant(10, 50, i );
        if(order<0){ order=aaprox.npoly; printf("(not converged)"); }
        printf("[%i] pow %g err %g coefs[%i]: ", i, aaprox.pows[i], aaprox.err, order );
        for(int j=0; j<order; j++ ) printf( "%g ", i, aaprox.coefs[j] );
        printf("\n");
    }
     DEBUG
    exit(0);

    plot1.render();

    p2.z=2.0;
    ogl=Draw::list(ogl);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LIGHTING);
        glShadeModel(GL_SMOOTH);
        glColor3f(0.0,0.4,0.8); Draw3D::drawIso( wf, 0.5, pmin, pmax,  0.01 );
        glColor3f(0.8,0.4,0.0); Draw3D::drawIso( wf, 0.5, pmin, pmax, -0.01 );
        glColor3f(0,0,0);
        Draw3D::drawBBox( pmin, pmax );

        //glEnd();
    glEndList();

}

void TestAppSp3Space::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);

    glCallList(ogl);

    Vec3d  h = p1-p2;
    double r = h.normalize();
    rot.fromDirUp(h,psi1.p);

    Vec3d c1,c2;
    rot.dot_to( psi1.p, c1 );
    rot.dot_to( psi2.p, c2 );

    glColor3f(0.0,0.0,0.0);

    glColor3f(0.0,0.0,0.0);
    Draw3D::drawPointCross(p1,0.1);
    Draw3D::drawPointCross(p2,0.1);
    //Draw3D::drawLine(p1,p2);

    //Draw3D::drawMatInPos( Mat3dIdentity, p1, psi1.p );
    //Draw3D::drawMatInPos( Mat3dIdentity, p2, psi2.p );
    glColor3f(1.0,1.0,1.0);
    Draw3D::drawVecInPos( psi1.p, p1 );
    Draw3D::drawVecInPos( psi2.p, p2 );

    Draw3D::drawMatInPos( rot, p1, c1 );
    Draw3D::drawMatInPos( rot, p2, c2 );

    //Draw3D::drawPoints(ps.size(),&ps[0],0.1);

    //Draw3D::drawAxis(1.5);


};


void TestAppSp3Space::drawHUD(){
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
	glTranslatef( 100.0,100.0,0.0 );
	glScalef    ( 20.0,300.00,1.0  );
	plot1.view();

}


void TestAppSp3Space::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppSp3Space* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppSp3Space( junk , 800, 600 );
	thisApp->loop( 1000000 );
	return 0;
}
















