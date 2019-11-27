
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
#include "quaternion.h"
#include "Mat3.h"
#include "Mat4.h"
//#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

#include "AOIntegrals.h"
#include "AOrotations.h"

#include "testUtils.h"

//#include "MMFF.h"
//#define R2SAFE  1.0e-8f

//#include "eFF.h"
//#include "e2FF.h"


double hrot_pp( const Vec3d& dh, const Vec3d& p1, const Vec3d& p2, const Vec2d& H ){
    Mat3d rot;
    rot.fromDirUp(dh,p1);
    Vec3d c1,c2;
    rot.dot_to( p1, c1 );
    rot.dot_to( p2, c2 );
    return (c1.x*c2.x + c1.y*c2.y)*H.x + c1.z*c2.z*H.y;
}

double hrot_sp( const Vec3d& dh, const Quat4d& c1, const Quat4d& c2, const Vec3d& H ){

    Mat3d rot;
    rot.fromDirUp(dh,c1.p);
    Vec3d p1,p2;
    rot.dot_to( c1.p, p1 );
    rot.dot_to( c2.p, p2 );

    return (p1.x*p2.x + p1.y*p2.y)*H.y + p1.z*p2.z*H.z + c1.s*c2.s*H.x;

}

inline double slater( Vec3d p, const Quat4d& c, double beta ){
    //double DEBUG_xy = p.y;
    //double DEBUG_z  = p.z;
    double r = p.normalize();
    //printf( "slater xy %g z %g %g \n", DEBUG_xy, DEBUG_z, r );
    double e = exp( -beta*r );
    return e*( c.s + c.p.dot(p) );
}

//template<double (*func)(const Vec3d& p) >
template< typename Func >
double integrateMidpoint3D( Func func, double h, Vec3d pmin, Vec3d pmax ){
    Vec3d dp = pmax-pmin;
    double invh = 1/h;
    Vec3i n = { round(dp.x*invh), round(dp.y*invh), round(dp.z*invh) };
    Vec3d d = { dp.x/n.x, dp.y/n.x, dp.z/n.x };
    pmin.add_mul( d, 0.5 );
    double sum = 0;
    for(int ix=0; ix<n.x; ix++){
        for(int iy=0; iy<n.y; iy++){
            for(int iz=0; iz<n.z; iz++){
                Vec3d p = pmin + (Vec3d){d.x*ix,d.y*iy,d.z*iz};
                sum +=  func(p);
            }
        }
    }
    return sum * d.totprod();
}


//template< double (*func)(const Vec3d& p)>
template< typename Func >
void drawIso( Func func, Vec3i n, Vec3d pmin, Mat3d dcell, double iso ){
    const int np = 4+6+1;
    Vec3f  ps  [np];
    double vals[np];
    const Vec3d p0s[np]{
        0.0,0.0,0.0,    // 0
        1.0,0.0,0.0,    // 1
        0.0,1.0,0.0,    // 2
        0.0,0.0,1.0,    // 3
        +0.5,+0.5,+0.5, // 4
        +0.5,-0.5,-0.5, // 5
        -0.5,+0.5,-0.5, // 6
        -0.5,-0.5,+0.5, // 7
        -0.5,+0.5,+0.5, // 8
        +0.5,-0.5,+0.5, // 9
        +0.5,+0.5,-0.5, // 10
    };
    const int nt = 12;
    const Quat4i tis[nt]{
        0,1, 4,9,
        0,1, 9,5,
        0,1, 5,10,
        0,1, 10,4,

        0,2, 4,10,
        0,2, 10,6,
        0,2, 6,8,
        0,2, 8,4,

        0,3, 4,8,
        0,3, 8,7,
        0,3, 7,9,
        0,3, 9,4
    };
    //glBegin(GL_LINES);
    glBegin(GL_TRIANGLES);
    for(int ix=0; ix<n.x; ix++){
        for(int iy=0; iy<n.y; iy++){
            for(int iz=0; iz<n.z; iz++){
                Vec3d p =  pmin + dcell.a*ix + dcell.b*iy + dcell.c*iz ;
                for(int i=0; i<np; i++){
                    Vec3d pi = p+dcell.dot( p0s[i] );
                    ps   [i] = (Vec3f)pi;
                    vals [i] = func(pi) - iso;
                    //printf( "(%i,%i,%i|%i) %g \n", ix,iy,iz, i, vals[i] );
                }
                Vec3f* pps[4];
                for(int i=0; i<nt; i++){
                    Quat4i ti = tis[i];
                    pps[0]=ps+ti.x; pps[1]=ps+ti.y; pps[2]=ps+ti.z; pps[3]=ps+ti.w;
                    Draw3D::drawTetraIso( pps, (Quat4d){vals[ti.x],vals[ti.y],vals[ti.z],vals[ti.w]} );
                    //if( (i>=8) && (i<12) )
                    /*
                    if( vals[ti.x]*vals[ti.y]*vals[ti.z]*vals[ti.w] < 0 ){
                        glColor3f(0.0,0.0,1.0);
                        Draw3D::drawSimplexLines( pps );
                    }
                    */
                    /*
                    for(int j=0; j<4; j++){
                        int jj = ti.array[j];
                        if( vals[jj]>0 ){ glColor3f(0.0,0.0,1.0); Draw3D::drawPointCross_bare( ps[jj], 0.1 ); }
                        else            { glColor3f(1.0,0.0,0.0); Draw3D::drawPointCross_bare( ps[ti.array[j]], 0.1 ); }
                    }
                    */
                }
            }
        }
    }
    glEnd();
}


template< typename Func >
double drawIso( Func func, double h, Vec3d pmin, Vec3d pmax, double  iso ){
    Vec3d dp = pmax-pmin;
    double invh = 1/h;
    Vec3i n = { round(dp.x*invh), round(dp.y*invh), round(dp.z*invh) };
    Mat3d dcell = Mat3dZero;
    dcell.xx = dp.x/n.x;
    dcell.yy = dp.y/n.x;
    dcell.zz = dp.z/n.x;
    drawIso( func, n, pmin, dcell, iso );
}



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
        //double Iref =0;
        double Iref = integrateMidpoint3D( wf, 0.2, pmin, pmax );
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
    plot1.render();

    p2.z=2.0;
    ogl=Draw::list(ogl);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LIGHTING);
        glShadeModel(GL_SMOOTH);
        glColor3f(0.0,0.4,0.8); drawIso( wf, 0.5, pmin, pmax,  0.01 );
        glColor3f(0.8,0.4,0.0); drawIso( wf, 0.5, pmin, pmax, -0.01 );
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
















