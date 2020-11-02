
#ifndef AOIntegrals_h
#define AOIntegrals_h

/// @file
/// \ingroup Molecular
/// @{

#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "spline_hermite.h"
#include "integration.h"

/// rotate p-p interaction
double hrot_pp( const Vec3d& dh, const Vec3d& p1, const Vec3d& p2, const Vec2d& H ){
    Mat3d rot;
    rot.fromDirUp(dh,p1);
    Vec3d c1,c2;
    rot.dot_to( p1, c1 );
    rot.dot_to( p2, c2 );
    return (c1.x*c2.x + c1.y*c2.y)*H.x + c1.z*c2.z*H.y;
}

/// rotate s-p interaction
double hrot_sp( const Vec3d& dh, const Quat4d& c1, const Quat4d& c2, const Vec3d& H ){

    Mat3d rot;
    rot.fromDirUp(dh,c1.p);
    Vec3d p1,p2;
    rot.dot_to( c1.p, p1 );
    rot.dot_to( c2.p, p2 );

    return (p1.x*p2.x + p1.y*p2.y)*H.y + p1.z*p2.z*H.z + c1.s*c2.s*H.x;

}

/// evaluate slater orbital
inline double slater( Vec3d p, const Quat4d& c, double beta ){
    //double DEBUG_xy = p.y;
    //double DEBUG_z  = p.z;
    double r = p.normalize();
    //printf( "slater xy %g z %g %g \n", DEBUG_xy, DEBUG_z, r );
    double e = exp( -beta*r );
    return e*( c.s + c.p.dot(p) );
}

double dot_shifted( int n, int ioff, const double* f1, const double* f2 ){
    double sum = 0;
    int n_ = n-ioff;
    for(int i=0;i<n_;i++){ sum+=f1[i+ioff]*f2[i]; }
    return sum;
}

double dot_rolled( int n, int ioff, const double* f1, const double* f2 ){
    // very often integrated function are either symmetric or antisymmetric, than it is enough to use just half of interval
    //printf( "dot_rolled %i %i (%g,%g) \n", n, ioff );
    const int n1 = n-ioff;
    double sum1 = 0; for(int i=0;i<n1  ;i++){ /*printf("sum1 %i %i \n", i+ioff, i    );*/ sum1+=f1[i+ioff]*f2[i   ]; }
    double sum2 = 0; for(int i=0;i<ioff;i++){ /*printf("sum2 %i %i \n", i     , n1+i );*/ sum2+=f1[i     ]*f2[n1+i]; }
    //printf( "dot_shifted_sym %i %i (%g,%g) %g(%g,%g,%g)\n", n, ioff, anti1, anti2, sum1+sum2+sum3, sum1,sum2,sum3 );
    return sum1+sum2;
}

double dot_shifted_sym( int n, int ioff, const double* f1, const double* f2, const double anti1, double anti2 ){
    // very often integrated function are either symmetric or antisymmetric, than it is enough to use just half of interval
    //printf( "dot_shifted_sym %i %i (%g,%g) \n", n, ioff, anti1, anti2 );
    const int n_ = n-ioff;
    //double sum1 = 0; for(int i=0;i<n_  ;i++){ /*printf("sum1 %i %i \n", i+ioff, i );*/ sum1+=f1[i+ioff]*f2[i     ]; }
    //double sum2 = 0; for(int i=0;i<ioff;i++){ /*printf("sum2 %i %i \n", i, ioff-i );*/ sum2+=f1[i     ]*f2[ioff-i]; } sum2*= anti2;
    //double sum3 = 0; for(int i=0;i<n_  ;i++){ /*printf("sum3 %i %i \n", i, ioff+i );*/ sum3+=f1[i     ]*f2[ioff+i]; } sum3*=(anti1*anti2);
    double sum1 = 0; for(int i=0;i<n_  ;i++){ /*printf("sum1 %i %i \n", i+ioff, i );*/ sum1+=f1[i+ioff]*f2[i     ]; }
    double sum2 = 0; for(int i=0;i<ioff;i++){ /*printf("sum2 %i %i \n", i, ioff-i );*/ sum2+=f1[i     ]*f2[ioff-i]; } sum2*= anti2;
    double sum3 = 0; for(int i=1;i<n_  ;i++){ /*printf("sum3 %i %i \n", i, ioff+i );*/ sum3+=f1[i     ]*f2[ioff+i]; } sum3*=(anti1*anti2);
    //printf( "dot_shifted_sym %i %i (%g,%g) %g(%g,%g,%g)\n", n, ioff, anti1, anti2, sum1+sum2+sum3, sum1,sum2,sum3 );
    return sum1+sum2+sum3;
}

double dot_shifted_sym_( int n, int ioff, const double* f1, const double* f2, const double anti1, double anti2 ){
    // very often integrated function are either symmetric or antisymmetric, than it is enough to use just half of interval
    //printf( "dot_shifted_sym %i %i (%g,%g) \n", n, ioff, anti1, anti2 );
    int DEBUG_n = 0;
    const int n_ = n-ioff;
    if(n_<=0){
        //const n__=n_+n;
        double sum2 = 0;
        for(int i=ioff-n+1;i<n;i++){
            int j=ioff-i-1;
            //printf( "[ %i * %i |%i] %g * %g \n", i, j, n, f1[i],f2[j] );
            sum2+=f1[i]*f2[j];
            //DEBUG_n++;
        }; sum2*= anti2;
        //printf(  "DEBUG dot_shifted_sym_ DEBUG_n %i | %i | n=%i \n", DEBUG_n, n*2-ioff-1, n );
        return sum2;
    }
    double sum1 = 0; for(int i=0;i<n_  ;i++){ sum1+=f1[i+ioff]*f2[i       ]; }
    double sum2 = 0; for(int i=0;i<ioff;i++){ sum2+=f1[i     ]*f2[ioff-i-1]; } sum2*= anti2;
    double sum3 = 0; for(int i=1;i<n_  ;i++){ sum3+=f1[i     ]*f2[ioff+i  ]; } sum3*=(anti1*anti2);
    int i,j;
    //double sum1 = 0; for(int ii=0;ii<n_  ;ii++){ DEBUG_n++;i=ii+ioff;j=ii      ; printf( "[ %i * %i ] %g * %g \n", i,j, n, f1[i],f2[j] ); sum1+=f1[i]*f2[j]; }
    //double sum2 = 0; for(int ii=0;ii<ioff;ii++){ DEBUG_n++;i=ii     ;j=ioff-i-1; printf( "[ %i * %i ] %g * %g \n", i,j, n, f1[i],f2[j] ); sum2+=f1[i]*f2[j]; } sum2*= anti2;
    //double sum3 = 0; for(int ii=1;ii<n_  ;ii++){ DEBUG_n++;i=ii     ;j=ioff+i  ; printf( "[ %i * %i ] %g * %g \n", i,j, n, f1[i],f2[j] ); sum3+=f1[i]*f2[j]; } sum3*=(anti1*anti2);
    //printf( "dot_shifted_sym %i %i (%g,%g) %g(%g,%g,%g)\n", n, ioff, anti1, anti2, sum1+sum2+sum3, sum1,sum2,sum3 );
    //printf(  "DEBUG dot_shifted_sym_ DEBUG_n %i | %i | n=%i \n", DEBUG_n, n*2-ioff-1, n );
    return sum1+sum2+sum3;
}

void intCyl_shift( int nr, int nz, int nint, const double* f1s, const double* f2s, const double* ws, double* Is, double anti1, double anti2 ){
    double sym = anti1*anti2;
    bool bSym = sym*sym>0.125;
    //  to Do - can be accelerated by FFT
    //int ir    = nr-1;
    //double wr = 1.0;
    for(int ir=0; ir<nr; ir++){
        double wr = ws[ir];
        int iroff = ir*nz;
        //if (bSym){ for(int ii=0;ii<nint;ii++){ Is[ii] += wr*dot_shifted_sym( nz, ii, f1s+iroff, f2s+iroff, anti1, anti2 ); } }
        if (bSym){ for(int ii=0;ii<nint;ii++){ Is[ii] += wr*dot_shifted_sym_( nz, ii, f1s+iroff, f2s+iroff, anti1, anti2 ); } }
        else     { for(int ii=0;ii<nint;ii++){ Is[ii] += wr*dot_shifted     ( nz, ii, f1s+iroff, f2s+iroff               ); } }
    }
}

void projectFr( int nr, int nz, int nrf, double dz, double frStep, double rs_scale, const double* rs, const double* fr, double* f, int im ){
    //printf( "===== projectFr ========== \n" );
    double invdr       = 1/frStep;
    double const rsafe = 1e-12;
    int i=0;
    double r2max = (nrf-3)*frStep; r2max*=r2max;
    for(int ir=0; ir<nr; ir++){
        double rxy = rs[ir]*rs_scale;
        for(int iz =0; iz<nz; iz++){
            double  z  = iz*dz;
            double Y;
            double  r2 = rxy*rxy + z*z;
            if(r2>r2max){
                f[i] = 0;
            }else{
                double  r = sqrt( r2 );
                if      (im==0){ Y=1;                       } // s
                else if (im==1){ Y=z            /(r+rsafe); } // pz
                else           { Y=M_SQRT1_2*rxy/(r+rsafe); } // py   // M_SQRT1_2 comes from angular integral:  sqrt(2) = sqrt( Integral{ cos(phi)^2 } )
                // interpolate
                //double  s  = r*invdr;
                //int    is  = (int)s+1;
                //double dr  = s - is;
                //f[i] = Spline_Hermite::val( dr, fr+is ) * Y;
                f[i] = Spline_Hermite::value( r*invdr, fr ) * Y;
                //printf( "projectFr[%i,%i] r %g f %g |Y %g \n", ir,iz,  r,  f[i], Y );
                //printf( "projectFr[%i,%i] xy %g z %g r %g \n", ir,iz,  rxy, z, r );
                //printf( "projectFr[%i,%i] %g,%g,%g   %g | %g %i \n", ir,iz,  rxy, z, r, f[i],Y,   s, is );
            }
            i++;
        }
    }
}

void projectFrLaplace( int nr, int nz, int nrf, double dz, double frStep, double rs_scale, const double* rs, const double* fr, double* f, double* Lf, int im ){
    //printf( "===== projectFr ========== \n" );
    double invdr       = 1/frStep;
    double const rsafe = 1e-12;
    int i=0;
    double r2max = (nrf-3)*frStep; r2max*=r2max;

    printf( "DEBUG projectFrLaplace: nz,nr,nrf %i %i %i \n", nz, nr, nrf );
    printf( "DEBUG projectFrLaplace: Rmax %g  frStep %g rs_scale %g \n", sqrt(r2max), frStep, rs_scale );
    double zmax = dz*nz;
    double V    = M_PI*rs_scale*rs_scale*zmax;
    printf( "DEBUG Cylinder R %g Z %g Volume projectFrLaplace: nz,nr,nrf %i %i %i \n", rs_scale, zmax, V );

    for(int ir=0; ir<nr; ir++){
        double rxy = rs[ir]*rs_scale;
        for(int iz =0; iz<nz; iz++){
            double  z  = iz*dz;
            double Y;
            double  r2 = rxy*rxy + z*z;
            if(r2>r2max){
                f[i] = 0;
            }else{
                double  r    = sqrt( r2 );
                double invr  = 1/(r+rsafe);
                if      (im==0){ Y=1;                  } // s
                else if (im==1){ Y=z            *invr; } // pz
                else           { Y=M_SQRT1_2*rxy*invr; } // py   // M_SQRT1_2 comes from angular integral:  sqrt(2) = sqrt( Integral{ cos(phi)^2 } )
                // https://en.wikipedia.org/wiki/Laplace%27s_equation#Forms_in_different_coordinate_systems
                // L{f} =  1/r^2 d(r^2 * d f ) = 1/r^2 ( 2*r*df + r^2*ddf ) = 2*df/r + ddf
                double y,dy,ddy;
                Spline_Hermite::valdd( r*invdr, fr, y, dy, ddy );
                if(f) f [i] = y*Y;                  // function         ( x,y )
                if(Lf)Lf[i] = (2*dy*invr + ddy)*Y;  // Laplace{function}( x,y )
                //f[i] = Spline_Hermite::value( r*invdr, fr ) * Y;
                //f[i] = exp(-r);
                //f[i] = 1.0/sqrt( V ); // DEBUG
            }
            i++;
        }
    }
}

template< typename Func >
void projectFr( Func func, int nr, int nz, double dz, double Rsc, const double* rs, double* f ){
    int i=0;
    for(int ir=0; ir<nr; ir++){
        double rxy = rs[ir]*Rsc;
        for(int iz =0; iz<nz; iz++){
            double  z  = iz*dz;
            //printf( "<>projectFr[%i,%i] xy %g z %g \n", ir,iz,  rxy, z );
            f[i]       = func( rxy, z );
            //printf( "<>projectFr[%i,%i] r %g f %g \n", ir,iz,  sqrt(rxy*rxy+z*z), f[i] );
            //if(ir==nr-1) printf( "projectFr[%i,%i] r %g z %g f %g \n", ir,iz, rxy, z, f[i] );
            i++;
        }
    }
}

/*
void mulYz( int nr, int nz, double dz, double* rs, double* fin, double* fout ){
    for(int ir=0; ir<nr; ir++){
        double r  = rs[ir];
        for(int iz=0; iz<nz; iz++){
            double  z  = iz*dz;
            double  rf = sqrt( r*r + z*z );

        }
    }
}
*/


void integrateS( int nrf, int order, int nint, double dz, double Rmax, double frStep, const double* fr1, const double* fr2, double* Is ){
    //constexpr const int nr = 8;
    //constexpr const double *ws_ = GaussQuadrature::ws_8;
    //constexpr const double *rs  = GaussQuadrature::xs_8;
    constexpr const int nr = 14;
    constexpr const double *ws_ = GaussQuadrature::ws_14;
    constexpr const double *rs  = GaussQuadrature::xs_14;
    double cw = Rmax*Rmax*(M_PI*2)*dz;
    double *ws=new double[nr];
    for(int i=0; i<nr; i++){ ws[i] = ws_[i]*rs[i]*cw; };
    const int nz  = nint+1;
    const int nrz = nr*nz;
    double * f1  = new double[nrz];
    double * f2  = new double[nrz];
    for(int i=0; i<nint; i++){ Is[i]=0; }
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr1, f1, 0 );
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr2, f2, 0 );
    intCyl_shift( nr, nz, nint, f1, f2, ws, Is,  1,  1 );
    delete [] f1; delete [] f2;
    delete [] ws;
}



//#include "Draw2D.h"

void integrateSK( int nrf, int order, int nint, double dz, double Rmax, double frStep, const double* fr1, const double* fr2, double* ISs, double* IKs ){


    //constexpr const int nr = 8;
    //constexpr const double *ws_ = GaussQuadrature::ws_8;
    //constexpr const double *rs  = GaussQuadrature::xs_8;

    //constexpr const int nr = 14;
    //constexpr const double *ws_ = GaussQuadrature::ws_14;
    //constexpr const double *rs  = GaussQuadrature::xs_14;
    //double cw = Rmax*Rmax*(M_PI*2)*dz;

    // midpoints
    const int nr       = Rmax/frStep;
    const double invnr = 1./nr;
    double ws_[nr];
    double rs [nr];
    double cw = Rmax*Rmax*(M_PI*2)*dz;
    for(int i=0; i<nr; i++){
        ws_[i] = invnr;
        rs [i] = (i+0.5)*invnr;
    }

    double *ws=new double[nr];
    for(int i=0; i<nr; i++){ ws[i] = ws_[i]*rs[i]*cw; }

    const int nz  = nint+1;
    const int nrz = nr*nz;
    double * f1   = new double[nrz];
    double * f2   = new double[nrz];
    //double * Lf1  = new double[nrz];
    double * Lf2  = new double[nrz];
    for(int i=0; i<nint; i++){ ISs[i]=0; IKs[i]=0; }
    projectFrLaplace( nr, nz, nrf, dz, frStep, Rmax, rs, fr1, f1, 0  , 0 );
    projectFrLaplace( nr, nz, nrf, dz, frStep, Rmax, rs, fr2, f2, Lf2, 0 );

    // { // DEBUG
    //     double dr= Rmax/nr;
    //     int ii=0;
    //     for(int ir=0; ir<nr; ir++){
    //         for(int iz=0; iz<nz; iz++){
    //             double x=iz*dz;
    //             double y=ir*dr;
    //             ii=ir*nz+iz;
    //             float f = f1[ii];
    //             //float c = 1/(1+exp(f));
    //             //glColor3f(f,f*5,f*25);
    //             //glColor3f(f,f*0.1,f*10);
    //             glColor3f(tanh(f*5),tanh(f),tanh(f*25));
    //             Draw2D::drawRectangle_d( {x-dz*0.5,y-dr*0.5},{x+dz*0.5,y+dr*0.5} , true);
    //         }
    //     }
    //     //double V   = (2*nz*dz)*Rmax*Rmax*M_PI;
    //     //double dwf = sqrt( 1/V );
    //     //for(int i=0; i<nrz; i++){
    //     //    f1[i]=dwf;
    //     //    f2[i]=dwf;
    //     //}
    // }

    intCyl_shift( nr, nz, nint, f1,  f2, ws, ISs,  1,  1 );
    intCyl_shift( nr, nz, nint, f1, Lf2, ws, IKs,  1,  1 );
    delete [] f1;  delete [] f2; delete [] Lf2; //delete [] Lf1;
    delete [] ws;
}


void integrateSP( int nrf, int order, int nint, double dz, double Rmax, double frStep, const double** fr1, const double** fr2, double** Is ){
    //constexpr const int nr = 8;
    //constexpr const double *ws_ = GaussQuadrature::ws_8;
    //constexpr const double *rs  = GaussQuadrature::xs_8;
    constexpr const int nr = 14;
    constexpr const double *ws_ = GaussQuadrature::ws_14;
    constexpr const double *rs  = GaussQuadrature::xs_14;

    double cw = Rmax*Rmax*(M_PI*2)*dz;
    double *ws=new double[nr];
    for(int i=0; i<nr; i++){ ws[i] = ws_[i]*rs[i]*cw; };

    const int nz  = nint+1;
    const int nrz = nr*nz;

    double * f1s  = new double[nrz];
    double * f2s  = new double[nrz];

    double * f1z  = new double[nrz];
    double * f2z  = new double[nrz];

    const double *fr1s=fr1[0], *fr1z=fr1[1], *fr1y=fr1[2],
                 *fr2s=fr2[0], *fr2z=fr2[1], *fr2y=fr2[2];

    double *Iss=Is[0], *Isz=Is[1], *Izs=Is[2], *Izz=Is[3], *Iyy=Is[4];

    for(int i=0; i<nint; i++){ Iss[i]=0; Isz[i]=0; Izs[i]=0; Izz[i]=0; Iyy[i]=0; }

    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr1s, f1s, 0 );
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr2s, f2s, 0 );
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr1z, f1z, 1 );
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr2z, f2z, 1 );

    intCyl_shift( nr, nz, nint, f1s, f2s, ws, Iss,  1,  1 );
    intCyl_shift( nr, nz, nint, f1s, f2z, ws, Isz,  1, -1 );
    intCyl_shift( nr, nz, nint, f1z, f2s, ws, Izs, -1,  1 );
    intCyl_shift( nr, nz, nint, f1z, f2z, ws, Izz, -1, -1 );

    // we dont need s-function anymore => reuse
    double * f1y = f1s;
    double * f2y = f2s;
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr1y, f1y, 2 );
    projectFr( nr, nz, nrf, dz, frStep, Rmax, rs, fr2y, f2y, 2 );
    intCyl_shift( nr, nz, nint, f1y, f2y, ws, Iyy,  1,  1 );

    delete [] f1s; delete [] f2s;
    delete [] f1z; delete [] f2z;
    //delete [] f1y; delete [] f2y;
    delete [] ws;
}

template<typename Func1,typename Func2>
void integrateCylFunc( Func1 func1, Func2 func2, int order, int nint, double dz, double Rmax, double* Is, double anti1, double anti2 ){
    constexpr const int nr = 8;
    constexpr const double *ws_ = GaussQuadrature::ws_8;
    constexpr const double *rs  = GaussQuadrature::xs_8;
    double *ws=new double[nr];
    double cw = Rmax*Rmax*(M_PI*2)*dz;
    for(int i=0; i<nr; i++){ ws[i] = ws_[i]*rs[i]*cw; };
    const int nz  = nint+1;
    const int nrz = nr*nz;
    double * f1s  = new double[nrz];
    double * f2s  = new double[nrz];
    projectFr( func1, nr, nz, dz, Rmax, rs, f1s);
    projectFr( func2, nr, nz, dz, Rmax, rs, f2s);
    for(int i=0; i<nint; i++){ Is[i]=0; }
    intCyl_shift( nr, nz, nint, f1s, f2s, ws, Is,  anti1, anti2 );
    delete [] f1s; delete [] f2s;
    delete [] ws;
}

///  @}

#endif



