
#ifndef basisIntegrals_h
#define basisIntegrals_h

#include <stdint.h>

#include "spline_hermite.h"
#include "integration.h"

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

void intCyl_shift( int nr, int nz, int nint, const double* f1s, const double* f2s, const double* ws, double* Is, double anti1, double anti2 ){
    double sym = anti1*anti2;
    bool bSym = sym*sym>0.125;
    //  to Do - can be accelerated by FFT
    //int ir    = nr-1;
    //double wr = 1.0;
    for(int ir=0; ir<nr; ir++){
        double wr = ws[ir];
        int iroff = ir*nz;
        if (bSym){ for(int ii=0;ii<nint;ii++){ Is[ii] += wr*dot_shifted_sym( nz, ii, f1s+iroff, f2s+iroff, anti1, anti2 ); } }
        else     { for(int ii=0;ii<nint;ii++){ Is[ii] += wr*dot_shifted    ( nz, ii, f1s+iroff, f2s+iroff               ); } }
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

void integrateSP( int nrf, int order, int nint, double dz, double Rmax, double frStep, const double** fr1, const double** fr2, double** Is ){
    constexpr const int nr = 8;
    constexpr const double *ws_ = GaussQuadrature::ws_8;
    constexpr const double *rs  = GaussQuadrature::xs_8;

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

#endif



