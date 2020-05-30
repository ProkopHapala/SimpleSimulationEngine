
#ifndef GaussianBasis_h
#define GaussianBasis_h

#include "fastmath.h"

namespace Gauss{

/*
######  Product of Gaussians 1D

gij(x) =  exp( wi*(x-xi)^2 ) * exp( wj*(x-xj)^2 )

wi*( x^2 + 2*xi*x + xi^2 ) +  wj*( x^2 + 2*xj*x + xj^2 )
(wi+wj)*x^2  +   (wi*xi + wj*xj)*2*x   +   wi*x^2 + wj*x^2

wij*( x^2   +  2*x*xij +   xij^2 )   +    C

wij = (wi+wj)
wij*xij =   (wi*xi + wj*xj)

xij =    (wi*xi + wj*xj) / (wi + wj)

C =   wi*xi^2 + wj*xj^2     -    wij*xij^2
C =   wi*xi^2 + wj*xj^2     -   (wi*xi + wj*xj)^2 / (wi + wj)
C = wi*wj*(xi**2 - 2*xi*xj + xj**2)/(wi + wj)
C = wi*wj*(xi-xj)**2/(wi+wj)

Derived & Tested in Python Here:
/home/prokop/git/SimpleSimulationEngine/python/pyGaussAtom/GaussProduct.py
*/





inline double product1D_w( double wi,double xi,   double wj,double xj,   double& W, double& X  ){
    W           =  wi+wj;
    double wxi  =  wi*xi;
    double wxj  =  wj*xj;
    double wx   =  wxi + wxj;
    X           =  wx/W;
    double logC =  wxi*xi + wxj*xj - wx*X;
    //double C   = np.exp(-logC) * Ci * Cj
    return logC;
}

inline double product3D_w( double wi, const Vec3d& pi, double wj, const Vec3d& pj,  double& wij, Vec3d& pij ){
    double junk;
    double logC =
        + product1D_w( wi,pi.x,  wj,pj.x,   wij ,  pij.x );
        + product1D_w( wi,pi.y,  wj,pj.y,   junk,  pij.y );
        + product1D_w( wi,pi.z,  wj,pj.z,   junk,  pij.z );

    /*
    Vec3d wxi   =  pi*wi;
    Vec3d wxj   =  pj*wj;
    double W    =  wi+wj;
    double wx   =  wxi + wxj;
    double X    =  wx/W;
    double logC =  wxi*xi + wxj*xj - wx*X;
    */
    return logC;
}


inline double product3D_s( double si, const Vec3d& pi, double sj, const Vec3d& pj,  double& sij, Vec3d& pij ){
    double junk;
    double wi = 1/(2*si*si);
    double wj = 1/(2*sj*sj);
    double wij;
    double logC =
        + product1D_w( wi,pi.x,  wj,pj.x,   wij ,  pij.x );
        + product1D_w( wi,pi.y,  wj,pj.y,   junk,  pij.y );
        + product1D_w( wi,pi.z,  wj,pj.z,   junk,  pij.z );
    sij = 1/sqrt(2*wij);
    //printf( "DEBUG product3D_s sij %g(%g,%g) wij %g(%g,%g)\n", sij,si,sj,  wij,wi,wj );
    return logC;
}


inline double norm3Dw( double w ){
    // https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function
    // n-dimansional Gaussian normalization C is :
    // G(x,y,z) = C*exp( -( x^2 + y^2 + z^2)/(2*sigma^2) ) =   C*exp( -w*( x^2 + y^2 + z^2) )
    // C  = 1/( sigma^n  (2*pi)^(n/2)  )
    // C3 = 1/( sigma^3 (2*pi)^(3/2) )
    // 1/sigma^2 = w
    // C3 = 1/( sqrt(2*sigma^2)^3/((2)^(3/2)) (2*pi)^(3/2) )
    // C3 = 1/( w^(3/2)/((2)^(3/2)) (2*pi)^(3/2) )
    // C3 = 1/( w*pi )^(3/2)
    double c2 = w/M_PI;
    double c  = sqrt( c2 );
    return c2*c;
}

inline double norm3Ds( double s ){
    const double sqrt_2pi = 2.50662827463;
    double c = 1/(sqrt_2pi*s);
    return c*c*c;
}

inline double sqnorm3Ds( double s ){
    const double sqrt_pi = 1.77245385091;
    // exp(-r2/(2*s^2))*exp(-r2/(2*s^2)) = exp(-r2/(s^2))
    // exp(-r2/(s^2)) = s*sqrt(pi)
    // ToDo: it may be better to store      1/sqrt(s) ... to avoid unnecessary sqrt() call
    double c = sqrt( 1/(sqrt_pi*s) ); // sqrt(sqrt(pi)*s)
    return c*c*c;
}

inline double sqnorm3Ds_sq( double s ){
    const double sqrt_pi = 1.77245385091;
    double c = 1/(sqrt_pi*s); // sqrt(sqrt(pi)*s)
    return c*c*c;
}

//inline double uy3Ds( double r, double s ){ return norm3Ds(s)* exp( (r*r)/(-2*s*s) ); }
//inline double uy3Ds2( double r2, double s ){ return norm3Ds(s)* exp( r2/(-2*s*s) ); }
inline double bas3D_r2( double r2, double s ){ return sqnorm3Ds(s)* exp( r2/(-2*s*s) ); }
inline double rho3D_r2( double r2, double s ){ return   norm3Ds(s)* exp( r2/(-2*s*s) ); }

inline double LoG3D_r2( double r2, double s ){ double wr2=r2/(-2*s*s); return   norm3Ds(s)* ( (1+wr2) )*exp( wr2 ); }



inline double kinetic_w( double w ){
    // !!!!! NOT NORMALIZED GAUSSIAN
    //    -((3*pi^(3/2))/(2^(3/2) )   /sqrt(w))
    return -5.90610372965 * sqrt(w);
}


inline double kinetic_w(  double r2, double w1, double w2 ){
    // Derived in Maxima: see : /home/prokop/Dropbox/MyDevSW/Maxima/Gauss_Kinetic-Cleaned.wxmx
    //                                                              Gauss_Kinetic-Polar-.wxmx
    // Cartes : (2*%pi^(3/2)*w1*w2*(2*w1*w2*x1^2-3*w2-3*w1)*%e^(-(w1*w2*x1^2)/(w2+w1)))/(w2+w1)^(7/2)
    // Polar  : (2*%pi^(3/2)*w1*w2*(2*w1*w2*x1^2-3*w2-3*w1)*%e^(-(w1*w2*x1^2)/(w2+w1)))/(w2+w1)^(7/2)
    //   (2*pi^(3/2) *   w1*w2*    (2*w1*w2* x1^2 -3*(w2+w1) )*  exp( -(w1*w2*x1^2)/(w2+w1))   )   /(w2+w1)^(7/2)
    //                   wprod*    (2 *wprod* x*x - 3*wsum )  *  exp( -( wprod *x*x )/wsum  )
    double C         = 11.1366559937; //(pi^(3/2))*2
    double wsum      = w2+w1 ;
    double iwsum     = 1/wsum;
    double sqrtiwsum = sqrt(iwsum);
    double wprod     = w1*w2 ;
    //double x2        = x*x;
    double g         = exp( -( wprod * r2 )*iwsum );
    return C*wprod * ( 2*wprod*r2 - 3*wsum )*g*iwsum*iwsum*iwsum*sqrtiwsum;
}

inline double kinetic( double s ){
    //     -(3*pi^(3/2)*s)/2
    //return -8.35249199525 * s * sqnorm3Ds(s) * sqnorm3Ds(s);
    return -8.35249199525 * s * sqnorm3Ds_sq(s);
}

inline double kinetic( double r2, double s1, double s2 ){
    return kinetic_w( r2, 1/(2*s1*s1), 1/(2*s2*s2) )  *  sqnorm3Ds(s1) * sqnorm3Ds(s2);
}


/// Boys Function
//  https://chemistry.stackexchange.com/questions/41214/boys-function-for-gaussian-integrals-in-ab-initio-calculations
//  BoysF is used for Coulomb integrals like :
//  Vpq = Integral{   exp(-p*(r1)^2) * exp(-q*(r2)^2)  /|r1-r2| }
//  Vpq = 2*pi^(5/2) / ( p*q*sqrt(p+q) ) * BoysF( a*r^2 )
//  F0() is related to error-function erf()
//  F0(x) = sqrt(pi/(4x)) * erf( sqrt(x) )
//  Vpq = (pi/sqrt(p*q))^3 * erf( sqrt(a) * r ) / r
// => It makes sense to tabulate rather F0(x^2)

}

#endif


