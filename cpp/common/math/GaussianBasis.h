
#ifndef GaussianBasis_h
#define GaussianBasis_h

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
    printf( "DEBUG product3D_s sij %g wij %g \n", sij, wij );
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
    double c = 1/(2.50662827463*s);
    return c*c*c;
}

inline double sqnorm3Ds( double s ){
    // exp(-r2/(2*s^2))*exp(-r2/(2*s^2)) = exp(-r2/(s^2))
    // ToDo: it may be better to store      1/sqrt(s) ... to avoid unnecessary sqrt() call
    double c = sqrt( 1/(1.77245385091*s) );
    return c*c*c;
}

//inline double uy3Ds( double r, double s ){ return norm3Ds(s)* exp( (r*r)/(-2*s*s) ); }
//inline double uy3Ds2( double r2, double s ){ return norm3Ds(s)* exp( r2/(-2*s*s) ); }
inline double bas3D_r2( double r2, double s ){ return sqnorm3Ds(s)* exp( r2/(-2*s*s) ); }


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



