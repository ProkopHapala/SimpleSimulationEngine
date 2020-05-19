
#ifndef GaussianBasis_h
#define GaussianBasis_h

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

inline double gaussProduct1D( double wi,double xi,   double wj,double xj,   double& wij, double& xij  ){
    double W    =  wi+wj;
    double wxi  =  wi*xi;
    double wxj  =  wj*xj;
    double wx   =  wxi + wxj;
    double X    =  wx/W;
    double logC =  wxi*xi + wxj*xj - wx*X;
    //double C   = np.exp(-logC) * Ci * Cj
    return logC;
}

inline double gaussProduct3D( double wi, const Vec3d& pi, double wj, const Vec3d& pj,  double& wij, Vec3d& pij ){
    double junk;
    double logC =
        + gaussProduct1D( wi,pi.x,  wj,pj.x,   wij ,  pij.x );
        + gaussProduct1D( wi,pi.y,  wj,pj.y,   junk,  pij.y );
        + gaussProduct1D( wi,pi.z,  wj,pj.z,   junk,  pij.z );

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

inline double gaussNorm3D( double w ){
    // https://math.stackexchange.com/questions/434629/3-d-generalization-of-the-gaussian-point-spread-function
    // n-dimansional Gaussian normalization C is :
    // G(x,y,z) = C*exp( -( x^2 + y^2 + z^2)/(2*sigma^2) ) =   C*exp( -w*( x^2 + y^2 + z^2) )
    // C  = 1/( sigma^n  (2*pi)^(n/2)  )
    // C3 = 1/( sigma^3 (2*pi)^(3/2) )
    // 1/sigma^2 = w
    // C3 = 1/( sqrt(2*sigma^2)^3/((2)^(3/2)) (2*pi)^(3/2) )
    // C3 = 1/( w^(3/2)/((2)^(3/2)) (2*pi)^(3/2) )
    // C3 = 1/( w*pi )^(3/2)
    double c2 = 1/(M_PI*w);
    double c  = sqrt( c2 );
    return c2*c;
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

#endif



