
#ifndef Forces_h
#define Forces_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"



#define COULOMB_CONST  14.3996448915f

#define RSAFE   1.0e-4f
#define R2SAFE  1.0e-8f
#define F2MAX   10.0f

inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double q ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*q*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceMorseQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double q, double alpha ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const double R2ELEC = 1.0;
    double r     = sqrt( dp.norm2()+R2SAFE );
    double expar = exp( alpha*(r-r0));
    //double E     = eps*( expar*expar - 2*expar );
    double ir    = 1/r;
    double fr    = eps*2*alpha*( expar*expar - expar ) + COULOMB_CONST*q/( r*r + R2ELEC );
    //printf( " %g -> %g | (%g,%g,%g) %g\n" , r, fr,  r0, eps,  q, alpha );
    f.add_mul( dp, fr*ir );
}

inline void addAtomicForceQ( const Vec3d& dp, Vec3d& f, double q ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double fr   = ( ir*q*-14.3996448915f )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceLJ( const Vec3d& dp, Vec3d& f, double r0, double eps ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceExp( const Vec3d& dp, Vec3d& f, double r0, double eps, double alpha ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double r    = sqrt(dp.norm2() + R2SAFE );
    double E    = eps*exp( alpha*(r-r0) );
    double fr   = alpha*E/r;
    f.add_mul( dp, fr );
    //f.add_mul( dp, 1/(dp.norm2()+R2SAFE) ); // WARRNING DEBUG !!!!
}

inline Vec3d REQ2PLQ( Vec3d REQ, double alpha ){
    double eps   = sqrt(REQ.y);
    double expar = exp(-alpha*REQ.x);
    double CP =    eps*expar*expar;
    double CL = -2*eps*expar;
    //printf( "REQ2PLQ: %g %g %g  ->  %g %g\n", REQ.x, eps, alpha,   CP, CL );
    return (Vec3d){ CP, CL, REQ.z };
}

inline Vec3d getForceSpringPlane( const Vec3d& p, const Vec3d& normal, double c0, double k ){
    double cdot = normal.dot(p) - c0;
    return normal * (cdot * k);
}

inline Vec3d getForceHamakerPlane( const Vec3d& p, const Vec3d& normal, double c0, double e0, double r0 ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    //printf(  " normal %g %g %g \n", normal.x, normal.y, normal.z );
    double cdot = normal.dot(p) - c0;
    double ir   = r0/cdot;
    double ir3  = ir*ir*ir;
    double f    = e0*(ir/r0)*ir3*(ir3-1);
    //printf( "%g %g %g %g %g %g %g \n", f, cdot, ir, ir3, e0, c0, r0  );
    return normal * f;
}

inline Vec3d getForceSpringRay( const Vec3d& p, const Vec3d& hray, const Vec3d& ray0, double k ){
    Vec3d dp; dp.set_sub( p, ray0 );
    double cdot = hray.dot(dp);
    dp.add_mul(hray,-cdot);
    return dp*k;
}

#endif
