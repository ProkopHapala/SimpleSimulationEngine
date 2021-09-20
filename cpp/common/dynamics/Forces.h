
#ifndef Forces_h
#define Forces_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

#define COULOMB_CONST  14.3996448915

#define RSAFE   1.0e-4f
#define R2SAFE  1.0e-8f
#define F2MAX   10.0f


void sum(int n, Vec3d* ps, Vec3d& psum){ for(int i=0;i<n;i++){ psum.add(ps[i]); } };

void sumTroq(int n, Vec3d* fs, Vec3d* ps, const Vec3d& cog, const Vec3d& fav, Vec3d& torq){
    for(int i=0;i<n;i++){  torq.add_cross(ps[i]-cog,fs[i]-fav);  }
    //for(int i=0;i<n;i++){  torq.add_cross(ps[i],fs[i]);  }
}

void checkForceInvariatns( int n, Vec3d* fs, Vec3d* ps, Vec3d& cog, Vec3d& fsum, Vec3d& torq ){
    cog =Vec3dZero;
    fsum=Vec3dZero;
    torq=Vec3dZero;
    double dw = 1./n;
    sum(n, ps, cog ); cog.mul(dw);
    sum(n, fs, fsum); //cog.mul(dw);
    sumTroq(n, fs, ps, cog, fsum*dw, torq );
}

inline double boxForce1D(double x, double xmin, double xmax, double k){
    double f=0;
    if(k<0) return 0;
    if(x>xmax){ f+=k*(xmax-x); }
    if(x<xmin){ f+=k*(xmin-x); }
    return f;
}

inline void boxForce(const Vec3d& p, Vec3d& f,const Vec3d& pmin, const Vec3d& pmax, const Vec3d& k){
    f.x+=boxForce1D( p.x, pmin.x, pmax.x, k.x);
    f.y+=boxForce1D( p.y, pmin.y, pmax.y, k.y);
    f.z+=boxForce1D( p.z, pmin.z, pmax.z, k.z);
}

inline double evalCos2(const Vec3d& hi, const Vec3d& hj, Vec3d& fi, Vec3d& fj, double k, double c0){
    double c    = hi.dot(hj) - c0;
    double dfc  =  k*-2*c;
    fi.add_mul(hj,dfc);
    fj.add_mul(hi,dfc);
    return k*c*c;
}

inline double evalCos2_o(const Vec3d& hi, const Vec3d& hj, Vec3d& fi, Vec3d& fj, double k, double c0){
    double c    = hi.dot(hj) - c0;
    double dfc  =  k*-2*c;
    double dfcc = -c*dfc;
    fi.add_lincomb( dfc,hj, dfcc,hi );
    fj.add_lincomb( dfc,hi, dfcc,hj );
    return k*c*c;
}

inline double evalCosHalf(const Vec3d& hi, const Vec3d& hj, Vec3d& fi, Vec3d& fj, double k, Vec2d cs ){
    Vec3d h; h.set_add( hi, hj );
    double c2 = h.norm2()*0.25d;               // cos(a/2) = |ha+hb|
    double s2 = 1-c2;
    //printf( "ang[%i] (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) c2 %g s2 %g \n", ig, ha.x,ha.y,ha.z,  hb.x,hb.y,hb.z,  h.x,h.y,h.z,   c2, s2 );
    double c = sqrt(c2);
    double s = sqrt(s2);
    cs.udiv_cmplx({c,s});
    double E         =  k*( 1 - cs.x );  // just for debug ?
    double fr        = -k*(     cs.y );
    // project to per leaver
    //c2 *=-2;
    //double lw     = 2*s*c;       //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    //double fra    = fr/(lbond[ib.a]*lw);
    //double frb    = fr/(lbond[ib.b]*lw);
    fr /= 2*c*s;  // 1/sin(2a)
    c2 *=-2*fr;
    Vec3d fa,fb;
    fi.set_lincomb( fr,h,  c2,hi );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    fj.set_lincomb( fr,h,  c2,hj );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );
    return E;
}



inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*qq*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceMorse( const Vec3d& dp, Vec3d& f, double r0, double eps, double beta ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const double R2ELEC = 1.0;
    double r     = sqrt( dp.norm2()+R2SAFE );
    double expar = exp ( beta*(r-r0) );
    //double E     = eps*( expar*expar - 2*expar );
    double fr    = eps*2*beta*( expar*expar - expar );
    //printf( " %g -> %g | (%g,%g,%g) %g\n" , r, fr,  r0, eps,  q, alpha );
    //printf( " r %g expar %g fr %g kqq %g a %g eps %g \n" , r, expar, fr, COULOMB_CONST*qq, alpha, eps );
    f.add_mul( dp, fr/r );
}

inline void addAtomicForceMorseQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq, double alpha ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const double R2ELEC = 1.0;
    double r     = sqrt( dp.norm2()+R2SAFE );
    double expar = exp( alpha*(r-r0));
    //double E     = eps*( expar*expar - 2*expar );
    double fr    = eps*2*alpha*( expar*expar - expar ) + COULOMB_CONST*qq/( r*r + R2ELEC );
    //printf( " %g -> %g | (%g,%g,%g) %g\n" , r, fr,  r0, eps,  q, alpha );
    //printf( " r %g expar %g fr %g kqq %g a %g eps %g \n" , r, expar, fr, COULOMB_CONST*qq, alpha, eps );
    f.add_mul( dp, fr/r );
}

inline double addAtomicForceQ( const Vec3d& dp, Vec3d& f, double qq ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double E    = COULOMB_CONST*qq*ir;
    double fr   = -E*ir2;
    f.add_mul( dp, fr );
    return E;
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

inline Vec3d REQ2PLQ( const Vec3d& REQ, double alpha ){
    //double eps   = sqrt(REQ.y);  // this is expected to be already done
    double eps   = REQ.y;
    double expar = exp(-alpha*REQ.x);
    double CP    =    eps*expar*expar;
    double CL    = -2*eps*expar;
    //printf( "REQ2PLQ: %g %g %g  ->  %g %g\n", REQ.x, eps, alpha,   CP, CL );
    return (Vec3d){ CP, CL, REQ.z };
}

inline Vec3d REnergyQ2PLQ( const Vec3d& REQ, double alpha ){
    return REQ2PLQ( {REQ.x, sqrt(REQ.y), REQ.z}, alpha );
}

inline Vec3d getForceSpringPlane( const Vec3d& p, const Vec3d& normal, double c0, double k ){
    double cdot = normal.dot(p) - c0;
    return normal * (cdot * k);
}

inline Vec3d getForceHamakerPlane( const Vec3d& p, const Vec3d& normal, double z0, double amp, double R ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    //printf(  " normal %g %g %g \n", normal.x, normal.y, normal.z );
    double cdot = normal.dot(p) - z0;
    double ir   = R/cdot;
    double ir3  = ir*ir*ir;
    double f    = amp*(ir/R)*ir3*(ir3-1);
    //printf( "%g %g %g %g %g %g %g \n", f, cdot, ir, ir3, e0, c0, r0  );
    return normal * f;
}

inline Vec3d getForceMorsePlane( const Vec3d& p, const Vec3d& normal, double amp, double R, double beta ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    //printf(  " normal %g %g %g \n", normal.x, normal.y, normal.z );
    double r       = normal.dot(p) - R;
    double expar   = exp( beta*r );
    double m_expar = 1-expar;
    double E       =  amp*m_expar*m_expar;
    double f       = -amp*m_expar*  expar * beta;
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
