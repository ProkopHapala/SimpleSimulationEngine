
#ifndef  Noise_h
#define  Noise_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

namespace Noise{

// ==== function declaration

void simplexNoise2D( const Vec2d& pos, Vec2d& dpos );
void warpNoise3R( const Vec2d& pos0, const Vec2d& rot, double fdown, double strenght, int n, Vec2d& dpos_ );

// ==== inline functions

inline double splineR2( double r2 ){
    return sq(1.0-r2);
    //if( (r2>0.01)|(r2>1.0) ){ return sq(1.0-r2); }else{ return 0;  }
}

inline double interR2( double ar2, double br2, double fa, double fb ){
    double wa = 1.0/(ar2+0.00001);
    double wb = 1.0/(br2+0.00001);
    return  ( fa*wa + fb*wb ) / (wa + wb);
}

inline double getR2( double dx, double dy ){
    return (dx*dx+dy*dy)*0.86602540378;
    //return (dx*dx+dy*dy)*1.5;
}

inline double getQuadruRand( Vec3d p, uint32_t seed ){
    //double h = ;
    double h = 0;
    h += (fhash_Wang( seed*53507 )-0.5) * p.x;
    h += (fhash_Wang( seed*55399 )-0.5) * p.y;
    h += (fhash_Wang( seed*56479 )-0.5) * p.z;
    h += (fhash_Wang( seed*57737 )-0.5) * p.x*p.x;
    h += (fhash_Wang( seed*61153 )-0.5) * p.x*p.y;
    h += (fhash_Wang( seed*56681 )-0.5) * p.x*p.z;
    h += (fhash_Wang( seed*59183 )-0.5) * p.y*p.x;
    h += (fhash_Wang( seed*60887 )-0.5) * p.y*p.y;
    h += (fhash_Wang( seed*61487 )-0.5) * p.y*p.z;
    h += (fhash_Wang( seed*65027 )-0.5) * p.z*p.x;
    h += (fhash_Wang( seed*64877 )-0.5) * p.z*p.y;
    h += (fhash_Wang( seed*64969 )-0.5) * p.z*p.z;
    return h;
}

inline double craterProfile( double r2, double fRidge, double hMin, double hMax ){
    double f2 = fRidge*fRidge;
    if( r2>f2 ){
        return ( (1/r2 - 1) )*f2;
    }else{
        return ( r2/f2 )-f2;
    }
}

inline double getCraterHeight( Vec3d p, int n, double scr, Vec3d* pos, double* sizes ){
    double h = 0.0;
    p.normalize();
    for(int i=0; i<n; i++ ){
        Vec3d d   = p - pos[i];
        double r2 = d.norm();
        double R2 = sq(scr*sizes[i]);
        if( r2<R2 ){
            //h += sizes[i]*(1-r2/R2);
            h += craterProfile( r2/R2, 0.5, 0.25, -0.75 ) * sizes[i];
        }
    }
    return h*4.0;
}

};

#endif
