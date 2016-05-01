
#ifndef  Noise_h
#define  Noise_h

#include "fastmath.h"
#include "Vec2.h"

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

};

#endif
