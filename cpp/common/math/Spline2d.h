
#ifndef  Spline2d_h
#define  Spline2d_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "geom2D.h"
#include "spline_hermite.h"

struct SplinePoint2d{ public:
    Vec2d p ; // position
    Vec2d dp; // derivative from left
    Vec2d dm; // derivative from right
};

class Spline2d{ public:
    std::vector<SplinePoint2d> CPs;

    inline Vec2d val( int i, double ds ){
        double c0,c1,d0,d1;
        Spline_Hermite::basis( ds, c0,c1,d0,d1 );
        //printf( "[%i,%g] c %g,%g d %g,%g \n", i, ds, c0,c1, d0, d1 );
        const SplinePoint2d& cp0 = CPs[i                              ];
        const SplinePoint2d& cp1 = CPs[wrap_index_fast(i+1,CPs.size())];
        return cp0.p *c0
             + cp1.p *c1
             + cp0.dp*d0
             + cp1.dm*d1;
    }

    inline bool val( double s, Vec2d& p ){
        if((s<0)||(s>CPs.size())){ return false; };
        int i=(int)s;
        p = val( i, s-i );
        return true;
    }

    void evalDerivs( double sc ){
        int n=CPs.size();
        for(int i=0; i<n; i++){
            Vec2d pm = CPs[wrap_index_fast(i-1,n)].p;
            Vec2d pp = CPs[wrap_index_fast(i+1,n)].p;
            pp.sub(pm);
            pp.mul(sc);
            CPs[i].dp=pp;
            CPs[i].dm=pp;
        }
    }

};


class Spline2dSampler{ public:
// ### To preload points and make it faster
//    prepare();
};

#endif



