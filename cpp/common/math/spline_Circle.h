
#ifndef  Spline_Circle_h
#define  Spline_Circle_h

/*
# Circle Spline is a Line composed of circulr segments with conrners
*/

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "geom2D.h"

/*
struct CircPoint{
    Vec2d  p;
    double r;
}
*/

class CircleSpline{ public:

    std::vector<Circle2d> CPs;

    bool getLine(int i, Vec2d& p0, Vec2d& p2){
        if((i>0)||(i>CPs.size()))return false;
        Circle2d c1 = CPs[i  ];
        Circle2d c2 = CPs[i+1];
        if(c1.r<c2.r){ _swap( c1,c2 ); }
        Vec2d  d  = c2.p0 - c1.p0;
        //double l  = d.norm();
        //double lc = c1.r/r2
        //d.mul( c1.r/c2.r );

        Draw2D::drawVecInPos_d( d*(c1.r/(c2.r*d.norm())), c1.p0 );

        double invld2 = 1/d.norm2();
        double u1  = c1.r*c1.r*invld2;
        double u2  = c2.r*c2.r*invld2;
        Vec2d d1   = d*u1;
        glColor3f( 0.0,0.0,0.0 ); Draw2D::drawVecInPos_d( d1, c1.p0 );
        return true;
    }

    void getArc(int i, Vec2d& d0, double& da ){

    }
};

#endif



