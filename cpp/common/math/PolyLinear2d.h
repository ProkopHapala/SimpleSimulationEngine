#ifndef  PolyLinear1d_h
#define  PolyLinear1d_h

// Various math routines for piecewise-linear functions (grids, splines) in 2D F(x,y)

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"

namespace Lin2Dcross{
//  symetrized split of square by both diagonals, centre point is put in the center with average of the boundary points  

double val( double x, double y, double f00, double f01, double f10, double f11 ){
    double fc = 0.25*(f00 + f01 + f10 + f11);
    if(x>y){
        if(x+y>1){ // right
            return 2*( (1-x)*fc + (x-0.5)( y*f11 + (1-y)*f01 ) );
        }else{     // down
            return 2*( y*fc + (0.5-y)( x*f01 + (1-x)*f00 ) );
        }
    }else{
        if(x+y>1){ // up
            return 2*( (1-y)*fc + (y-0.5)( x*f11 + (1-x)*f10 ) );
        }else{     // left
            return 2*( x*fc + (0.5-x)( y*f10 + (1-y)*f00 ) );
        }
    }
}

};

#endif
