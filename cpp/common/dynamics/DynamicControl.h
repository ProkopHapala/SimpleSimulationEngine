
#ifndef DynamicControl_h
#define DynamicControl_h

#include <cstddef>
#include <math.h>

//typedef void (*ForceFunction)( int n, double * xs, double * dfs );

class DynamicControl{ public:
    //double * x   = NULL; // pointer to input
    //double * y   = NULL; // pointer to output
    double   y0  = 0.0;  //target value of output
    double   oy  = 0.0;  // previous value o output

    double vmax = 1.0;   // maximal rate input change ... max(dx/dt)
    double dydx = 1.0;   // this could eventually be recalculated for each target y0

    double dx_O1( double y, double dt ){ // first order dynamical controler
        double dydt = (y-oy);
        // TODO : we should modify this to consider observed output change rate  "dydt"  probely
        double dxmax = vmax*dt;
        double dymax = dxmax*dydx;

        return _clamp( y0-y, -dymax,dymax ) / dydx;
    };

};

#endif
