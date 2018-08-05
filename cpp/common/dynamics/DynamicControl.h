
#ifndef DynamicControl_h
#define DynamicControl_h

#include <cstddef>
#include <math.h>

//typedef void (*ForceFunction)( int n, double * xs, double * dfs );

class DynamicControl{ public:
    //double * x   = NULL; // pointer to input
    //double * y   = NULL; // pointer to output

    double x        =  0.0; // current value of control
    double dxdt_max =  1.0; // maximim speed of control movement
    double xmin     = -1.0; // minimum value of control
    double xmax     =  1.0; // maximum value of control

    double   y0 = 0.0;  // target value of output
    double   oy = 0.0;  // previous value o output

    double dydx = 1.0;  // derivative of output with respect to constrol;   this could eventually be recalculated for each target y0
    double T    = 0.0; // for second order dynamics

    double K    = 1.0;
    double damp = 0.8;

    // axuliary
    double vy,ovy;

    // ==== functions

    void setup( double y0_, double xmin_, double xmax_, double K_ ){ y0=y0_; xmin=xmin_; xmax=xmax_; K=K_; }


    double dx_O1( double y, double dt ){ // first order dynamical controler

        // dydt = 0;
        // TODO : we should modify this to consider observed output change rate  "dydt"  probely
        double dxmax = dxdt_max*dt;
        //double dymax = dxmax*dydx;

        //double dydt = (y-oy);  oy = y; // this is for second order

        double dy  =  y0-y;
        double dyT = (y-oy)*T;
        if( (dy*dyT)>0 ){
            dy-=dyT;
            if((dy*dyT)<0) dy=0.0;
        };

        double dx = _clamp( dy/dydx, -dxmax, dxmax);
        oy=y;
        //double dx = _clamp( dydt, -dymax,dymax ) / dydx;

        x += dx;
        //double x_ = _clamp( x+dx, xmin,xmax );
        //dx = x_-x;
        //x  = x_;

        return dx;
    };


    double dx_O2( double y, double dt ){ // first order dynamical controler
        //   y_ = y + vy*t + k*x*t**2
        //::::|  y(t+T) = y0
        //   vy*T + (y-y0) + k*x*T**2 = 0
        //   x = ((y0-y)-vy*T)/(k*T**2)

        double dy = y0-y;
        vy  = (y-oy)/dt;
        dy -= (vy*T)*damp;

        double x_ = dy/(T*T*K);
        double dx = x_-x;
        oy=y;ovy=vy;
        x+=dx;
        return dx;
    };

    void x_O1( double y, double dt ){ // first order dynamical controler
        //   y_ = y + vy*t + k*x*t**2
        //::::|  y(t+T) = y0
        //   vy*T + (y-y0) + k*x*T**2 = 0
        //   x = ((y0-y)-vy*T)/(k*T**2)

        double dy = y0-y;
        ovy=(y-oy)/dt;
        oy=y;
        //x=K*dy;

        x=_clamp( K*dy, xmin, xmax );


    };


};

#endif
