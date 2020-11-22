
#ifndef CombatModelsCommon_h
#define CombatModelsCommon_h

#include <math>

double const_G = 9.81;

double timeOfFlightParabolic( double distance, double vMuzzle ){
    // parabolic trajectory
    // x = v0*ca*t
    // y = v0*sa*t - 0.5*g*t^2    = 0
    // (v0*sa)/(0.5*g) = t
    // x = v0*ca*(v0*sa)/(0.5*g)
    // x = (v0^2/g) * (2*ca*sa)
    // x/(v0^2/g) = sin(b) 
    double sin_beta = distance*const_G/vMuzzle;
    double alpha    = asin( sin_beta )/2;
    double  
}

#endif