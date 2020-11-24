
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



/*

Simple Damage Model

- Vehicle is composed o modules
- Modules are esential for function ( performance of whole system is typically product of function of modules )
- Modules have certain 
    1] crossection and 
    2] certain armor
    both of which can by in simple terms projected into cartesina direction x,y,z
      =>  hit corssection = area = dot( S, dir );        armor = dot( A, dir )
-If the shot pass through it destroys certain volume (proportional to cross area) of the system
   - some systems react well to damage ( e.g. radiators )
   - some systems react badly (even small damage makes it disfunctional)
     - this can be expressed by some function - damage exponent
     decrease of performace scales as (1-x)^alha  or 1-x^alpha 
 


*/


#endif