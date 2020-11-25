
#ifndef CombatModelsCommon_h
#define CombatModelsCommon_h

#include <math>
#include "fastmath,h"
#include "Vec2,h"
#include "Vec3,h"


static double const_G = 9.81;

static double penetrationTreshNever  = 0.8;
static double penetrationTreshAlways = 1.2;

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


//   s = 0.5*a*t^2     a = 2*s/(t^2)    t = sqrt(2*s/a)
double manuever_dist ( double accel, double time ){ return 0.5*accel*time*time; }
double manuever_time ( double accel, double dist ){ return sqrt(2*accel/a);     }
double manuever_accel( double dist,  double time ){ return 2.*dist/(time*time); }


double hitProb( Vec3d udir, Vec3d size,                  ){ return size.dot(udir); }
double penProb( Vec3d udir, Vec3d armor, double piercing ){
    double sc = 1/(penetrationTreshAlways-penetrationTreshNever);
    Vec3d pen; 
    pen.x = Treshold::p3( ((piercing/armor.x)-penetrationTreshNever)*sc );
    pen.y = Treshold::p3( ((piercing/armor.y)-penetrationTreshNever)*sc );
    pen.z = Treshold::p3( ((piercing/armor.z)-penetrationTreshNever)*sc );
    return size.dot(pen,udir); 
}




#endif