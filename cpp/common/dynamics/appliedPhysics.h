
#ifndef appliedPhysics_h
#define appliedPhysics_h

#include "fastmath.h"

inline double gunEnergy( double d, double l, double mPowder, double mProjectile, double Vchamber ){
    // http://www.stardestroyer.net/Armour/ShepStuff/Website/Statistics/MilitaryPropellants.htm
    // https://en.wikipedia.org/wiki/Adiabatic_process
    // W = alpha * p1*V1*( (V2/V1)^(1-gamma) - 1 )
    // W = alpha * p1*V1*( (V2/V1)^(1/alpha) - 1 )      // alpha is number of degrees of freedom  gamma = (alpha+1)/alpha
    const double alpha = 1.5;
    const double energyDensity = 4.6e+6;
    double U          = mPowder * energyDensity;
    //double p0         = U/(alpha*Vchamber);
    double V2         = Vchamber + d*l;
    //double Wtot       = alpha*p0*Vchamber*( pow( (V2/Vchamber), 1/alpha ) - 1 );
    double workRatio  = pow( (V2/Vchamber), 1/alpha ) - 1;
    double Wtot       = U*workRatio;
    double efficiency = mProjectile/(mProjectile+mPowder);
    return Wtot*efficiency;
};

#endif

