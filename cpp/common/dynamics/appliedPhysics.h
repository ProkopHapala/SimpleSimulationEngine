
#ifndef appliedPhysics_h
#define appliedPhysics_h

#include "fastmath.h"
#include "Vec3.h"

static const double  const_Graviational = 6.674e-11;




inline double kineticEnergy( double v, double mass ){ return 0.5*mass*v*v; }

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

inline double armorThicnessFactor_MomentumModel( double velocityTangent, double massFactor ){
//    # https://panzerworld.com/relative-armor-thickness, NIKO HOLKKO, MECHANISMS OF ARMOUR PENETRATION, Bachelor's Thesis, 2015
//    ### Model 1 - momentum refraction
//    * Distance Traveled by projectil in armor can be computed by pyctagorean theorem from velocity vII paralel and vT perpendicular to armor plat ( vec{v}=(vII,vT) )
//       D = D0 * sqrt( 1 + (vII/vT)^2 )
//    where D0 is nominal(perpendicular) thinkcness of armor
//    * We assume that upon impact the perpendicular component vT is modified accoding to redistibution of momentum, while vII is unchanged
//       vT = vT0 * ( m/(m+M) )
//    where m is mass of projectile and M is effective mass of relevant armor segment
//    * The effective armor mass M can be calculated assuming some effective contact area S
//      M = rho * D0 * S
//    NOTE 1): contact area S can be interpreted as area of plate deformed (bulged) by impacting projectile
//    NOTE 1): Geometric thicness of sloped armor ( well known cosine rule ) correspond to M=0
    return sqrt( 1 + sq( (massFactor+1) * velocityTangent ) );
}


inline Vec3d centralGravityForce( const Vec3d& d, double Mm ){
    double r2 = d.norm2();
    //printf( "centralGravityForce %f %f \n", r2, Mm );
    return d * (  Mm * const_Graviational / ( r2 * sqrt(r2) ) );
}



#endif

