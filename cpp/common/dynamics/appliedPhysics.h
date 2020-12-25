
#ifndef appliedPhysics_h
#define appliedPhysics_h

#include "fastmath.h"
#include "Vec3.h"

static const double const_GravAccel              = 9.81;
static constexpr const double const_Graviational = 6.6743015e-11;
static constexpr const double const_Rgas         = 8.31446261815324;

static constexpr const double const_AU = 149.597870700e+9; // [m] Astronomic Unit
// second in year   3.17098e-8 kW year

// Specific energy units : https://en.wikipedia.org/wiki/Energy_density
static constexpr const double const_EkgAnihil = 89.875517874e+15; // Energy from anihilation
static constexpr const double const_EkgDT     = 337.387388e+12;   // Energy of one kg of Deuterium-Tritium fusion
static constexpr const double const_EkgU      = 80e+12;           // Energy of one kg of Uranium (Plutonium, Thorium) completely fissioned //
static constexpr const double const_EktTNT    = 4.184e+12;        // kiloton of TNT
static constexpr const double const_EMWy      = 31.54e+12;        // MegaWatt Year

static constexpr const double const_heatCapacityRatio_monoatimic = 1.6666666;
static constexpr const double const_heatCapacityRatio_diatomic   = 1.4;
static constexpr const double const_heatCapacityRatio_watter     = 1.3333333;

static constexpr const double const_SolarRadEarth =  1366.1; // [W/m^2]


static constexpr const double const_hour  =  3600; // [W/m^2]
static constexpr const double const_day   =  86400; // [W/m^2]
static constexpr const double const_month =  2592000; // [W/m^2]
static constexpr const double const_year  =  31536000; // [W/m^2]

char* timeInfo(char* s, double t_sec){
    if(t_sec<const_hour  ) return s+sprintf(s,"%g s"     ,t_sec);
    if(t_sec<const_day   ) return s+sprintf(s,"%g hours" ,t_sec/const_hour );
    if(t_sec<const_month ) return s+sprintf(s,"%g days"  ,t_sec/const_day  );
    if(t_sec<const_year  ) return s+sprintf(s,"%g months",t_sec/const_month);
    return s+sprintf(s,"%g years" ,t_sec/const_year);
}





inline double solarRadDist_SI( double r ){ return const_SolarRadEarth * sq( const_AU/ r ); }


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

double tsielkovsky_speed( double payload, double vexh ){
    return vexh*log( 1/payload );
}

double tsielkovsky_payload( double deltaV, double vexh ){
    return exp( -deltaV/vexh );
}


double jetEfficiency( double expansionRatio, double kappa=const_heatCapacityRatio_monoatimic ){
    double efficiency =  pow( expansionRatio, (1-kappa)/kappa ) - 1 ;
    //double a = 1/(kappa-1);
    //return a * efficiency;
    return efficiency;
}

double exhaustVelocity( double T, double molarMass=1., double efficiency=1., double kappa=const_heatCapacityRatio_monoatimic  ){
    //https://en.wikipedia.org/wiki/Adiabatic_process
    // (a+1)/a = kappa   =>   1/a = kappa - 1 =>  a = 1/(kappa-1)
    // W = - alpha * n*R*T1 * ( expansionRatio^(1-kappa) - 1 );
    double alpha = 1/(kappa-1);
    double W = alpha * efficiency*const_Rgas*T/molarMass;
    double v = sqrt(W*2);
    return v;
}

#endif

