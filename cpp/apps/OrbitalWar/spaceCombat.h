
#ifndef  spaceCombat_h
#define  spaceCombat_h

#include <math.h>
//#include "fastMath.h"
#include "appliedPhysics.h"

/*

Brain-Storm
===========

Weapon Types
 - Mass-Drivers
    solid-state projectiles (typically made of heavy metal) with high penetration and high impact on target. However preccision is much worse than lasers flight time is very long, and require very long acceleration line
    - Magnetic - provide optimal energy efficieny (60-80%) and optimal accuracy, however the acceleration force is limited by magnetic fields
    - Ablation - Ablation of projectile by neutralized ion beam provides lower efficieny, limited by yeld strenght
    - Electro-spray - microscopic coloide particles ()
 - Plasma
    - plasmoides
 - Neutralized Ion Beam -
 - Nuclear missiles


*/


//double kineticEnergy( double force, double length ){ return force * length; }
//double kineticEnergy( double dist, double time ){ return force * length; }



/*

Ablation Accelerated projectile
 - limited by yield-strength of projectile (2GPa)

*/



// ================ Laser

double const_LightSpeed = 3e+8; //[m/s]

double difractionLimit_spot( double wavelenght, double aperture, double distance ){
    return distance*wavelenght/aperture;
}

double difractionLimit_intensity( double wavelenght, double aperture, double distance ){
    //double r = difractionLimit_spot( wavelenght, aperture, distance );
    //return 1/(r*r);
    double invr = aperture/(distance*wavelenght);
    return invr*invr;
}

// ================ Kinetic


double kineticDispersion( double v0, double m0, double mShield ){
    // p = v0*m0 =    v1*(m0 + mShield)
    double m1 = m0+mShield;
    double mratio = m0/m1;
    double v1 = v0*mratio;
    //double dE = 0.5*( m0*v0*v0 - m1*v1*v1 );
    //double dv = sqrt( ( m0*v0*v0 - m1*v1*v1 )/m1 );
    //double dv = sqrt( ( m0/m1)*v0*v0 - v1*v1 );
    //double dv = sqrt( mratio*v0*v0 - v0*v0*mratio*mratio );
    double dv = sqrt( (v0 - v1)*v1 );
    printf( "v0 %g v1 %g dv %g\n", v0, v1, dv );
    double tg = dv/v1;
    return tg;
}


double diskArea(double radius){
    return M_PI * radius * radius;
}

double diskVolume(double radius, double thick){
    return diskArea(radius) * thick;
}

double accelTime( double accel, double length ){
    return sqrt( 2*length / accel );
}

double pressureLimitedAccelerator( double massPerArea, double maxPressure, double length ){
    double accel = maxPressure / massPerArea; // (Force/area) / (mass/area)
    // s = 0.5 * a * t^2   ;   t = sqrt( 2*s/a )
    double t     = accelTime( accel, length );
    return accel * t;
}

// some limit on gas exhoust velocity (from adiabatic expansion equation)
//double ( double vexh ){
//    return k/vexh;
//}

double rocketEquation( double vexh, double massRatio ){
    return vexh * log( massRatio );
}

double rocketEquation( double massRatio ){
    return exp( massRatio ); // v/vexh
}

//double powerLimitedAcceleration( double d ){}


const double const_Rgas_SI = 8.3; // J/(K*mol)

// pV  = nRT
// p   = nRT/V = (n/S)*R*T/l
// P   = F*v =
// P/S = (F/S)*v = p * v
// v   = sqrt(2*Ek/M) = sqrt( 2*(3/2)RT/M ) = sqrt(3*RT/M)


// dp = vech * dm
// F  = dp/dt = m*dv/dt = vexh*dm/dt
// m*dv = vexh*dm
// p  = (m/S)*dv  =  vexh*(dm/S)

double adiabaticExpansionEnergy( double pressureRatio, double T0, double kappa ){
    // ToDo Check this !!!!!!!!!!!!!!!!!!!!!
    double e = (kappa-1)/kappa;
    return const_Rgas_SI * T0 * ( 1-pow( pressureRatio, e ) )/e; // energy per mol
}

//double adiabaticExpansionEnergy( double pressureRatio, double T0, double kappa ){
//    return pow( pressureRatio, e );
//}

double meanThermalVelocity(double temparature, double molarMass){
    return sqrt( 3 * const_Rgas_SI * temparature / molarMass );
}

double rocketMassFlow( double force, double vexh ){
// F  = dp/dt
//    = m*dv/dt = vexh*(dm/dt)
// p  = (m/S)*dv  =  vexh*(dm/S)
//  dm = F
    return force/vexh;
}







double timeToTarget(double velocity, double dist){
    return dist / velocity;
}

double accelToSpread(double accel, double time){
    return 0.5 * accel * ( time * time );
}

double dist( double accel, double velocity, double dist ){
    return accelToSpread( accel, timeToTarget( velocity, dist ) );
}


#endif
