
#ifndef balistics_h
#define balistics_h

#include <math.h>
#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "appliedPhysics.h"

class Material{ public:
    double strength;
    double cShear;
    //double strength_shear;
    double density;
};

inline double caliber2area( double caliber ){ return M_PI*0.25*caliber; }
//inline double  projectileMass( double L, double caliber, double density );

//double flatIncidence( double v, double thickNor, caliber=100e-3, mass=None, density=7.8e+3, strength=2e+9, cShear=0.75, cDisp=1.0, debug=False ){

inline double flatIncidence( double v, double thick, double caliber, double mass, const Material& mat ){
    const double cDisp = 0.5;   // TODO : this should be espemated more carefully; e.g. from COG movement https://en.wikipedia.org/wiki/Reduced_mass
    double Ek      = kineticEnergy( v, mass );
    double S       = caliber2area(caliber);
    double Sshear  = M_PI*caliber*thick*mat.cShear;  if(Sshear<S) S = Sshear;
    double Estatic        = thick*S*mat.strength;
    double displaced_mass = thick*S*mat.density;
    double Edynamic       = cDisp * kineticEnergy( v, displaced_mass );
    double Eout           = Ek - Estatic - Edynamic;   // if( Ek<0 ) Ek = 0;
    return Eout;
}

inline double ricochetAngle_1( double v, double ld, double proj_strenght, double rho_p, double rho_t ){
    // A simple estimate of the minimum target obliquity required for the ricochet of a high speed long rod projectile
    // 1979 J. Phys. D: Appl. Phys. 12 1825; http://iopscience.iop.org/0022-3727/12/11/011
    double sqrtDensRatio = sqrt( rho_p/rho_t );
    double cDens  = 1.0 + sqrtDensRatio;
    double cShape = (ld*ld + 1)/(ld*1.0);
    double rhs    = (2.0/3.0)*((rho_p * v*v)/proj_strenght)*cShape*cDens;
    double tgA = pow( rhs, 1.0/3.0 );
    return tgA;
}

inline double surfaceDeflection( double v, double tgA, double caliber, double mass, double normalStress ){
    double S     = caliber2area(caliber);
    double aT    = (normalStress * S * tgA )/mass;
    double tdent = (caliber*tgA)/v;
    double vT    = tdent * aT;
    return vT;
}

inline double obliqueIncidence( double vr, Vec2d& rot, double thick, double caliber, double mass, const Material& mat ){
    // deflection of projectile by sloped armor ( increase angle )
    double normalStress = mat.strength;
    double tgA = rot.y/rot.x;
    double vT  = surfaceDeflection( vr, tgA, caliber, mass, normalStress );
    double dsA = vT/vr;
    dsA = 1-1/(dsA+1); // TODO : this is just stupid hack to make sure it is not too high
    Vec2d drot; drot.fromSin( dsA );
    rot.mul_cmplx( drot ); // modified rotation
    // calculate energy which pass the armor of given thickness
    double effective_thickness = thick/rot.x;
    double Eout = flatIncidence( vr, effective_thickness, caliber, mass, mat );
    return Eout;
}

inline Vec3d penetrate( Vec3d v, Vec3d nor, double thick, double caliber, double mass, const Material& mat ){
    double vr   = v.norm();
    double cdot = nor.dot(v);
    Vec2d rot; rot.fromCos(cdot/vr);
    //double Eout = flatIncidence( vr, thick, caliber, mass, mat );
    double Eout = obliqueIncidence( vr, rot, thick, caliber, mass, mat );
    if( Eout > 0 ){ // calculate velocity vector after penetration
        double vro = sqrt(2*Eout/mass);
        Vec3d vout = nor*(rot.x*vro) + (v-nor*cdot)*(rot.y*vro/vr);
    }
    return (Vec3d){0.0,0.0,0.0};
}

#endif
