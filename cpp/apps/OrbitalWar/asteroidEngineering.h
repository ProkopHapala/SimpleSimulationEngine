
#ifndef  asteroidEngineering_h
#define  asteroidEngineering_h

/*

This deals with economy of asteroides

*/

#include <vector>
#include <unordered_map>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "geom3D.h"
//#include "Body.h"

#include "appliedPhysics.h"
#include "CommodityNetwork.h"
#include "SpaceBodies.h"

static double lambda0=0;
static double spectrumStep=0;

struct SpectralProperties{
    int n;
    double* relfectance;

    double getReflectance(double logLambda){
        double x = (logLambda-lambda0)*spectrumStep;
        int i=(int)x;
        x-=i;
        return relfectance[i+1]*x-relfectance[i]*(1-x);
    }
};

struct MineralType{
    std::string name;
    std::unordered_map<int,double> elements;  // elemntal composition
};

struct RockType{
    std::string name;
    std::unordered_map<MineralType*,double> minerals;
    double hardness;        // How hard it is to cut
    double heatOfMelting;
    double ablationEnergy;  // How hard it is to evaporate
    SpectralProperties spectal;
};

struct Deposit{
    //RockType* type;
    double amount;
};


double manuever_planeChange( double angle, double mass, Orbit* orbit, double& mProp, double& deltaV ){
    // --- deltaV change needed
    double sa     = sin(angle/2);
    double r      = orbit->apoapsis();
    deltaV = sa*orbit->L/r;
    // ---- engine performance
    double T         = 3000.0; // K
    double molarMass = 0.018;  // water vapor
    double vexh      = exhaustVelocity( T, molarMass, 0.8, const_heatCapacityRatio_watter );
    double payload   = tsielkovsky_payload( deltaV, vexh );
    //printf( "payload %g fuel %g deltaV %g vexh %g \n", payload, 1-payload, deltaV, vexh );
    mProp            = mass * ( 1 - payload );
    double nMol      = mProp/molarMass;
    double E         = nMol*const_Rgas*T;
    return E;
}

class Asteroide{
    Orbit* orbit;
    double mass;
    std::unordered_map<RockType*,Deposit> rocks;
    std::unordered_map<RockType*,Deposit> minerals;

    double planeChange( double angle, double& mProp, double& deltaV ){
        return manuever_planeChange( angle, mass, orbit, mProp, deltaV );
    }

};


#endif
