
#ifndef AirCombatModel_h
#define AirCombatModel_h

#include  "AeroSurf.h"

/*

Fd = Cd * (rho/2) * v^2
Fl = Cl * (rho/2) * v^2
Fg = m * g;


Fd + 

area[m^2] * v[m/s] * rho[kg/m^3]  =  mass_flow [kg/s]
E = (v2^2 - v1^2)*m
P = 0.5*(v2^2 - v1^2)*mass_flow
P = 0.5*area*rho *  v1*(v2^2-v1^2)
P = 0.5*area*rho *  v2*(v2^2-v1^2)
P = 0.5*area*rho * ( v0 + v )*( v0^2 + 2*v0*v + v^2 - v0^2 )
P = 0.5*area*rho * ( v0 + v )*( 2*v0*v + v^2 )
P/(area*rho)

F = (v2-v1)*mass_flow

*/

inline double getPoperTrust(double power, double v0, double& Dv, double area_dens ){
    // From change of kinetic energy of flowing air
    // P = 0.5*area*rho *  v0*( (v0+Dv)^2 - v0^2 )
    // P = 0.5*area*rho *  v0*( v0^2 + 2*v0*Dv + Dv^2 - v0^2 )
    // P = 0.5*area*rho *  v0*( 2*v0*Dv + Dv^2  )
    // P = area*rho *  v0*( v0*Dv + 0.5*Dv^2  )
    double dm = area_dens*v0;
    double Dv_junk;
    quadratic_roots( 0.5*dm, dm*v0, -power,     Dv_junk, Dv );
    return dm * Dv;
}

inline void makeClimbCurves( int n, double amin, double amax, const AeroPolar& polar, double area, double mass, double speed0, double power, double airDens, double gravConst ){
    double area_dens = area*airDens;
    double dAng = (amax-amin)/n; 
    for(int i=0; i<n; i++){
        double ang = amin + dAng*i;
        double Cl,Cd;
        double ca = cos(ang);
        double sa = sin(ang);
        polar.getLD( ca, sa, Cd, Cl );
        double Fr = area * airDens * speed0*speed0;
        //double Ft = power/speed0;
        double Dv;
        double Ft = getPoperTrust( power, speed0, Dv, area_dens );
        double Fd = Cd * Fr + Ft; // = 0  
        double Fl = Cl * Fr;
        double Fx = Fd*ca - Fl*sa;
        double Fy = Fd*sa + Fl*ca + gravConst*mass;  // =0

        
    }
}


#endif