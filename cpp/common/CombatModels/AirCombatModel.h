
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

static const double const_airDensGround = 1.22;
static const double const_GravAccel     = 9.81;

struct CombatAirCraft{
    double mass;        // [kg]
    double power;       // [W]   engine power
    double area_wing;   // [m^2] wing area - makes
    double area_hull;   // [m^2] hull area - makes only drag
    double aspect;      // [1]   aspec ratio of wing - affetcs predominantly lift-induced drag (Oswald efficiency factor included in)
    double CD0;         // [1    min Drag Coef at zero lift
    double CLmax;       // [1]   max Lift coef at stall
    double accelMax;    // [m/s^2]  maximum acceleration

    double area_total()const{
        return area_wing + area_hull;
    }

    double maxSpeed_simple()const{
        // max speed ignoring gravity / lift
        // P = F*v
        // F = 0.5 * air_dens * CD0 * area_tot * v^2
        // P = v^3 * (  0.5*air_dens * CD0 * area_tot  )
        double area = area_total();
        return pow( power / ( 0.5 * const_airDensGround * area * CD0 ) ,   1./3.  );
    }

    double trunRate_simple( double& CL, double& speed, double& accel, double& Rturn )const{
        // FDrag = FLift
        // CD = CD0 +    CL^2 * ( pi*e/aspect );   // https://en.wikipedia.org/wiki/Aspect_ratio_(aeronautics)#Details
        // P =  v^3 * (  0.5*air_dens * area_tot * CD(CL)  )
        // P =  v^3 * ( CL^2 * ( pi*e/aspect ) * area_wing + CD0 *area_tot ) * 0.5*air_dens
        // FLift = v^2 *  CL * area_wing * 0.5*air_dens
        // FLift = ( P/(( CL^2 * ( pi*e/aspect ) * area_wing + CD0 *area_tot) * 0.5*air_dens ) )^(2/3) * CL * area_wing * 0.5*air_dens
        // FLift = c1*Cl*( c2 /( c3*Cl^2 + c4 )^(2/3)
        // c1 =    0.5*air_dens
        // c2 = P/(0.5*air_dens)
        // c3 =  ( pi*e/aspect ) * area_wing
        // c4 =  CD0             * area_tot
        // solve (c1 c2^(2/3) (3 c4 - c3 x^2))/(3 (c4 + c3 x^2)^(5/3)) = 0
        // x12 = +- sqrt(3*c4/c3)
        double area_tot = area_total();
        double c3     = area_wing/(M_PI*aspect);
        double c4     = CD0 * area_tot;
        CL     = sqrt(3*c4/c3);   // Best Lift
        CL     = fmin( CL, CLmax );
        double P_Fd = power/( 0.5*const_airDensGround * ( CL*CL*c3 + c4 ) );
        speed  = pow( P_Fd,  1./3. );
        double Flift  = speed*speed * CL * area_wing * 0.5 * const_airDensGround;
        accel  = Flift/mass;
        //printf( "CL %g c3 %g c4 %g speed %g P_Fd %g Fl %g[kN] accel %g[g] \n", CL, c3, c4, speed, P_Fd, Flift*1e-3, accel/const_GravAccel );
        Rturn  = (speed*speed)/accel;
        return speed/Rturn;  // Omega
    }

    double climbRate_simple( double& sa )const{
        // ToDo : We need to consider dependence of Propeler Efficiency on speed to get more realistic results
        // ToDo : Consider limited lift coeff CL (stall)
        // A1  = 0.5 * air_dens * area_wing
        // A2  = 0.5 * air_dens * area_tot
        // Fg  = m*g = CL * A1 * v^2
        // CL  =  m*g / ( v^2 * A1 )
        // P   =       v^2 * ( CD0*A2 + ( CL^2             /(pi*a) )*A1 )   +   m*g * vy
        // P   =       v^2 * ( CD0*A2 + ( (m*g /(v^2*A1))^2/(pi*a) )*A1 )   +   m*g * vy
        // vy  = ( P - v^2 * ( CD0*A2 - ( (m*g /(v^2*A1))^2/(pi*a) )*A1 )  ) / (m*g);
        // vy  =  c1 - c2*(v^2) - c3/(v^2)
        // c1  = P/(m*g)
        // c2  = (CD0*A2*)/(m*g)
        // c3  = (m*g)/(A1*pi*a)
        double A1 = 0.5 * const_airDensGround * area_wing;
        double A2 = 0.5 * const_airDensGround * area_total();
        double mg = mass * const_GravAccel;
        double c1 = power/mg;
        double c2 = CD0*A2/mg;
        double c3 = mg/(A1*M_PI*aspect);
        //double vy =  (2 (c3 - c2 v^4))/v^3;
        double v  = pow( c3/c2,  0.25 );
        double v2 = v*v;
        double vy = c1 - c2*v2 - c3/v2;
        sa = vy/v;
        //double ca = sqrt( 1 - sa*sa );
        return vy;  // Omega
    }

    double climbRate_CLmax( double& sa )const{
        // ToDo : We need to consider dependence of Propeler Efficiency on speed to get more realistic results
        // ToDo : Consider limited lift coeff CL (stall)
        // A1  = 0.5 * air_dens * area_wing
        // A2  = 0.5 * air_dens * area_tot
        // Fg  = m*g = CL * A1 * v^2
        // CL  =  m*g / ( v^2 * A1 )
        // v =  sqrt( m*g / ( CL * A1 ) )
        // P   =       v^2 * ( CD0*A2 + ( CL^2             /(pi*a) )*A1 )   +   m*g * vy
        // P   =       v^2 * ( CD0*A2 + ( (m*g /(v^2*A1))^2/(pi*a) )*A1 )   +   m*g * vy
        // vy  = ( P - v^2 * ( CD0*A2 - ( (m*g /(v^2*A1))^2/(pi*a) )*A1 )  ) / (m*g);
        // vy  =  c1 - c2*(v^2) - c3/(v^2)
        // c1  = P/(m*g)
        // c2  = (CD0*A2*)/(m*g)
        // c3  = (m*g)/(A1*pi*a)
        double A1 = 0.5 * const_airDensGround * area_wing;
        double A2 = 0.5 * const_airDensGround * area_total();
        double mg = mass * const_GravAccel;
        double v2 = mg/(CLmax*A1);
        double v  = sqrt( v2 );
        double c1 = power/mg;
        double c2 = A2*CD0/mg;
        double c3 = A1*(CLmax*CLmax)/(M_PI*aspect);
        double vy = c1 - c2*v2 - c3/v2;
        sa = vy/v;
        //double ca = sqrt( 1 - sa*sa );
        return vy;  // Omega
    }


};


//  Slip Away
// ------------------
//   Consider situation when fast aircaft try to shoot on slower but more maneuverable airceft
//   The target aircraft (pray) can turn inside envelope of the faster aircraft
//   Consider the slower aircaft turns at its maximum turn speed
//   While the faster airfact turns does this as well to compensate
//   The trajectory of the two is
//  x1 = R1*sin(O1*t)         y1 = R1*(1-cos(O1*t))
//  x2 = R2*sin(O2*t) + L0    y2 = R2*(1-cos(O2*t))
//  we should care about the angle of attack in which the hunger see the pray, this is limited
//  alpha < O1*t
//  tg(alpha) =   ( R1*(1-cos(O1*t)) - R2*(1-cos(O2*t)) ) /  ( L0 + R2*sin(O2*t) - R1*sin(O1*t) )
// So the final condition is :
// sin(O1*t)/cos(O1*t) <   ( R1*(1-cos(O1*t)) - R2*(1-cos(O2*t)) ) /  ( L0 + R2*sin(O2*t) - R1*sin(O1*t) )




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
