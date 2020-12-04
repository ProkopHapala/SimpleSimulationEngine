
#ifndef  spaceCombat_h
#define  spaceCombat_h

#include <math.h>
//#include "fastMath.h"
#include "appliedPhysics.h"

//#include "Combat.h"

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


// =========== Classes


struct whippleShieldType{
    int n;
    double layerDens;  // [kg/m^2]
    double spacing;    // [m]
    double critEdens;  // [J/m^2]

    void impact( double& m, double& R, double& v, double& vT ){
        for(int i=0; i<n; i++){
            // --- scatter
            double Ek    = m*( 0.5*v*v + vT*vT*0.2 );
            double area  = (M_PI*R*R*0.25);
            double Ecrit = area*critEdens;
            printf( "impact[%i] Ek %g[MJ] Ecrit %g[MJ] | area %g[m^2] R %g[m] \n", i, Ek*1e-6, Ecrit*1e-6, area, R );
            if(  Ek < Ecrit ){
                v=0; vT=0; break;
            }else{
                double efac = sqrt( 1-Ecrit/Ek );
                v *=efac;
                vT*=efac;
            };
            double M = layerDens*area;
            // v1*m = (m+M)*v2
            double m_   = m+M;
            printf( "impact[%i] m %g[kg] M %g[kg] m/(m+M) %g \n", i, m, M, m/m_ );
            double invm_ = 1./m_;
            double v_ = v*m*invm_;
            double Eside = Ek - 0.5*v_*v_*m_;
            double vT_   = sqrt( 4*Eside*invm_ );
            printf( "impact[%i] Ek %g[MJ] Eside %g[MJ] vT_ %g[km/s] v_ %g[km/s] m_ %g[kg]\n", i, Ek*1e-6, Eside*1e-6, vT_*1e-3, v_*1e-3, m_ );
            vT=vT_;
            v=v_;
            m=m_;
            // --- fly
            R += spacing*(vT/v);
        }
    }

    void fromString(const char* s){
        sscanf( s,                     "%i %lf %lf %lf", &n, &layerDens, &spacing, &critEdens );
        printf(    "whippleShieldType : %i %g %g %g \n", n, layerDens, spacing, critEdens );
    }
};


struct ProjectileType{

    double mass;    // [kg]
    double caliber; // [m]

    double getAreaDensity(){ return mass/(M_PI*caliber*caliber*0.25); } // [kg/m^2]

    void fromString(const char* s){
        sscanf( s,                 "%lf %lf"  , &mass, &caliber );
        printf(    "ProjectileType : %g %g \n", mass, caliber );
    }

};

struct SpaceGunType{
    int    kind;      // 0=laser, 1=railgun
    double length;    // [m]
    double maxForce;  // [N]
    double maxPower;  // [W]
    double scatter;   // tg(a)   - angular uncertainity of flight direction
    double fireRate;  // [Hz,1/s] // --- ToDo : fireRate should be limited by delay time (time spend by projectile in the barrel)

    void fromString(const char* s){
        sscanf( s,                "%lf %lf %lf %lf %lf"  , &length, &maxForce, &maxPower, &scatter, &fireRate  );
        printf(    "SpaceGunType : %g  %g  %g  %g  %g \n", length,  maxForce,   maxPower,  scatter,  fireRate );
    }

    double getMuzzleVelocity( ProjectileType* shotType, double& t ){
        printf( "getMuzzleVelocity \n" );
        // s = 0.5*a*t^2
        // v = a*t
        // E = F*s
        double m = shotType->mass;
        // --- part #1 - force limited acceleration
        double a1 = maxForce/m;
        double v1 = maxPower/maxForce;
        double t1 = v1/a1;
        double s1 = 0.5*a1*t1*t1;
        printf( "getMuzzleVelocity.#1 a1 %g[m/s^2] v1 %g[m/s] t1 %g[s] s1 %g[m] | m %g[kg]  maxPower %g[W] maxForce %g[N] \n", a1, v1, t1, s1,    m,  maxPower, maxForce );
        if(s1>length){
            t = sqrt(2*length/a1);
            printf( "getMuzzleVelocity: length(%g)<s1(%g) => t %g[%s] v %g[km/s] \n", length, s1, t, t*a1*1e-3  );
            return t*a1;
        }
        double s2 = length - s1;
        // --- part #2 - power limited acceleration
        // P = F*v  = m*(dv/dt)*v
        // (P/m)*dt = v*dv
        //  2*t*P/m = v^2 - v0^2
        //   v = sqrt( v0^2 +  2*t*P/m )
        //   v = v0 * sqrt( 1 + (P/(0.5*m*v0^2))*t  )
        //   v = v0 * sqrt( 1 + K*t  )   .... K = P/(0.5*m*v0^2)
        //   s/v0 = interal   (1+K*t)^0.5   dt
        //   s/v0 = (2/3)*(1+K*t)^(3/2) + C
        //   0    = (2/(3K)) + C       .... C = -2/(3K)
        //    (s/v0 - 2/3)^2 = (1+K*t)^3
        //   t  = (-K^2 v^2 + (K^6 v^4 (s - C v)^2)^(1/3))/(K^3 v^2)
        double Ek1  = 0.5*m*v1*v1;
        double invK = Ek1/maxPower;
        double C    = (2./3.)*invK;
        // double v2 = v0 * np.sqrt( 1. + K*t )
        //double  s2 = v0*(   C*(1.+K*t)**(3./2.) - C  );   // Distance
        double vfac  = pow( 1. + s2/( C * v1 ), 1./3. );
        double t2  = ( vfac*vfac - 1. )*invK;
        printf( "getMuzzleVelocity.#2 Ek1 %g[MJ] K %g vfac %g t2 %g[s] \n", Ek1*1e-6, 1/invK, vfac, t2 );
        t = t1+t2;
        double v = v1 * vfac; // DEBUG
        printf( "getMuzzleVelocity.end v %g[km/s] t %g[s] \n", v*1e-3, t );
        return v;
    }

    //ProjectileType* shotType;
};

struct SpaceGun{
    int n;         // number of guns
    int nDisabled;

    ProjectileType* shotType;
    SpaceGunType*   gunType;

    //double vMuzzle;     // muzzle velocity of projectile
    //double muzzleDelay; // time between trigger and projectile releasing gun

    SpaceGun()=default;
    SpaceGun(int n_, SpaceGunType* gunType_, ProjectileType* shotType_ ){ n=n_; gunType=gunType_; shotType=shotType_; }

};

struct SpaceSalvo{

    int n;
    //Vec3d pos,vel;
    ProjectileType* shotType;
    // derived properties
    double speed;     // [m/s],   299792458 m/s for laser
    //double Energy;
    double scatter0;  //
    double delay;

    inline double timeToTarget(double dist){ return dist/speed; }
    double getSpread(double dist, double aDelta){
        double tof = timeToTarget( dist );
        double t   = delay + tof;
        double scatter = scatter0*dist  +  aDelta*t*t*0.5;
        printf( "getSpread : t %g(tof %g, t0 %g)[s] scatter %g(d %g ,d^2 %g)[m] \n",  t, tof, delay,      scatter, scatter0*dist, aDelta*t*t*0.5 );
        return scatter;
    }

    inline void fromGun( SpaceGun& guns ){
        n        = guns.n;
        scatter0 = guns.gunType->scatter;
        speed    = guns.gunType->getMuzzleVelocity( guns.shotType, delay );
        shotType = guns.shotType;
    }
    SpaceSalvo(SpaceGun& guns){ fromGun(guns); };
    SpaceSalvo()=default;
};

struct OpticalMaterial{
    int n;
    double albedol;
};

struct ProjectedTarget{  // target projected in particular direction
    double area;
    double areaMass;           // for gause-gun  [kg/m^2]
    OpticalMaterial*   mat;  // for  laser     [in various wavelenghths]
    whippleShieldType* wshield;

    double HPs;        // [Joule] hitpoints
    double HPexponent;

    double health; // perentage of function
    double damage;

    void fromString(const char* s){
        sscanf( s,                   "%lf %lf %lf"  , &area, &HPs, &HPexponent    );
        printf(    "projectedTarget : %g  %g  %g \n",  area,  HPs,  HPexponent );
    }

    void reset(){ damage=0; };

    void getDamage( double E_Damage ){
        // ToDo : consider co-locality of hits, not just damage
        damage += E_Damage;
        health  = pow( fmax( 0., 1.-(damage/HPs) ), HPexponent );
        printf( "getDamage   E_Damage %g[MJ] damage %g[MJ] HPs %g[MJ] damage/HPs %g health %g HPexponent %g \n", E_Damage*1e-6, damage*1e-6, HPs*1e-6,   damage/HPs, health, HPexponent  );
    };
};

struct CombatAssembly{
    Vec3d dir,pMe,pHi;
    //std::vector<SpaceGun> guns;           // visible along certain direction

    std::vector<SpaceSalvo>      salvos;
    std::vector<ProjectedTarget> targets; // targets projected in direction, ordered by distance

    void fireGun( SpaceGun& guns ){
        SpaceSalvo s(guns);
        salvos.push_back( s );
    }

    void addTarget( ProjectedTarget& t  ){ targets.push_back( t ); }

    void reset(){ for( ProjectedTarget& t : targets){ t.reset(); }; };
    void colide( double dist, double aDelta ){
        int DEBUG_i = 0;
        for( SpaceSalvo& s : salvos){
            printf( "colide s[%i] n %i delay %g[s] \n", DEBUG_i, s.n, s.delay );
            double R0         = s.getSpread( dist, aDelta );
            double spreadArea = M_PI*R0*R0;
            printf( "colide spread %g[m] area %g[m^2] \n", R0, spreadArea );
            //double cumHitProb =   // TODO : occlusion
            for( ProjectedTarget& t : targets){
                double hitProb = t.area/(spreadArea + t.area);
                // TODO : there should be perhaps some randomness
                printf( "hitProb %g | t.area %g[m^2] spreadArea %g[m^2] R %g[m]\n", hitProb, t.area, spreadArea, R0 );
                double R  = s.shotType->caliber;
                double v  = s.speed;
                double m  = s.shotType->mass;
                double vT = 0;
                t.wshield->impact( m, R, v, vT );
                double E  = (vT*vT*0.25 + v*v*0.5)*m;
                printf( "E %g[MJ] v %g[km/s] vT %g[km/s] m %g[kg] \n",  E*1e-6,  v*1e-3, vT*1e-3, m );
                if(E>1e-8){
                    t.getDamage( E );
                }
            }
            DEBUG_i++;
        }
    }

};






#endif
