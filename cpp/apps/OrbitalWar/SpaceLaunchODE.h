#ifndef  SpaceLaunchODE_h
#define  SpaceLaunchODE_h

// TODO : in futere maybe warp it into obejct ??

namespace SpaceLaunchODE{

typedef struct{
    int n    = 10;
    double * thetaCPs ;
    double * thrustCPs;
    Vec3d  * dirCPs   ;
    double uT     = 50.0; // duration of each CP interval
    double inv_uT = 1.0/uT;
} LaunchSpline;

typedef struct{
    int n;
    double dvM;
    double inv_dvM = 1/0.25;      // [1/Mach]
    double vMmax   =     9.0d;  // [Mach]
    double * Cd_CPs  = NULL;
} AeroSpline;

typedef struct{
    int n;
    double dh;
    double inv_dh  = 1/2000.0;  // [1/m]
    double hmax    = 84000.0d;  // [m]
    double * rho_CPs = NULL;
    //double rho0  = 1.34;    // [kg/m3]
    //double zrate = 0.14e-3; // [1/m]
} AtmosphereSpline;

// rewrite to c
typedef struct{
    // rocket parameters - Saturn V
    double  vexhaust      = 4500.0; // [km/s]
    double  dm_F          = 1.0/vexhaust;
    double  mass_initial  = 2.970e+6; // put there zero to switch off the engine
    double  mass_empty    = 130e+3 + 140e+3 + 40e+3 + 13.5e+3;
    double  thrust_full   = 35e+6;
    double  AeroArea      = 320.0; // diameter 10.1 m
} Rocket;

typedef struct{
    // planet parameters
    double R     = 6371e+3;
    double GM    = 3.9860044189e+14; //https://en.wikipedia.org/wiki/Standard_gravitational_parameter
    Vec3d  pos   = (Vec3d){0.0,-6371e+3,0.0};
} Planet;

typedef struct{
    double on       = false;
    double tmax     = 1000.0;
    double dt_trig  = 1.0;
    double t_trig   = 0.0;
} LogTrig;

typedef struct{
    double t;
    double Cd;
    double rho;
    double FD;
    double m;
    double h;
    double v;
    double vx;
    double vy;
    double T;
    double Tx;
    double Ty;
    double G;
} LogVars;

// ==== Globals

LaunchSpline      launch;
AeroSpline        aero;
AtmosphereSpline  atmosphere;

Rocket rocket;
Planet planet;

LogTrig logTrig;
LogVars logVars;
void (*logFunc)();

// ==== Functions

inline void spline_sample( double t, double inv_dt, double& u, int& icp ){
    double tcp = t*inv_dt;
    icp = (int) tcp;
    u   = tcp - icp;
}

inline double getSpline( double u, double* CPs ){
    double cp0 = CPs[0];
    double cp1 = CPs[1];
    double cp2 = CPs[2];
    double cp3 = CPs[3];
    return Spline_Hermite::val( u, cp1, cp2, (cp2-cp0)*0.5d, (cp3-cp1)*0.5d );
}

inline Vec3d getSpline3d( double u, Vec3d* CPs ){
    //double b0,b1,b2,b3;
    //basis( TYPE x, TYPE& c0, TYPE& c1, TYPE& d0, TYPE& d1 );
    Vec3d cp0 = CPs[0];
    Vec3d cp1 = CPs[1];
    Vec3d cp2 = CPs[2];
    Vec3d cp3 = CPs[3];
    return { // TODO Optimize this
        Spline_Hermite::val( u, cp1.x, cp2.x, (cp2.x-cp0.x)*0.5d, (cp3.x-cp1.x)*0.5d ),
        Spline_Hermite::val( u, cp1.y, cp2.y, (cp2.y-cp0.y)*0.5d, (cp3.y-cp1.y)*0.5d ),
        Spline_Hermite::val( u, cp1.z, cp2.z, (cp2.z-cp0.z)*0.5d, (cp3.z-cp1.z)*0.5d )
    };
}

inline double addGravity( const Vec3d& dp, Vec3d& G, double GM ){
    double r2  = dp.norm2();
    double r   = sqrt(r2);
    //double ir2 = 1/r2;
    //double ir  = sqrt( ir2 );
    //double ir3 = 1ir2*ir;
    G.add_mul( dp, -GM/(r*r2));
    return r;
}

inline double getAtmosphereDensity( double h ){
    // http://farside.ph.utexas.edu/teaching/sm1/lectures/node56.html
    // http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
    // return plante_rho0 * exp( planet_zrate * h ); // TODO : read from spline
    double u;
    int    icp;
    spline_sample( h, atmosphere.inv_dh, u, icp );
    //printf( "%i %g %g getAtmosphereDensity\n", icp, u, h );
    return getSpline( u, atmosphere.rho_CPs+icp );
}

inline double getDragCoef( double M ){
    // https://www.orekit.org/
    // http://www.braeunig.us/apollo/saturnV.htm
    // https://en.wikipedia.org/wiki/Normal_shock_tables
    // http://www.aerorocket.com/AeroCFD/Instructions.html
    // https://www.grc.nasa.gov/www/k-12/airplane/normal.html
    // http://what-when-how.com/space-science-and-technology/rocket-propulsion-theory/
    if(M<aero.vMmax){
        double u;
        int    icp;
        spline_sample( M, aero.inv_dvM, u, icp );
        //printf( "%i %g %g getDragCoef\n", icp, u, M );
        return getSpline( u, aero.Cd_CPs+icp );
    }else{
        return aero.Cd_CPs[aero.n-1];
    }
}

inline Vec3d getAirDragForce( const Vec3d& vel, double h ){
    double v2    = vel.norm2();
    //printf( " v2 %g \n", v2 );
    if(v2>1e-16){
        double v     = sqrt(v2);
        // http://aero.stanford.edu/stdatm.html
        double rho   = getAtmosphereDensity( h );
        double vMach = v/340.0;
        double Cd    = getDragCoef(vMach);
        double DD    = 0.5*rocket.AeroArea*rho*Cd*v*v;

        //printf( " rho %g  Cd %g \n", rho, Cd );

        if(logTrig.on){
            logVars.Cd    = Cd;
            logVars.rho   = rho;
            logVars.FD    = DD;
            //log_v     = v;
            //log_h     = h;
        }
        return vel*( -DD/v );
    }else{
        return {0.0,0.0,0.0};
    }
}

void getODEDerivs( double t, int n, double * Ys, double * dYs ){

    if( ( t>logTrig.t_trig )&( t<logTrig.tmax ) ) logTrig.on = true;

    Vec3d pos = *(Vec3d*) Ys;
    Vec3d vel = *(Vec3d*)(Ys+3);
    double m  = Ys[6];

    // Vec3d acc; acc.set(0.0d);
    // double r = addGravity( pos-planet_pos, acc, planet_GM );

    // --- Gravity
    Vec3d dp; dp.set_sub(pos, planet.pos );
    double r2  = dp.norm2();
    double r   = sqrt(r2);
    dp.mul(1/r);
    double G = -planet.GM/r2;
    Vec3d acc; acc.set_mul( dp, G );

    double h = r - planet.R;
    if( h < atmosphere.hmax ){
        acc.add_mul( getAirDragForce( vel, h ), 1/m );
        //exit(0);
    }else{
        if(logTrig.on){ logVars.Cd= 0; logVars.rho= 0; logVars.FD= 0; }
    }

    if(logTrig.on){
        double cr   = dp.dot(vel);
        Vec3d velT  = vel - dp*cr;
        double cl   = velT.norm();
        logVars.v  = vel.norm();
        logVars.vx = cl;
        logVars.vy = cr;
        logVars.G  = G;
        logVars.h  = h;
        logVars.m  = m;
    }

    // engine thrust
    double dm = 0.0d;
    if( m > rocket.mass_empty ){   //  check if propelant exhausted
        // --- interpolate from inputs
        //double tcp = t*inv_uT;
        //int    icp = (int) tcp;
        int icp; double u;
        spline_sample( t, launch.inv_uT, u, icp );
        if((icp+3)<launch.n){
            //double u   = tcp - icp;

            // TODO : this can be very much optimized if we evaluate spline basis instead !!!!
            Vec3d  dir    = getSpline3d( u, launch.dirCPs   +icp );
            double thrust = getSpline  ( u, launch.thrustCPs+icp );

            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawVecInPos( dir*5.0, pos*view_scale );
            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( pos*view_scale, pos*view_scale + dir );
            //Vec3d dir_ = {1.0,dir.y,0.0};
            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( pos*view_scale, pos*view_scale + dir_ );
            //printf( "t %g tcp %g icp %i u %g (%g,%g,%g) \n", t, tcp, icp, u, dir.x,dir.y,dir.z  );

            dm =     -rocket.dm_F * thrust; // propelant consumption

            if(logTrig.on){
                double cr   = dp.dot(dir);
                Vec3d dirT  = dir - dp*cr;
                double cl   = dirT.norm();
                logVars.Tx      = thrust*cl;
                logVars.Ty      = thrust*cr;
                logVars.T       = thrust;
            }

            //printf( "T %g T/m %g acc %g m %g dm %g \n", thrust,thrust/m, acc.norm(), m, dm  );
            acc.add_mul( dir, thrust/m );
        }else{
            if(logTrig.on){ logVars.Tx= 0; logVars.Ty= 0; logVars.T= 0; }
        }
        //exit(0);
    };

    if(logTrig.on){
        logVars.t = t;
        logFunc();
        logTrig.on      = false;
        logTrig.t_trig += logTrig.dt_trig;
    }

    // TODO : air drag

    ((Vec3d*) dYs   )->set( vel );
    ((Vec3d*)(dYs+3))->set( acc );
    dYs[6] = dm;
};

}; // namespace SpaceLaunchODE

#endif
