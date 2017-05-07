#ifndef  SpaceLaunchODE_h
#define  SpaceLaunchODE_h

// TODO : in futere maybe warp it into obejct ??

constexpr static int nCP  = 10;
double  thetaCPs [nCP];
double  thrustCPs[nCP];
Vec3d   dirCPs   [nCP];

double  uT    = 50.0; // duration of each CP interval
double inv_uT = 1.0/uT;

// rocket parameters - Saturn V
double  vexhaust      = 4500.0; // [km/s]
double  dm_F          = 1.0/vexhaust;
double  mass_initial  = 2.970e+6; // put there zero to switch off the engine
double  mass_empty    = 130e+3 + 140e+3 + 40e+3 + 13.5e+3;
double  thrust_full   = 35e+6;
double  AeroArea      = 320.0; // diameter 10.1 m

// planet parameters
const double planet_R   = 6371e+3;
double       planet_GM  = 3.9860044189e+14; //https://en.wikipedia.org/wiki/Standard_gravitational_parameter
Vec3d        planet_pos = {0.0,-planet_R,0.0};

double plante_rho0  = 1.34;    // [kg/m3]
double planet_zrate = 0.14e-3; // [1/m]

int nrho;
int nCds;
double inv_dvM = 1/0.25;      // [1/Mach]
double inv_dh  = 1/2000.0;  // [1/m]
double hmax    = 84000.0d;  // [m]
double vMmax   =     9.0d;  // [Mach]
double * planet_rho_CPs = NULL;
double * rocket_Cd_CPs  = NULL;

// ---- log

void (*logFunc)();
double bLog         = false;
double log_tmax     = 1000.0;
double log_dt_trig  = 1.0;
double log_t_trig   = 0.0;

double log_t;
double log_Cd;
double log_rho;
double log_FD;
double log_m;
double log_h;
double log_v;
double log_vx;
double log_vy;
double log_T;
double log_Tx;
double log_Ty;
double log_G;

// ==== Globals

double view_scale = 1e-6;

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
    spline_sample( h, inv_dh, u, icp );
    //printf( "%i %g %g getAtmosphereDensity\n", icp, u, h );
    return getSpline( u, planet_rho_CPs+icp );
}

inline double getDragCoef( double M ){
    // https://www.orekit.org/
    // http://www.braeunig.us/apollo/saturnV.htm
    // https://en.wikipedia.org/wiki/Normal_shock_tables
    // http://www.aerorocket.com/AeroCFD/Instructions.html
    // https://www.grc.nasa.gov/www/k-12/airplane/normal.html
    // http://what-when-how.com/space-science-and-technology/rocket-propulsion-theory/
    if(M<vMmax){
        double u;
        int    icp;
        spline_sample( M, inv_dvM, u, icp );
        //printf( "%i %g %g getDragCoef\n", icp, u, M );
        return getSpline( u, rocket_Cd_CPs+icp );
    }else{
        return rocket_Cd_CPs[nCds-1];
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
        double DD    = 0.5*AeroArea*rho*Cd*v*v;

        //printf( " rho %g  Cd %g \n", rho, Cd );

        if(bLog){
            log_Cd    = Cd;
            log_rho   = rho;
            log_FD    = DD;
            //log_v     = v;
            //log_h     = h;
        }
        return vel*( -DD/v );
    }else{
        return {0.0,0.0,0.0};
    }
}

void getODEDerivs( double t, int n, double * Ys, double * dYs ){

    if( ( t>log_t_trig )&( t<log_tmax ) ) bLog = true;

    Vec3d pos = *(Vec3d*) Ys;
    Vec3d vel = *(Vec3d*)(Ys+3);
    double m  = Ys[6];

    // Vec3d acc; acc.set(0.0d);
    // double r = addGravity( pos-planet_pos, acc, planet_GM );

    // --- Gravity
    Vec3d dp; dp.set_sub(pos, planet_pos );
    double r2  = dp.norm2();
    double r   = sqrt(r2);
    dp.mul(1/r);
    double G = -planet_GM/r2;
    Vec3d acc; acc.set_mul( dp, G );

    double h = r - planet_R;
    if( h < hmax ){
        acc.add_mul( getAirDragForce( vel, h ), 1/m );
        //exit(0);
    }else{
        if(bLog){ log_Cd= 0; log_rho= 0; log_FD= 0; }
    }

    if(bLog){
        double cr   = dp.dot(vel);
        Vec3d velT  = vel - dp*cr;
        double cl   = velT.norm();
        log_v  = vel.norm();
        log_vx = cl;
        log_vy = cr;
        log_G  = G;
        log_h  = h;
        log_m  = m;
    }

    // engine thrust
    double dm = 0.0d;
    if( m > mass_empty ){   //  check if propelant exhausted
        // --- interpolate from inputs
        //double tcp = t*inv_uT;
        //int    icp = (int) tcp;
        int icp; double u;
        spline_sample( t, inv_uT, u, icp );
        if((icp+3)<nCP){
            //double u   = tcp - icp;

            // TODO : this can be very much optimized if we evaluate spline basis instead !!!!
            Vec3d  dir    = getSpline3d( u, dirCPs   +icp );
            double thrust = getSpline  ( u, thrustCPs+icp );

            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawVecInPos( dir*5.0, pos*view_scale );
            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( pos*view_scale, pos*view_scale + dir );
            //Vec3d dir_ = {1.0,dir.y,0.0};
            //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( pos*view_scale, pos*view_scale + dir_ );
            //printf( "t %g tcp %g icp %i u %g (%g,%g,%g) \n", t, tcp, icp, u, dir.x,dir.y,dir.z  );

            dm =     -dm_F * thrust; // propelant consumption

            if(bLog){
                double cr   = dp.dot(dir);
                Vec3d dirT  = dir - dp*cr;
                double cl   = dirT.norm();
                log_Tx      = thrust*cl;
                log_Ty      = thrust*cr;
                log_T       = thrust;
            }

            //printf( "T %g T/m %g acc %g m %g dm %g \n", thrust,thrust/m, acc.norm(), m, dm  );
            acc.add_mul( dir, thrust/m );
        }else{
            if(bLog){ log_Tx= 0; log_Ty= 0; log_T= 0; }
        }
        //exit(0);
    };

    if(bLog){
        log_t = t;
        logFunc();
        bLog  = false;
        log_t_trig += log_dt_trig;
    }

    // TODO : air drag

    ((Vec3d*) dYs   )->set( vel );
    ((Vec3d*)(dYs+3))->set( acc );
    dYs[6] = dm;
};





#endif
