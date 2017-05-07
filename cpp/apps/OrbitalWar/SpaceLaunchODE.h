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

// ==== Globals

double view_scale = 1e-6;

// ==== Functions

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

inline void addGravity( const Vec3d& dp, Vec3d& G, double GM ){
    double ir2 = 1.0d/dp.norm2();
    double ir  = sqrt( ir2 );
    double ir3 = ir2*ir;
    G.add_mul( dp, -ir3*GM);
}


inline double getAtmosphereDensity( double h ){
    // http://farside.ph.utexas.edu/teaching/sm1/lectures/node56.html
    // http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
}

inline double getDragCoef( double M ){
    // http://www.braeunig.us/apollo/saturnV.htm
}

inline void addDrag( const Vec3d& vel, Vec3d& D, double h ){
    double v2    = vel.norm2();

    double v     = sqrt(v2);
    // http://aero.stanford.edu/stdatm.html
    double rho   = getAtmosphereDensity( h );
    double vMach = v/340.0;
    double Cd    = getDragCoef(vMach);
    double DD    = 0.5*AeroArea*rho*Cd*v*v;

    D.add_mul( vel, -DD/v );

}

inline void spline_sampe( double t, double& u, int& icp ){
    double tcp = t*inv_uT;
    icp = (int) tcp;
    u   = tcp - icp;
}

void getODEDerivs( double t, int n, double * Ys, double * dYs ){
    Vec3d pos = *(Vec3d*) Ys;
    Vec3d vel = *(Vec3d*)(Ys+3);
    double m  = Ys[6];

    Vec3d acc; acc.set(0.0d);
    addGravity( pos-planet_pos, acc, planet_GM );

    // engine thrust
    double dm = 0.0d;
    if( m > mass_empty ){   //  check if propelant exhausted
        // --- interpolate from inputs
        //double tcp = t*inv_uT;
        //int    icp = (int) tcp;
        int icp; double u;
        spline_sampe( t, u, icp );
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

            //printf( "T %g T/m %g acc %g m %g dm %g \n", thrust,thrust/m, acc.norm(), m, dm  );
            acc.add_mul( dir, thrust/m );
        }
        //exit(0);
    };

    // TODO : air drag

    ((Vec3d*) dYs   )->set( vel );
    ((Vec3d*)(dYs+3))->set( acc );
    dYs[6] = dm;
};





#endif
