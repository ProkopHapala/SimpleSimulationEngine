
#ifndef  SpaceBodies_h
#define  SpaceBodies_h

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include <string>
#include "Body.h"
#include "ODEintegrator.h"

#include "spline_hermite.h"

#include "appliedPhysics.h"


constexpr const double DEG_2_RAD = M_PI/180;


/*
class SpaceCraftBody : public RigidBody{ public:
}
*/

inline double findNextInd( double t, int imax, const double* ts, int& i ){
    double tleft;
    for(;i<imax;i++){
        //printf("%i %f %f \n", i, t, ts[i] );
        if(t<ts[i])break;
    };
    double ot =  ts[i-1];
    double dt = (ts[i]-ot);
    //printf( "findNextInd i=%i t=%f ot=%f ts=%f dt=%f u=%f \n", i, t, ot, ts[i], dt, (t-ot)/dt );
    return (t-ot)/dt;
}

void nonUni2spline( double t0, double dt, int n, const double* ts, const Vec3d* ps, int nout, Vec3d* out ){
    int j=1;
    for(int i=0; i<nout; i++){
        double t = dt*i + t0;
        double u = findNextInd( t, n-1, ts, j);
        out[i] = ps[j-1]*(1-u) + ps[j]*u;
        //printf( "%i %i %f %f (%f,%f,%f) \n", i, j, t, u, out[i].x, out[i].y, out[i].z );
    };
    //exit(0);
}

double julian_date( double day, double month, double year ){
    //double m,y,a,b,c,d,jd;
    double m = month;
    double y = year;
    if ( m < 3 ){
        y=y-1;
        m=m+12;
    }
    double a  = floor    ( y * 0.01 );
    double b  = 2-a+floor( a * 0.25 );
    double c  = floor(365.25*y);
    double d  = floor(30.6001*(m+1));
    double jd = b+c+d+day+1720994.5;
    return jd;
}

double string2julian_date( char* s ){
    double d = atof(s+6);   s[6]=0;
    double m = atof(s+4);   s[4]=0;
    double y = atof(s);
    return julian_date(d,m,y);
}

double true_anomaly( double ma, double ec, double errConv ){
    // https://en.wikipedia.org/wiki/True_anomaly
    // fast approximation https://en.wikipedia.org/wiki/True_anomaly#From_the_mean_anomaly
    // https://en.wikipedia.org/wiki/Equation_of_the_center#Series_expansion
    //double errConv = 1e-6;
    double e    = ma;
    double err = 1.0;
    while( err > errConv ){
        e   = ma + ec*sin(e);
        err = abs( e - ec*sin(e) - ma );
    }
    return 2*atan( sqrt( (1+ec)/(1-ec) ) * tan(e*0.5) );
}

double semi_major2period( double semi_major ){
    return  365.25 * sqrt( semi_major*semi_major*semi_major );
}


double mean_anomaly_at_epoch( double period, double mean_anomaly, double epoch_diff ){
    //mean_anomaly = mean_anomaly + (new_epoch-epoch)*2*M_PI/period;
    double  ma   = mean_anomaly/(2*M_PI) + epoch_diff/period;
    double ima;
    return modf(ma,&ima)*2*M_PI;
}

struct OrbitalElements{
    //int id; // identify to which body it belongs
    // dimensions
    double semi_major;
    double eccentricity;
    // orbit rotation angles (euler?)
    double inclination;
    double arg_periapsis;
    double node_longitude;
    // position on orbit and time
    double mean_anomaly=0;
    double period;
    double epoch=0;
    //double true_anomaly;

    OrbitalElements() = default;
    OrbitalElements(double semi_major, double eccentricity, double inclination, double node_longitude, double arg_periapsis, double mean_anomaly ):
    semi_major(semi_major),eccentricity(eccentricity),inclination(inclination),node_longitude(node_longitude),arg_periapsis(arg_periapsis),mean_anomaly(mean_anomaly){ period=semi_major2period(semi_major); };

    void shift_epoch( double new_epoch ){
        mean_anomaly = mean_anomaly_at_epoch( period, mean_anomaly, (new_epoch-epoch) );
        epoch = new_epoch;
    }

};


struct Orbit{
    //int id; // identify to which body it belongs
    //Vec3d apsis;
    //Vec3d node;
    Mat3d rot;
    double eccentricity;
    double semi_major;
    // position on orbit and time
    double mean_anomaly;
    double period;
    double epoch;

    Orbit() = default;
    Orbit(const OrbitalElements& el){fromElements(el);}

    inline void shift_epoch( double new_epoch ){
        mean_anomaly = mean_anomaly_at_epoch( period, mean_anomaly, (new_epoch-epoch) );
        epoch = new_epoch;
    }

    inline double true_anomaly_at_epoch(double at_epoch)const{
        double mean_anomaly_epoch = mean_anomaly_at_epoch( period, mean_anomaly, at_epoch-epoch );
        return true_anomaly( mean_anomaly_epoch, eccentricity, 1e-6 );
    }

    void fromElements( const OrbitalElements& el ){
        //Mat3d rot;
        // https://en.wikipedia.org/wiki/Orbital_elements#Euler_angle_transformations
        //rot.fromEuler( el.node_longitude, el.inclination+M_PI_2, el.arg_periapsis );
        rot.fromEuler_orb( el.inclination, el.node_longitude, el.arg_periapsis );
        semi_major   = el.semi_major;
        eccentricity = el.eccentricity;

        mean_anomaly = el.mean_anomaly;
        period       = el.period;
        epoch        = el.epoch;
    }

    Vec3d pointAtEpoch(double at_epoch)const{
        double anomaly =  true_anomaly_at_epoch(at_epoch);
        double ca   = cos(anomaly);
        double sa   = sin(anomaly);
        double slr  = semi_major*( 1 - sq(eccentricity) );   // semi-latus rectum   https://en.wikipedia.org/wiki/Ellipse#Polar_form_relative_to_focus
        double r    = slr/(1+eccentricity*ca);
        return rot.a*(r*ca) + rot.b*(r*sa);
    }

    void toPoints( double ang0, double dang, int n, Vec3d* ps )const{
        double slr  = semi_major*( 1 - sq(eccentricity) );   // semi-latus rectum   https://en.wikipedia.org/wiki/Ellipse#Polar_form_relative_to_focus
        double ang = ang0;
        for(int i=0; i<n; i++){
            double ca   = cos(ang);
            double sa   = sin(ang);
            double r    = slr/(1+eccentricity*ca);
            ps[i] = rot.a*(r*ca) + rot.b*(r*sa);
            //printf( "%i %g %g %g \n", i, ps[i].x, ps[i].y, ps[i].z, ang );
            ang+=dang;
        }
    }

    void toPoints_epochs( double epoch0, double epoch1, int n, Vec3d* ps )const{
        double ang0 =  true_anomaly_at_epoch(epoch0);
        double ang1 =  true_anomaly_at_epoch(epoch1);
        double dang = (ang1-ang0)/(n-1);
        toPoints( ang0, dang, n, ps );
    }

};


//class SpaceBody : public PointBody  { public:
class SpaceBody : public RigidBody  { public:

    std::string name;
    char Stype[16];
    double radius;
    Vec3d sizes = (Vec3d){1.0,0.5,0.25};

    Vec3d * trjPos    = 0; // positions in some times
    //Vec3d * trjVel    = 0; // velocities in some times
    Vec3d * trjThrust = 0; // vector of thrust in time

    Orbit*     orbit    =0;
    SpaceBody* orbCenter=0;

    SpaceBody()=default;
    SpaceBody( std::string name, SpaceBody* orbCenter, Orbit* orbit, double radius ):name(name),orbCenter(orbCenter),orbit(orbit),radius(radius){};
    SpaceBody( std::string name, SpaceBody* orbCenter, double a, double e, double inc, double lan, double apa, double ma ):name(name),orbCenter(orbCenter){
        orbit = new Orbit( OrbitalElements(a,e,inc,lan, apa, ma) );
    };

    Vec3d pointAtEpoch(double epoch)const{
        if(orbit){
            return orbit->pointAtEpoch( epoch );
        }else if (trjPos){
            return (Vec3d){NAN,NAN,NAN};
            // ToDo
        }
        return (Vec3d){NAN,NAN,NAN};
    }

    inline Vec3d getTrjPos( int iTrj, double du )const {
        Vec3d p; Spline_Hermite::curve_point( du,trjPos[iTrj-1],trjPos[iTrj],trjPos[iTrj+1],trjPos[iTrj+2], p); return p;
    };

    inline Vec3d getThrust(int itrj, double du ){
        //if( trjThrust ){
            //printf( "%i %f   (%f,%f,%f)   (%f,%f,%f) \n", itrj, du, trjThrust[itrj].x, trjThrust[itrj].y, trjThrust[itrj].z,  trjThrust[itrj+1].x, trjThrust[itrj+1].y, trjThrust[itrj+1].z );
            return trjThrust[itrj]*(1-du) + trjThrust[itrj+1]*du;
        //}else{
        //    return (Vec3d){0.0,0.0,0.0};
        //}
    }
    //Vec3d getThrust(double t){};


    void fromString_astorb(char* s){
        // from ftp://ftp.lowell.edu/pub/elgb/astorb.html
        //double lo,so;file:///home/prokop/git/SimpleSimulationEngine/cpp/apps/OrbitalWar/data/asteroids_small_legend.dat

        OrbitalElements el;
        s[26]  = 0; // behind name
        s[65]  = 0; // behind diameter
        s[183] = 0; // behind semimajor-axis

        name   = s+7;
        radius = atof( s+59 ) * 0.5 * 1000.0;  // [m]

        double density = 2000.0;                         // [kg/m^3]
        double volume = radius*radius*radius*4./3.*M_PI; // [m^3]
        mass = volume * density;                         // [kg]

        //sscanf(s+65,"%s", &Stype ); // spectral type
        //strncpy( Stype, s+66, 16 );
        char c;
        //for(int i=0;i<16;i++){ Stype[i]=s[66+i]; }; Stype[15]=0;
        //c = s[65];
        printf( "%s \n", s );
        printf( "Stype '%c' '%c' '%c' '%c' '%c' '%c' '%c' '%c' '%c' \n", s[60],s[61],s[62],s[63],s[64],s[65],s[66], s[67], s[68]  );
        //printf( "Stype '%c' '%c' '%c' '%c' '%c' '%c' '%c' '%c' '%c' \n", s[60],s[61],s[62],s[63],s[64],s[65],s[66]s[67] );
        //printf( "Stype `%s` \n", Stype );
        //strncpy( Stype, s+66, 16 );
        //exit(0);

        char sdate[16];
        // 13 Mean anomaly, deg.
        // 14 Argument of perihelion, deg (J2000.0).
        // 15 Longitude of ascending node, deg (J2000.0).
        // 16 Inclination, deg (J2000.0).
        // 17 Eccentricity.
        // 18 Semimajor axis, AU.

        double mean_anomaly;
        sscanf( s+106, "%s %lf %lf %lf %lf %lf %lf", sdate, &el.mean_anomaly, &el.arg_periapsis, &el.node_longitude, &el.inclination, &el.eccentricity, &el.semi_major );
        //printf( "%s R %g date{%s} ma %g so %g lo %g inc %g e %g a %g \n", name.c_str(), radius, sdate, el.mean_anomaly, el.arg_periapsis, el.node_longitude, el.inclination, el.eccentricity, el.semi_major );
        //printf( "%s R %g date{%s} a %g e %g so %g lo %g inc %g ma %g T %g   \n", name.c_str(), radius, sdate,  &el.eccentricity, &el.semi_major,   el.arg_periapsis, el.node_longitude, el.inclination, el.mean_anomaly, el.epoch, el.period );


        //fgetc(inp);
        //double q=a*(1-e);
        //if(q<min){ // objects which are too far out won't be included
        //double y,m,d,epoch,period,m_now
        el.arg_periapsis *=DEG_2_RAD;
        el.node_longitude*=DEG_2_RAD;
        el.inclination   *=DEG_2_RAD;
        el.mean_anomaly  *=DEG_2_RAD;
        el.epoch  = string2julian_date( sdate );
        //period = pow(semi_major,1.5) * 365.25; // Kepler law
        el.period =  semi_major2period( el.semi_major ); // Kepler law
        orbit = new Orbit();
        orbit->fromElements( el );

        printf( "%s R,m %g,%g date{%s} e %g a %g so,lo,inc (%g,%g,%g) ma %g epoch %g T %g \n", name.c_str(), radius, mass, sdate,  el.eccentricity, el.semi_major,   el.arg_periapsis, el.node_longitude, el.inclination,   el.mean_anomaly, el.epoch, el.period );

    }

};


class BodyInteraction{ public:
    int i,j,kind;
};


class SpaceBodyIntegrator : public ODEderivObject, public ODEintegrator_RKF45 { public:
    SpaceBody * o;
    int ncenters;
    SpaceBody ** centers;

    //void bind( SpaceBody * o, ncenters ){}

    virtual void getDerivODE( double t, int n, double * Ys, double * dYs ){
        Vec3d p = *((Vec3d*)(Ys  ));
        Vec3d v = *((Vec3d*)(Ys+3));
        Vec3d f; f.set(0.0);
        for(int i=0; i<ncenters; i++){
            Vec3d d;  // d = ( interpolate spline )
            d.sub( p );
            double r2 = d.norm2();
            f.add( centralGravityForce( d, o->mass * centers[i]->mass ) );
        }
        if(o->trjThrust){
            //f.add(  );   interpolate trjThrust
        }
        (*(Vec3d*)(dYs  )) = v;
        (*(Vec3d*)(dYs+3)).set_mul(f,1.0/o->mass);
    }

    void evalTrj( int nt, double dt, int ncent ){
        derivObj = this;
        (*(Vec3d*)(Y  )) = o->pos;
        (*(Vec3d*)(Y+3)) = o->vel;
        for( int i=0; i< nt; i++ ){
            step_RKF45( dt );
            save_step();
            o->trjPos[i] = *((Vec3d*)Y);
        }
    }

};

#endif
