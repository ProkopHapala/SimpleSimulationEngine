
#ifndef  SpaceWorld_h
#define  SpaceWorld_h

#include <vector>

#include "appliedPhysics.h"
#include "SpaceBodies.h"


namespace Time{
    static const double
        minute=60,
        hour  =3600,
        day   =86400,
        week  =day*7,
        month =2629746,
        year  =31556952;

    static const int nNice = 9;
    static const double niceUnits[nNice] = {1,10,minute,minute*10,hour,day,week,month,year};

    /*
    int toStr(double nsec,char* out){
        double absval = fabs(nsec);
        if     ( absval>year   ){ return sprintf(out, "%g y",  nsec/year   ); }
        else if( absval>month  ){ return sprintf(out, "%2.2f months", nsec/month  ); }
        else if( absval>week   ){ return sprintf(out, "%2.2f weeks",  nsec/week   ); }
        else if( absval>hour   ){ return sprintf(out, "%2.2f hours",  nsec/hour   ); }
        else if( absval>minute ){ return sprintf(out, "%2.2f minutes",nsec/minute ); }
        else                    { return sprintf(out, "%g sec",    nsec        ); }
    }
    */

    int toStr(double nsec,char* out){
        double absval = fabs(nsec);
        if     ( absval>(100*year)){ return sprintf(out, "%gy",    nsec/year   ); }
        if     ( absval>year      ){ return sprintf(out, "%2.2fy", nsec/year   ); }
        else if( absval>day       ){ return sprintf(out, "%2.2fd", nsec/day    ); }
        else if( absval>hour      ){ return sprintf(out, "%2.2fh", nsec/hour   ); }
        else if( absval>minute    ){ return sprintf(out, "%2.2fm", nsec/minute ); }
        else if( absval>0.1       ){ return sprintf(out, "%2.2fs", nsec        ); }
        else                       { return sprintf(out, "%gs",    nsec        ); }
    }

    double niceUnitAbove(double val){
        int inice = 0;
        if     (val<niceUnits[0      ] ){ return  pow(10.0,(int)log10(val));              }
        else if(val>niceUnits[nNice-1] ){ return (pow(10.0,(int)log10(val/year)+1)*year); }
        else{
            for(int i=0;i<nNice;i++){ if(val<niceUnits[i]){ inice=i; break; } }
            return niceUnits[inice];
        }
    };

};


class SpaceWorld : public ODEderivObject { public:
    std::vector<SpaceBody> planets;
    std::vector<SpaceBody> ships;
    //std::vector<SpaceBody> projectiles;

    //ODEintegrator_RKF45 odePlanets;
    //ODEintegrator_RKF45 odeShip;
    ODEintegrator_RKF45 ode;

    int    trj_n    = 0;
    int    trj_i    = 0;
    double trj_t0   = 0;
    double trj_dt   = 0.1;
    //Vec3d ** trj_Pos = 0;
    double * masses = 0;

    inline double trjTime2index(double t, int& itrj){ double u=(t-trj_t0)/trj_dt; itrj=(int)u; return u-itrj; };
    inline double ind2time     ( double ind        ){ return trj_dt*ind+trj_t0; };

    void intertialTansform( Vec3d p0, Vec3d v0 ){
        for(SpaceBody& b : planets ){ b.pos.add(p0); b.vel.add(v0); };
        for(SpaceBody& b : ships   ){ b.pos.add(p0); b.vel.add(v0); };
    };

    void addPlanet( std::string name, double mass, double radius, Vec3d pos, Vec3d vel ){
        planets.emplace_back();
        SpaceBody& p = planets.back();
        p.name = name;
        p.mass = mass;
        p.radius = radius;
        p.pos = pos;
        p.vel = vel;
    }

    void addShip( std::string name, double mass, double radius, Vec3d pos, Vec3d vel ){
        ships.emplace_back();
        SpaceBody& p = ships.back();
        p.name = name;
        p.mass = mass;
        p.radius = radius;
        p.pos = pos;
        p.vel = vel;
    }

    virtual void getDerivODE( double t, int n, double * Ys, double * dYs ){
        int nplanet = planets.size();
        int nship   = ships.size();
        int nobj    = nplanet + nship;
        //Vec3d* Yplanet  = (Vec3d*)(Ys           );
        //Vec3d* Yship    = (Vec3d*)(Ys+6*nplanet );
        //Vec3d* dYplanet = (Vec3d*)(dYs          );
        //Vec3d* dYship   = (Vec3d*)(dYs+6*nplanet);

        int itrj; double du = trjTime2index(t,itrj);
        //printf("%i %f \n", itrj, du );

        for( int i=0; i<nobj; i++ ){
            Vec3d f = (Vec3d){0.0,0.0,0.0};
            int i6 = i*6;
            Vec3d  p   = *(Vec3d*)(Ys+i6);
            Vec3d* pvs =  (Vec3d*)(Ys   );
            //double m = (i<nplanet)? planets[i].mass : ships[i-nplanet].mass;
            //double m = masses[i];
            for( int j=0; j<nplanet; j++ ){
                if(i!=j){;
                //gravity( pvs-p, m*planets[j].mass );
                //printf( " %i %i (%f,%f,%f)\n", i, j,  p.x, p.y, p.z );
                //f.add( centralGravityForce( (*pvs)-p, m*masses[j] ) );
                f.add( centralGravityForce( (*pvs)-p, masses[j] ) );
                //printf( " %i %i (%f,%f,%f)\n", i, j,  f.x, f.y, f.z );
                }
                pvs+=2;
            }
            if(i>=nplanet){
                f.add( ships[i-nplanet].getThrust(itrj,du) );
                //Vec3d T = ships[i-nplanet].getThrust(itrj,du); f.add( T );
                //printf( "%s T=( %f %f %f ) \n", ships[i-nplanet].name.c_str(), T.x,T.y,T.z );
            }
            (*(Vec3d*)(dYs+i6  )) = (*(Vec3d*)(Ys+i6+3));
            //(*(Vec3d*)(dYs+i6+3)) = f * (1/masses[i]);
            (*(Vec3d*)(dYs+i6+3)) = f;

        }

        /*
        printf( "getDerivODE  Y : " );
        for( int i=0; i<n; i++ ) printf(" %f ",  Ys[i]);
        printf( "\n" );

        printf( "getDerivODE dY : " );
        for( int i=0; i<n; i++ ) printf(" %f ",  dYs[i]);
        printf( "\n" );
        */

    }

    void allocateODE(){
        ode.derivObj = this;
        int n = planets.size() + ships.size();
        ode.reallocate( n*6 );
        //if(masses) delete [] masses; masses[]=new ;
        _realloc(masses,n);
    }

    void allocateTrjs(int n){
        trj_n = n;
        for( SpaceBody& b : planets ){
            //if(b.trjPos) delete [] b.trjPos; b.trjPos = new Vec3d[n];
            _realloc(b.trjPos,n+3);
        }
        for( SpaceBody& b : ships ){
            //if(b.trjPos   ) delete [] b.trjPos;    b.trjPos    = new Vec3d[n];
            //if(b.trjThrust) delete [] b.trjThrust; b.trjThrust = new Vec3d[n];
            _realloc(b.trjPos,n+3);
            _realloc(b.trjThrust,n+3);
            for(int i=0; i<n*3; i++) ((double*)(b.trjThrust))[i] = 0;
        }
    }

    void objects2ode(){
        Vec3d  * pvs =  (Vec3d*)(ode.Y);
        double * ms  = masses;
        for( SpaceBody& b : planets ){
            pvs[0] = b.pos;
            pvs[1] = b.vel;
            *ms     = b.mass;
            pvs+=2; ms++;
        }
        for( SpaceBody& b : ships ){
            pvs[0] = b.pos;
            pvs[1] = b.vel;
            *ms     = b.mass;
            pvs+=2; ms++;
        }
    }

    void ode2objects(){
        Vec3d* pvs =  (Vec3d*)(ode.Y);
        for(int i=0; i<planets.size(); i++){
            planets[i].pos = pvs[0];
            planets[i].vel = pvs[1];
            pvs+=2;
        }
        for(int i=0; i<ships.size(); i++){
            ships[i].pos = pvs[0];
            ships[i].vel = pvs[1];
            pvs+=2;
        }
    }

    void ode2trjs(int j){
        Vec3d* pvs =  (Vec3d*)(ode.Y);
        for(int i=0; i<planets.size(); i++){
            planets[i].trjPos[j] = pvs[0];
            //printf( "step %i planets[%i] (%f,%f,%f) \n", j, i, planets[i].trjPos[j].x, planets[i].trjPos[j].y, planets[i].trjPos[j].z );
            pvs+=2;
        }
        for(int i=0; i<ships.size(); i++){
            ships[i].trjPos[j]   = pvs[0];
            pvs+=2;
        }
    }

    void predictTrjs( int n, double dt ){
        trj_dt = dt;
        ode.t = trj_t0;
        if(n>trj_n ) n=trj_n;
        int nplanet = planets.size();
        int nship   = ships.size();
        // --- initialize
        objects2ode();
        // --- run simulation
        for(int i=0; i<n; i++){
            //printf( "ODE step %i \n ", i );
            ode.step( dt );
            ode2trjs(i);
        }
    }

    int load_astorb(char* fname, int n_reserve){
        //printf( "load_astorb \n" );
        FILE * pFile;
        const int nbuff = 1024;
        char str[nbuff];
        pFile = fopen ( fname, "r");
        if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
        planets.reserve( n_reserve );
        int n=0;
        //printf( "TO FGETS \n" );
        while ( fgets( str , nbuff, pFile) != NULL ){
            //printf( "%s", str );
            SpaceBody b;
            b.fromString_astorb(str);
            planets.push_back(b);
            n++;
        }
        fclose(pFile);
        return n;
    }

    int pickPlanet( const Vec3d& ro, const Vec3d& rd, double epoch ){
        double t_ray;
        double r2min = 1e+300;
        int imin     = -1;
        for(int i=0; i<planets.size(); i++){
            if( !planets[i].orbit ) continue;
            Vec3d p = planets[i].orbit->pointAtEpoch( epoch );
            double r2 = rayPointDistance2( ro, rd, p, t_ray );
            if(r2<r2min){ imin=i; r2min=r2; }
        }
        printf( "pickPlanet %i %g  ro(%g,%g,%g) rd(%g,%g,%g) \n", imin, r2min,   ro.x,ro.y,ro.z,    rd.x,rd.y,rd.z  );
        return imin;
    }

};

#endif
