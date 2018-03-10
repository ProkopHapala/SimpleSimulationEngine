
#ifndef  SpaceWorld_h
#define  SpaceWorld_h

#include <vector>

#include "appliedPhysics.h"
#include "SpaceBodies.h"

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

    void addPlanet( std::string name, double mass, double radius, Vec3d pos, Vec3d vel ){
        planets.emplace_back();
        SpaceBody& p = planets.back();
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

        double u  =(t-trj_t0)/trj_dt;
        int itrj  = (int)u;
        double du = u-itrj;

        for( int i=0; i<nobj; i++ ){
            Vec3d f = (Vec3d){0.0,0.0,0.0};
            int i6 = i*6;
            Vec3d  p   = *(Vec3d*)(Ys+i6);
            Vec3d* pvs =  (Vec3d*)(Ys   );
            //double m = (i<nplanet)? planets[i].mass : ships[i-nplanet].mass;
            double m = masses[i];
            for( int j=0; j<nplanet; j++ ){
                if(i!=j){;
                //gravity( pvs-p, m*planets[j].mass );
                //printf( " %i %i (%f,%f,%f)\n", i, j,  p.x, p.y, p.z );
                    f.add( centralGravityForce( (*pvs)-p, m*masses[j] ) );
                //printf( " %i %i (%f,%f,%f)\n", i, j,  f.x, f.y, f.z );
                }
                pvs+=2;
            }
            if(i>nplanet){
                f.add( ships[i-nplanet].getThrust(itrj,du) );
            }
            (*(Vec3d*)(dYs+i6  )) = (*(Vec3d*)(Ys+i6+3));
            (*(Vec3d*)(dYs+i6+3)) = f * (1/masses[i]);

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
            if(b.trjPos) delete [] b.trjPos; b.trjPos = new Vec3d[n];
        }
        for( SpaceBody& b : ships ){
            if(b.trjPos   ) delete [] b.trjPos;    b.trjPos    = new Vec3d[n];
            if(b.trjThrust) delete [] b.trjThrust; b.trjThrust = new Vec3d[n];
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
        if(n>trj_n ) n=trj_n;
        int nplanet = planets.size();
        int nship   = ships.size();
        // --- initialize
        ode.t = trj_t0;
        objects2ode();
        // --- run simulation
        for(int i=0; i<n; i++){
            //printf( "ODE step %i \n ", i );
            ode.step( dt );
            ode2trjs(i);
        }
    }

};

#endif
