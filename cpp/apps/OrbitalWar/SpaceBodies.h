
#ifndef  SpaceBodies_h
#define  SpaceBodies_h

#include <vector>
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Body.h"
#include "ODEintegrator.h"

#include "appliedPhysics.h"

/*
class SpaceCraftBody : public RigidBody{ public:
}
*/

class SpaceBody : public PointBody  { public:
    Vec3d * trjPos    = 0; // positions in some times
    //Vec3d * trjVel    = 0; // velocities in some times
    Vec3d * trjThrust = 0; // vector of thrust in time

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
            f.add( gravity( d, o->mass * centers[i]->mass ) );
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
