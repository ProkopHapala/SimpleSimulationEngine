
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

#include "appliedPhysics.h"

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

class SpaceBody : public PointBody  { public:

    std::string name;
    double radius;

    Vec3d * trjPos    = 0; // positions in some times
    //Vec3d * trjVel    = 0; // velocities in some times
    Vec3d * trjThrust = 0; // vector of thrust in time

    SpaceBody* orbCenter=0;

    inline Vec3d getThrust(int itrj, double du ){
        //if( trjThrust ){
            //printf( "%i %f   (%f,%f,%f)   (%f,%f,%f) \n", itrj, du, trjThrust[itrj].x, trjThrust[itrj].y, trjThrust[itrj].z,  trjThrust[itrj+1].x, trjThrust[itrj+1].y, trjThrust[itrj+1].z );
            return trjThrust[itrj]*(1-du) + trjThrust[itrj+1]*du;
        //}else{
        //    return (Vec3d){0.0,0.0,0.0};
        //}
    }
    //Vec3d getThrust(double t){};

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
