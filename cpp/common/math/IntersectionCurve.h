
#ifndef  IntersectionCurve_h
#define  IntersectionCurve_h

#include "Vec2.h"
#include "Vec3.h"


typedef  void (*FieldGrad)(const Vec3d& p, Vec3d& f);

namespace FIRE{
// "Fast Inertial Realxation Engine" according to
// Bitzek, E., Koskinen, P., Gähler, F., Moseler, M. & Gumbsch, P. Structural relaxation made simple. Phys. Rev. Lett. 97, 170201 (2006).
// Eidel, B., Stukowski, A. & Schröder, J. Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods. Pamm 11, 509–510 (2011).
// http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf

    // parameters
    double finc    = 1.1;             // factor by which time step is increased if going downhill
    double fdec    = 0.5;             // factor by which timestep is decreased if going uphill
    double falpha  = 0.99;            // rate of decrease of damping when going downhill
    double dtmax   = RELAX::dt;       // maximal timestep
    double acoef0  = RELAX::damping;  // default damping

    // variables
    double dt      = dtmax;           // time-step ( variable
    double acoef   = acoef0;          // damping  ( variable

    inline void setup( double dt_, double damping){
        dtmax   = RELAX::dt;
        acoef0  = RELAX::damping;
        dt      = dtmax;
        acoef   = acoef0;
    }

    inline void  moveMD( const Vec3d& f, Vec3d& r, Vec3d& v ){
        v.mul( 1 - damping );
        v.add_mul( f, dt );
        r.add_mul( v, dt );
    }

    // relaxation step using FIRE algorithm
    inline void move( const Vec3d& f, Vec3d& r, Vec3d& v ){
        double ff = f.norm2();
        double vv = v.norm2();
        double vf = f.dot(v);
        if( vf < 0 ){ // if velocity along direction of force
            v.set( 0.0d );
            dt    = dt * fdec;
              acoef = acoef0;
        }else{       // if velocity against direction of force
            double cf  =     acoef * sqrt(vv/ff);
            double cv  = 1 - acoef;
            v.mul    ( cv );
            v.add_mul( f, cf );	// v = cV * v  + cF * F
            dt     = fmin( dt * finc, dtmax );
            acoef  = acoef * falpha;
        }
        // normal leap-frog times step
        v.add_mul( f , dt );
        r.add_mul( v , dt );
    }

};


class Opt3d{
    double dt;
    Vec3d p;
    Vec3d v;


};


class IntersectionCurve{
    FieldGrad field1;
    FieldGrad field2;

    int maxRelaxStep=100;
    Vec3d p;
    Vec3d op;  // dimer direction

    inline moveMD( Vec3d& p, Vec3d& p, const Vec3d& f ){

    }

    void relax( ){
        Vec3d v=Vec3dZero;
        for(int i=0; i<maxRelaxStep; i++){
            Vec3d f1,f2;
            field1(p,f1);
            field2(p,f2);


        }
    }


    inline trace( int nmax, double Fconv=1e-6 ){

    }


}

#endif


