
#ifndef  Opt3d_h
#define  Opt3d_h

//#include "Vec2.h"
#include "Vec3.h"

class Opt3d{ public:
// "Fast Inertial Realxation Engine" according to
// Bitzek, E., Koskinen, P., Gähler, F., Moseler, M. & Gumbsch, P. Structural relaxation made simple. Phys. Rev. Lett. 97, 170201 (2006).
// Eidel, B., Stukowski, A. & Schröder, J. Energy-Minimization in Atomic-to-Continuum Scale-Bridging Methods. Pamm 11, 509–510 (2011).
// http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf

    // parameters
    double fTinc = 1.1;   // factor by which time step is increased if going downhill
    double fTdec = 0.5;   // factor by which timestep is decreased if going uphill
    double fDamp = 0.95;  // rate of decrease of damping when going downhill
    //
    double dtmax   = 0.2;   // maximal timestep
    double dtmin   = 0.02;
    double dampMax = 0.5;   // default damping

    // variables
    double dt     = dtmax;    // time-step ( variable
    double damp   = dampMax;  // damping  ( variable

    inline void setup( double dtmax_, double dtmin_, double dampMax_=0.5){
        dtmax   = dtmax_;
        dtmin   = dtmin_;
        dampMax = dampMax_;
        dt      = dtmax;
        damp    = dampMax;
    }

    inline void  moveMD( const Vec3d& f, Vec3d& p, Vec3d& v ){
        v.mul( 1 - damp  );
        v.add_mul( f, dt );
        p.add_mul( v, dt );
    }

    // relaxation step using FIRE algorithm
    inline void move( const Vec3d& f, Vec3d& p, Vec3d& v ){
        double ff = f.norm2();
        double vv = v.norm2();
        double vf = f.dot(v);
        if( vf < 0 ){ // if velocity along direction of force
            v.set( 0.0d );
            dt   = fmax( dt * fTdec, dtmin );
            damp = dampMax;
        }else{       // if velocity against direction of force
            double cf  =     damp * sqrt(vv/ff);
            double cv  = 1 - damp;
            v.mul    ( cv    );
            v.add_mul( f, cf );	// v = cV * v  + cF * F
            dt    = fmin( dt * fTinc, dtmax );
            damp  = damp     * fDamp;
        }
        // normal leap-frog times step
        v.add_mul( f , dt );
        p.add_mul( v , dt );
    }

};

#endif


