
#ifndef DynamicOpt_h
#define DynamicOpt_h

//#include <cstddef>
#include <math.h>
#include "macroUtils.h"

typedef void (*ForceFunction)( int n, double * xs, double * dfs );

class DynamicOpt{ public:
    // variables
    int n=0;
    double * pos       = 0;
    double * vel       = 0;
    double * force     = 0;
    double * invMasses = 0;

    // parameters
    double dt           = 0.1d;
    double damping      = 0.2d;

    double f_limit      = 10.0;
    double v_limit      = 10.0;
    //double l_limit      = 0.1;
    double dr_limit     = 0.1;
    //double fscale_safe  = 1;
    double scale_dt  = 1;
    double ff_safety    = 1e-32;

    // FIRE
    int    minLastNeg   = 5;
    double finc         = 1.1d;
    double fdec         = 0.5d;
    double falpha       = 0.98d;
    double kickStart    = 1.0d;

    double dt_max       = dt;
    double dt_min       = 0.1 * dt;
    double damp_max     = damping;

    int    lastNeg      = 0;

    // other
    int method    = 2;
    int stepsDone = 0;
    double t      = 0.0d;

    ForceFunction getForce = 0;

    // ==== function declarations

    void   move_LeapFrog( double dt_loc );
    //void   move_LeapFrog_vlimit();
    void   move_GD      ( double dt_loc );
    //double move_GD_safe ( double dt_loc );
    void   move_MD( double dt_loc, double damp);
    //double move_MD_safe ( double dt_loc );
    void   move_MDquench(){move_MD(dt,damping);};
    double move_FIRE();
    double optStep();
    bool   optimize( double convF, int nMaxSteps );

    double getFmaxAbs( );
    double getFsqSum( );

    // ==== inline functions

    inline void setInvMass(double invM){  if(invMasses==0){ _realloc(invMasses,n);}  for(int i=0;i<n;i++){invMasses[i]=invM;} };

    inline void bindArrays( int n_, double * pos_, double * vel_, double * force_, double * invMasses_ ){
        n = n_; pos=pos_;  vel=vel_; force=force_; invMasses=invMasses_;
    }

    inline void bindOrAlloc( int n_, double * pos_, double * vel_, double * force_, double * invMasses_ ){
        n = n_;
        if(pos_  ==0){ _realloc(pos  ,n); }else{ pos   = pos_;   };
        if(vel_  ==0){ _realloc(vel  ,n); }else{ vel   = vel_;   };
        if(force_==0){ _realloc(force,n); }else{ force = force_; };
        //if(invMasses_==0) { _realloc(invMasses,n); setInvMass(1.0); }else{ invMasses=invMasses_; }
        if(invMasses_==0) setInvMass(1.0);
    }

    inline void realloc( int n_ ){
        n = n_;
        _realloc(pos    ,n);
        _realloc(vel    ,n);
        _realloc(force  ,n);
        _realloc(invMasses,n);
    }

    inline void dealloc( ){
        _dealloc(pos);
        _dealloc(vel);
        _dealloc(force);
        _dealloc(invMasses);
    }

    inline void cleanForce( ){  for(int i=0; i<n; i++){ force[i]=0; } }
    inline void cleanVel  ( ){  for(int i=0; i<n; i++){ vel  [i]=0; } }

    // f ~ sqrt(k/m)    dpos ~ f*dt*dt
    inline double limit_dt_x2 (double xx,double xmax){ double sc=1.0; if( xx > (xmax*xmax) ){ sc= fmin( sc, sqrt(xmax/sqrt(xx)) ); }; return sc;       }
    //inline double limit_dt_x2 (double xx,double xmax){ double sc=1.0; if( xx > (xmax*xmax) ){ sc= fmin( sc, xmax/sqrt(xx) ); }; return sc;       }
    inline double limit_dt_vf2(double ff, double vv ){ scale_dt=fmin(limit_dt_x2(ff,f_limit),limit_dt_x2(vv,v_limit));         return scale_dt; }

    inline void initOpt( double dt_, double damp_ ){
        dt      = dt_max   = dt_;  dt_min=0.1*dt_max;
        damping = damp_max = damp_;
        cleanForce( );
        cleanVel  ( );
    }

};

#endif
