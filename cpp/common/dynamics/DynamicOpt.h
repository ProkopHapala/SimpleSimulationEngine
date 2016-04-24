
#ifndef DynamicOpt_h
#define DynamicOpt_h

#include <math.h>

typedef void (*ForceFunction)( int n, double * xs, double * dfs );

/*
const double RELAX_damping   = 0.1;
const double RELAX_dt        = 0.1;
const double RELAX_convF2    = 0.000001;
const int    RELAX_maxIters  = 1000;
*/

class DynamicOpt{
	public:
	// variables
	int n;
	double * pos;
	double * vel;
	double * force;

	// parameters
	double dt           = 0.05d;
	double damping      = 0.1d;

	// FIRE
	int    minLastNeg   = 5;
	double finc         = 1.1d;
	double fdec         = 0.5d;
	double falpha       = 0.98d;
	double kickStart    = 1.0d;

	double dt_max       = dt;
	double damp_max     = damping;

	int    lastNeg      = 0;

	// other
	int method    = 2;
	int stepsDone = 0;
	double t      = 0.0d;

	ForceFunction getForce;

	// ==== function declarations

	void   move_LeapFrog();
	void   move_MDquench();
	void   move_FIRE();
	double optStep();
	bool   optimize( double convF, int nMaxSteps );

	double getFmaxAbs( );
	double getFsqSum( );

	// ==== inline functions

	inline void bindArrays( int n_, double * pos_, double * vel_, double * force_ ){
		n = n_;
		pos   = pos_;
		vel   = vel_;
		force = force_;
	}

	inline void allocate( int n_ ){
		n = n_;
		pos   = new double[n];
		vel   = new double[n];
		force = new double[n];
	}

	inline void deallocate( ){
		delete pos;
		delete vel;
		delete force;
	}

	inline void cleanForce( ){  for(int i=0; i<n; i++){ force[i]=0; } }
	inline void cleanVel  ( ){  for(int i=0; i<n; i++){ vel  [i]=0; } }

	inline void initOpt( double dt_, double damp_ ){
		dt      = dt_max   = dt_;
		damping = damp_max = damp_;
		cleanForce( );
		cleanVel  ( );
	}

};

#endif
