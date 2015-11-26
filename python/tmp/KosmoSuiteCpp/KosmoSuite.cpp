#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

#include "cpp/fastmath.h"
#include "cpp/Vec3.h"
#include "cpp/spline_hermite.h"
#include "cpp/optimizerDerivs.h"
#include "cpp/ODEintegrator.h"

const double GRAV_CONTS = 6.67384e-11;

void gravity_force( const Vec3d& r12, double GM12, Vec3d& f )   {
	double ir2 = 1.0d/r12.norm2( );
	double fr =  ir2*sqrt(ir2)*GM12;
	f.add_mul( r12, fr ); 
	//printf( "fr %e ir2 %e G %e ma %e mb %e \n",  fr,  ir2, GRAV_CONTS, ma, mb );
	//printf( " ir2  %e fr %e  f %e %e %e r12  %e %e %e \n", ir2,  fr,    f.x, f.y, f.z,    r12.x, r12.y, r12.z );
}

void lorenz_force( double q, const Vec3d& v, const Vec3d& B, Vec3d& f ){
	Vec3d df;
	df.set_cross( v,  B );
	f.add_mul   ( df, q );
}

void biot_savart_element( const Vec3d& R, const Vec3d& dI, Vec3d& B ){
	Vec3d dB;
	dB.set_cross( dI, R );
	double r2 = R.norm2();
	B.add_mul( dB, 1e-7 / ( r2 * sqrt(r2) ) ); 
}

#include "cpp/OrbitalUtils.h"

#include "cpp/Nbody.h"
#include "cpp/ShipAccel.h"
#include "cpp/elmag.h"


#include "cpp/fissionPulse.h"

extern "C"{

// ========= Nbody

void nbody_setup( int n, double * mass, double * poss, double * vs, double * errs ){
	NBody::setup( n, mass, poss, vs, errs );
}

void nbody_run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
	NBody::run( dt_start, dt_min, dt_max, nstep, tsIn, tsOut, poss, vs );
}

// ========= ShipAccel

void shipAccel_setup( double GM, double accel, double ivErr, double irErr, double * poss, double * vs ){
	ShipAccel::setup( GM, accel, irErr, ivErr, poss, vs );
}

void shipAccel_run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
	ShipAccel::run( dt_start, dt_min, dt_max, nstep, tsIn, tsOut, poss, vs );
}

// ========= ElMag

void elmag_sample ( int m, double * wheres, double * Bs, int n, double * ps, double * dIs ){
	ElMag::sample_field (  m, (Vec3d*)wheres, (Vec3d*)Bs, n, (Vec3d*)ps, (Vec3d*)dIs );
}

void elmag_setup( double charge, double mass,  double ivErr, double irErr, double * poss, double * vs, int ncoil, double * coil_ps, double * coil_dIs ){
	ElMag::setup( charge, mass, irErr, ivErr, poss, vs, ncoil, (Vec3d*)coil_ps, (Vec3d*)coil_dIs );
}

void elmag_run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
	ElMag::run( dt_start, dt_min, dt_max, nstep, tsIn, tsOut, poss, vs );
}

// ========= FissionPulse

void FissionPulse_set_params( double mass,double alphaDOF,double total_crossection,double fission_share, double back_scatter_share, double generation_rate,  double spontaneous_rate, double E_fission, double nmult, double nmult_spontal){
	FissionPulse::set_params( mass, alphaDOF, total_crossection, fission_share, back_scatter_share, generation_rate, spontaneous_rate, E_fission, nmult, nmult_spontal );
}

void FissionPulse_set_initial( double R0, double v0, double Nf0, double Nn0, double Q0 ){
	FissionPulse::set_initial( R0, v0, Nf0, Nn0, Q0 );
}

void FissionPulse_run_fixStep( double dt, int nstep, double * buff ){
	FissionPulse::run_fixStep( dt, nstep, buff );
}

void FissionPulse_set_trigger( double time_trigger, double R_trigger, double Nn_spark ){ 
	FissionPulse::set_trigger( time_trigger, R_trigger, Nn_spark );
}

}
