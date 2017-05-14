#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cstdint>

/*
#include "../../common/math/fastmath.h"
#include "../../common/math/Vec3.h"
#include "../../common/math/spline_hermite.h"
#include "../../common/optimization/optimizerDerivs.h"
#include "../../common/dynamics/ODEintegrator.h"
#include "../../common/dynamics/Shock1D.h"
*/

#include "fastmath.h"
#include "Vec3.h"
#include "spline_hermite.h"
//#include "optimizerDerivs.h"
#include "ODEintegrator.h"
//#include "Shock1D.h"

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
#include "cpp/SpaceLaunchODE.h"
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

// ========= SpaceLaunchODE

ODEintegrator_RKF45  odeint;

void SpaceLaunchODE_init( SpaceLaunchODE::Planet *planet_, SpaceLaunchODE::Rocket *rocket_, SpaceLaunchODE::Aerodynamics *aero_, SpaceLaunchODE::Atmosphere *atmosphere_){
    using namespace SpaceLaunchODE;
    planet = *planet_; rocket = *rocket_; aero = *aero_; atmosphere = *atmosphere_;
    /*
    printf( "%i \n", planet_ );
    printf( "%i \n", rocket_ );
    printf( "%i \n", aero_ );
    printf( "%i \n", atmosphere_ );


    printf( "%i \n", atmosphere_->n );
    printf( "%i %g %g %g %g\n", atmosphere_->n, atmosphere_->dh, atmosphere_->inv_dh, atmosphere_->hmax, atmosphere_->rho_CPs[2] );
    printf( "%i %g %g %g %g\n", atmosphere.n, atmosphere.dh, atmosphere.inv_dh, atmosphere.hmax, atmosphere.rho_CPs[1] );

    printf( "HEY  atmosphere.n %i \n", atmosphere.n );

    for(int i=0; i<atmosphere.n; i++){
        printf( "%i %f\n", i,atmosphere.rho_CPs[i] );
    }
    */

    odeint.reallocate( 7 );
    odeint.dt_max    = 0.005;
    odeint.dt_min    = 0.0001;
    odeint.dt_adapt  = 0.001;
    odeint.getDerivs = getODEDerivs;

}

int SpaceLaunchODE_run( int nMax, int nMaxIters, SpaceLaunchODE::Launch *launch_, SpaceLaunchODE::LogTrig *logTrig_, double * outbuff ){
    using namespace SpaceLaunchODE;
    launch = *launch_; logTrig = *logTrig_;
    logbuff = (LogVars*)outbuff;
    logFunc = logfunc_default;
    //logTrig.tmax   = 0;
    //logTrig.t_trig = 1e+300;
    //logTrig.on     = false;

    ((Vec3d*)(odeint.invMaxYerr  ))->set(1/1e-1);
    ((Vec3d*)(odeint.invMaxYerr+3))->set(1/1e-2);
    odeint.invMaxYerr[6] = 1/1e+1;

    Vec3d  * pos  = (Vec3d*)(odeint.Y   );
    Vec3d  * vel  = (Vec3d*)(odeint.Y+3 );
    Vec3d  * acc  = (Vec3d*)(odeint.dY+3);
    double * mass = odeint.Y+6;

    pos->set(0.0,0.0,0.0); vel->set(0.0,0.0,0.0);
    //pos->set(0.0,100.0e+3,0.0); vel->set(7.9e+3,0.0,0.0); satelite
    *mass = rocket.mass_initial;

    /*
    printf( "rocket %g %g %g %g %g %g \n",         rocket.AeroArea, rocket.dm_F, rocket.mass_empty, rocket.mass_initial, rocket.thrust_full, rocket.vexhaust );
    printf( "planet %g %g (%g,%g,%g) \n",          planet.GM, planet.R, planet.pos.x, planet.pos.y, planet.pos.z );
    printf( "launch %g %g %g %g \n",               launch.uT, launch.inv_uT, launch.hTarget, launch.vTarget  );
    printf( "logTrig %i %i %g %g %g \n",           logTrig.i, logTrig.on, logTrig.dt_trig, logTrig.t_trig, logTrig.tmax  );

    printf( ">>pos<< (%g,%g,%g) vel (%g,%g,%g) mass %g \n", pos->x,pos->y,pos->z, vel->x,vel->y,vel->z, *mass );
    for(int i=0; i<launch.n; i++){
        printf( "%i : %g %g (%g,%g,%g) \n", i, launch.thetaCPs[i], launch.thrustCPs[i], launch.dirCPs[i].x, launch.dirCPs[i].y, launch.dirCPs[i].z );
    }
    printf( "xMax %i nMaxIters %i \n", nMax, nMaxIters );
    //exit(0);
    */

    int nstep = 0;
    for(int i=0; i<nMax; i++){
        //odeint.step( 0.1d );
        //odeint.adaptive_step_RKF45( );
        //nstep = odeint.integrate_adaptive( odeint.dt_adapt, odeint.t+0.2d );
        nstep += odeint.integrate_adaptive( odeint.dt_adapt, odeint.t+5.0d );
        if( ( logVars.vx >= launch.vTarget ) && ( logVars.h  >= launch.hTarget ) ) break;
        if( ( logVars.t > launch.tMax ) || ( nstep > nMaxIters ) ){ logTrig.i*=-1; break; }
        if( logTrig.i >= nMax ) break;
        //logbuff[i].t = i;
    }
    return logTrig.i;
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
