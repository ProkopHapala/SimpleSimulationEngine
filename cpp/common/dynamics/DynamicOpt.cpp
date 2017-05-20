
#include "DynamicOpt.h" // THE HEADER

#include <cstdio>// DEBUG

// ===============  MoveSteps

/*
void DynamicOpt::move_LeapFrog(){
    double dtv = dt*fscale_safe;
	for ( int i=0; i<n; i++ ){
        vel[i] += invMasses[i]*dtv*force[i];
		pos[i] += dt*vel[i];
	}
	stepsDone++;
	t += dt;
}
*/

void DynamicOpt::move_LeapFrog(double dt_loc){
    //double dt_ = dt*fscale_safe;
	for ( int i=0; i<n; i++ ){
        vel[i] += invMasses[i]*dt_loc*force[i];
		pos[i] += dt_loc*vel[i];
	}
	stepsDone++;
	t += dt_loc;
}

void DynamicOpt::move_LeapFrog_vlimit(){
    double dtv = dt*fscale_safe;
    double vmax = 0.0d;
	for ( int i=0; i<n; i++ ){
        double v = vel[i] + invMasses[i]*dtv*force[i];
        vmax = fmax( fabs(v), vmax );
        vel[i]=v;
    }
    double dtp = dt;
    if( vmax>v_limit ) dtp=v_limit/vmax;
    //printf("vmax %g dtp  %g dtv %g\n", vmax);
    for ( int i=0; i<n; i++ ){
		pos[i] += dtp*vel[i];
	}
	stepsDone++;
	t += dt;
}

void DynamicOpt::move_MDquench(){
	double cdamp = 1 - damping;
	double dtv = dt*fscale_safe;
	for ( int i=0; i<n; i++ ){
		vel[i]  = cdamp*vel[i] + dtv*force[i]*invMasses[i];
		pos[i] += dt*vel[i];
	}
	stepsDone++;
	t += dt;
}

double DynamicOpt::move_FIRE(){
	double ff=0,vv=0,vf=0;
	for(int i=0; i<n; i++){
		double fi = force[i];
		double vi = vel[i];
		ff += fi*fi;
		vv += vi*vi;
		vf += vi*fi;
	}
	if( vf < 0.0 ){
		//dt       = dt * fdec;
		dt       = fmax( dt * fdec, dt_min );
	  	damping  = damp_max;
		lastNeg  = 0;
		//cleanVel  ( );
		for(int i=0; i<n; i++){ vel[i] = kickStart*dt*force[i]; }
		//for(int i=0; i<n; i++){ vel[i] = dmax*force[i]*sqrt(1/ff)/dt_var; }
	}else{
		double cf  =      damping * sqrt(vv/(ff+ff_safety));
		//double cf     =     damping * sqrt(vv/ff);
		double cv     = 1 - damping;
		for(int i=0; i<n; i++){
			vel[i]    = cv * vel[i]  + cf * force[i];
		}
		if( lastNeg > minLastNeg ){
			dt        = fmin( dt * finc, dt_max );
			damping   = damping  * falpha;
		}
		lastNeg++;
	}

	double dt_=dt;
    if( ff > f_limit*f_limit ){
        double f = sqrt(ff);
        dt_*=sqrt(f_limit/f);
        printf( "force too large: %g => limit dt: %g \n", f, dt_ );
    };
    move_LeapFrog( dt_ );
	//move_LeapFrog();
	//move_LeapFrog_vlimit();  // this does not seem to help

	//printf( " %i f v vf  %f %f %f   dt damp  %f %f \n",  stepsDone,   sqrt(ff), sqrt(vv), vf/sqrt(vv*ff),   dt_var, damp_var  );
	//stepsDone++;
	return ff;
}

double DynamicOpt::optStep(){
	if(getForce){
        cleanForce( );
        getForce( n, pos, force );
	}
	switch( method ){
		case 0: move_LeapFrog(dt);
		case 1: move_MDquench();
		case 2: move_FIRE();
	}
	return getFmaxAbs( );
}

bool DynamicOpt::optimize( double convF, int nMaxSteps ){
	for( int i=0; i<nMaxSteps; i++ ){
		double f = optStep();
		if( f < convF ) return true;
	}
	return false;
}

// =============== common rutines

double DynamicOpt::getFmaxAbs( ){
	double fmax = 0;
	for(int i=0; i<n; i++){
		double fi = fabs( force[i] );
		fmax=(fi>fmax)?fi:fmax;
	}
	return fmax;
}

double DynamicOpt::getFsqSum( ){
	double ff = 0;
	for(int i=0; i<n; i++){
		double fi = force[i];
		ff += fi*fi;
	}
	return ff;
}

/*
void DynamicOpt::autoTimeStep( double dt0, double v0, double f0 ){
	double fmax = 0;
	double vmax = 0;
	for(int i=0; i<n; i++){
		double fi = fabs( force[i] );  fmax=(fi>fmax)?fi:fmax;
		double vi = fabs( vel[i]   );  vmax=(vi>vmax)?vi:vmax;
	}
	//  f*dt < dvmax
	//  v*dt < dxmax
	//  dt = 1/sqrt( k )
	//  k  = f/x
	dt = dt0/( 1 + vmax/v0 + fmax/f0 );
	printf( " %f %f %f \n", dt, vmax, fmax );
}
*/


