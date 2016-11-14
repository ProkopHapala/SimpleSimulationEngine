
#include "DynamicOpt.h" // THE HEADER

// ===============  MoveSteps

void DynamicOpt::move_LeapFrog(){
    double dtv = dt*fscale_safe;
	for ( int i=0; i<n; i++ ){
        vel[i] += invMasses[i]*dtv*force[i];
		pos[i] += dt*vel[i];
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

void DynamicOpt::move_FIRE(){
	double ff=0,vv=0,vf=0;
	for(int i=0; i<n; i++){
		double fi = force[i];
		double vi = vel[i];
		ff += fi*fi;
		vv += vi*vi;
		vf += vi*fi;
	}
	if( vf < 0.0 ){
		dt       = dt * fdec;
	  	damping  = damp_max;
		lastNeg  = 0;
		//cleanVel  ( );
		for(int i=0; i<n; i++){ vel[i] = kickStart*dt*force[i]; }
		//for(int i=0; i<n; i++){ vel[i] = dmax*force[i]*sqrt(1/ff)/dt_var; }
	}else{
		//double cf  =      damping * sqrt(vv/(ff+1e-8));
		double cf     =     damping * sqrt(vv/ff);
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

    if( ff > fmax*fmax ){ fscale_safe=fmax/sqrt(ff); }else{ fscale_safe=1; }
	move_LeapFrog();

	//printf( " %i f v vf  %f %f %f   dt damp  %f %f \n",  stepsDone,   sqrt(ff), sqrt(vv), vf/sqrt(vv*ff),   dt_var, damp_var  );
	//stepsDone++;
}

double DynamicOpt::optStep(){
	//cleanForce( );
	getForce( n, pos, force );
	switch( method ){
		case 0: move_LeapFrog();
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


