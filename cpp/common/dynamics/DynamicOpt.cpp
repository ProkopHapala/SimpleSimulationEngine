
#include "DynamicOpt.h" // THE HEADER

//#include <cstdio>// DEBUG
//#include "VecN.h"

#include "fastmath.h"

// ===============  MoveSteps

void DynamicOpt::move_LeapFrog(double dt_loc){
    //double dt_ = dt*fscale_safe;
    for ( int i=0; i<n; i++ ){
        //printf( "i %i v %g f %g p %g iM %g \n", i, vel[i],force[i],pos[i],invMasses[i]  );
        vel[i] += force[i]*invMasses[i]*dt_loc;
        pos[i] += vel[i]*dt_loc;
    }
    stepsDone++;
    t += dt_loc;
}

/*
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
*/

void DynamicOpt::move_GD(double dt_loc){
    //double dt_ = dt*fscale_safe;
    for ( int i=0; i<n; i++ ){
        pos[i] += force[i]*dt_loc;
    }
    stepsDone++;
    t += dt_loc;
}

/*
double DynamicOpt::move_GD_safe(double dt_loc){
    double fmax = VecN::absmax(n,force);
    scale_dt = fmin( 1, f_limit/fmax );
    dt_loc*=scale_dt;
    move_GD(dt_loc);
    stepsDone++;
    t += dt_loc;
    return fmax;
}
*/

void DynamicOpt::move_MD(double dt_loc,double damp){
    double cdamp = 1 - damp;
    for ( int i=0; i<n; i++ ){
        vel[i]  = cdamp*vel[i] + force[i]*invMasses[i]*dt_loc;
        pos[i] += vel[i]*dt_loc;
    }
    stepsDone++;
    t += dt;
}

/*
double DynamicOpt::move_MD_safe(double dt_loc){
    double fmax = VecN::absmax(n,force);
    scale_dt = fmin(1,fmin( v_limit/VecN::absmax(n,vel), f_limit/fmax ));
    dt_loc*=scale_dt;
    move_MD(dt_loc);
    return fmax;
}
*/


/*
double DynamicOpt::move_FIRE(){
    double ff=0,vv=0,vf=0;
    for(int i=0; i<n; i++){
        double fi = force[i];
        double vi = vel[i];
        ff += fi*fi;
        vv += vi*vi;
        vf += vi*fi;
    }

    bool bOverLimit = ff>(f_limit*f_limit);
    if( (vf<0.0)||bOverLimit){
        //dt       = dt * fdec;
        dt       = fmax( dt * fdec, dt_min );
        damping  = damp_max;
        lastNeg  = 0;
        cleanVel();
        //move_GD (dt*scale_dt*scale_dt);
        //for(int i=0; i<n; i++){ vel[i] = kickStart*dt*force[i]; }
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

    if(bOverLimit){
        move_GD (dt_min*sqrt(f_limit*f_limit/ff));
    }else{
        move_LeapFrog( dt*scale_dt );
    }

    //move_LeapFrog( dt*limit_dt_vf2(ff,vv) );
    //move_LeapFrog();
    //move_LeapFrog_vlimit();  // this does not seem to help
    //printf( " %i f v vf  %f %f %f   dt damp  %f %f \n",  stepsDone,   sqrt(ff), sqrt(vv), vf/sqrt(vv*ff),   dt_var, damp_var  );
    //stepsDone++;
    return ff;
}

*/


double DynamicOpt::move_FIRE(){
	double ff=0,vv=0,vf=0;
	//printf( "DEBUG 5.5.1: %i\n", n  );
	for(int i=0; i<n; i++){
		double fi = force[i];
		double vi = vel[i];
		ff += fi*fi;
		vv += vi*vi;
		vf += vi*fi;
        //printf( "move_FIRE %i f %g v %g p %g \n", i, force[i], vel[i], pos[i] );
	}
	//printf( "DEBUG 5.5.2 \n" );
	if( vf < 0.0 ){
	//if( (vf<0.0)||bOverLimit){
		//dt       = dt * fdec;
		dt       = fmax( dt * fdec, dt_min );
	  	damping  = damp_max;
		lastNeg  = 0;
		cleanVel  ( );
		//for(int i=0; i<n; i++){ vel[i] = kickStart*dt*force[i]; }
		//for(int i=0; i<n; i++){ vel[i] = force[i] * 0.5*sqrt(vv/(ff+ff_safety)); }
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
	//printf( "DEBUG 5.5.3 \n" );

	double dt_=dt;

    if( ff>(f_limit*f_limit )){
        dt_*=sqrt(f_limit/sqrt(ff));
        //printf( "force too large: %g => limit dt: %g \n", f, dt_ );
    };


    /*
    // dr_limit
    // dr = (v+f*dt)*dt
    {
        double v=sqrt(vv);
        double f=sqrt(ff+ff_safety);
        double x1,x2;
        if( ((v+f*dt)*dt) > dr_limit*4.0 ){
            quadratic_roots( f+ff_safety, v, -dr_limit, x1, x2 );

            printf( "v %g f %g  dt %g dr %g  | roots %g %g \n", v, f, dt, ((v+f*dt)*dt), x1,x2 );
            dt_=x2;
        }
    }
    */
    move_LeapFrog( dt_ );
	//move_LeapFrog();
	//move_LeapFrog_vlimit();  // this does not seem to help

	//printf( " %i f v vf  %f %f %f   dt damp  %f %f \n",  stepsDone,   sqrt(ff), sqrt(vv), vf/sqrt(vv*ff),   dt_var, damp_var  );
	//stepsDone++;
	return ff;
}

double DynamicOpt::optStep(){
    //cleanForce( );
    getForce( n, pos, force );
    switch( method ){
        //case 0: move_LeapFrog(dt);
        case 0: move_GD(dt);
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
