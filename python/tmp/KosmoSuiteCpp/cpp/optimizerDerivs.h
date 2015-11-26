
typedef void (*FunctionDerivatives)( int n, double * xs, double * dfs );

/*
const double RELAX_damping   = 0.1;
const double RELAX_dt        = 0.1;
const double RELAX_convF2    = 0.000001;
const int    RELAX_maxIters  = 1000;
*/

class OptimizerDerivs{
	public:
		int n;
		int stepsDone=0;
		double * pos;
		double * vel;
		double * force;

		double dt      = 0.2;
		double damping = 0.2;

		FunctionDerivatives func;

	OptimizerDerivs(int n_, double* pos_, double* vel_, double* force_, FunctionDerivatives func_  ){
		n = n_;
		pos   = pos_;
		vel   = vel_;
		force = force_;
		func  = func_;
		for(int i=0; i<n; i++){
			vel  [i]=0;
			force[i]=0;
		}
	};

	double getFmaxAbs( ){
		double fmax = 0;
		for(int i=0; i<n; i++){
			double fi = fabs( force[i] );  
			fmax=(fi>fmax)?fi:fmax;
		}
		return fmax;
	}

	double getFsqSum( ){
		double ff = 0;
		for(int i=0; i<n; i++){
			double fi = force[i];
			ff += fi*fi;
		}
		return ff;
	}

	void autoTimeStep( double dt0, double v0, double f0 ){
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

	virtual void move(){
		double cdamp = 1 - damping;
		for ( int i=0; i<n; i++ ){
			vel[i]  = cdamp*vel[i] + dt*force[i];
			pos[i] += dt*vel[i];
		}
		stepsDone++;
	}

	virtual void step(){
		func( n , pos, force );
		move();
	};
};


class OptimizerFIRE : public OptimizerDerivs {
	public:

		int minLastNeg  = 5;
		double dt       = 0.05;
		double damping  = 0.1;
		double finc     = 1.1; 
		double fdec     = 0.5;
		double falpha   = 0.98;
		double kickStart = 0.0;

		int    lastNeg  = 0; 
		double dt_var   = dt; 
		double damp_var = damping;


	OptimizerFIRE(int n_, double* pos_, double* vel_, double* force_, FunctionDerivatives func_  ):
	OptimizerDerivs  ( n_, pos_, vel_, force_, func_   ){
	};

	virtual void move(){
		double ff=0,vv=0,vf=0;
		for(int i=0; i<n; i++){
			double fi = force[i];
			double vi = vel[i];
			ff += fi*fi;
			vv += vi*vi;
			vf += vi*fi;
		}
		if( vf < 0.0 ){
			dt_var   = dt_var * fdec;
		  	damp_var = damping;
			lastNeg = 0;
			//for(int i=0; i<n; i++){ vel[i] = 0.0d; }
			for(int i=0; i<n; i++){ vel[i] = kickStart*dt_var*force[i]; }
			//for(int i=0; i<n; i++){ vel[i] = dmax*force[i]*sqrt(1/ff)/dt_var; }
		}else{
			//double cf  =     damp_var * sqrt(vv/(ff+1e-8));
			double cf     =     damp_var * sqrt(vv/ff);
			double cv     = 1 - damp_var;
			for(int i=0; i<n; i++){
				vel[i]    = cv * vel[i]  + cf * force[i];
			}
			if( lastNeg > minLastNeg ){
				dt_var    = fmin( dt_var * finc, dt );
				damp_var  = damp_var * falpha;
			}
			lastNeg++;
		}
		for ( int i=0; i<n; i++ ){
			vel[i] += dt_var*force[i];
			pos[i] += dt_var*vel[i];
		}

		//printf( " %i f v vf  %f %f %f   dt damp  %f %f \n",  stepsDone,   sqrt(ff), sqrt(vv), vf/sqrt(vv*ff),   dt_var, damp_var  );
		stepsDone++;
	}

	virtual void step(){
		func( n , pos, force );
		move();
	};

};

/*
class OptimizerFIRE : public OptimizerDerivs {
	public:
		double dt       = 0.1;
		double damping  = 0.1;
		double finc     = 1.1; 
		double fdec     = 0.5;
		double falpha   = 0.99;
		double dt_var   = dt; 
		double damp_var = damping;

	OptimizerFIRE(int n_, double* pos_, double* vel_, double* force_, FunctionDerivatives func_  ):
	OptimizerDerivs  ( n_, pos_, vel_, force_, func_   ){
	};

	virtual void move(){
		double ff=0,vv=0,vf=0;
		for(int i=0; i<n; i++){
			double fi = force[i];
			double vi = vel[i];
			ff += fi*fi;
			vv += vi*vi;
			vf += vi*fi;
		}
		if( vf < 0 ){
			for(int i=0; i<n; i++){ vel[i] = 0.0d; }
			dt_var   = dt_var * fdec;
		  	damp_var = damping;
		}else{
			//double cf  =     damp_var * sqrt(vv/(ff+1e-8));
			double cf  =     damp_var * sqrt(vv/ff);
			double cv  = 1 - damp_var;
			for(int i=0; i<n; i++){
				vel[i]      = cv * vel[i]  + cf * force[i];
			}
			dt_var     = fmin( dt_var * finc, dt );
			damp_var  *= damp_var * falpha;
		}
		for ( int i=0; i<n; i++ ){
			vel[i] += dt_var*force[i];
			pos[i] += dt_var*vel[i];
		}
		stepsDone++;
	}

	virtual void step(){
		func( n , pos, force );
		move();
	};
};
*/



