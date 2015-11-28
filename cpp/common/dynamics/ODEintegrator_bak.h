
typedef void (*ODEderivFunc)( double t, int n, double * Ys, double * dYs );

// =================================================
// ============ ODEintegrator        ===============
// =================================================

class ODEintegrator{
	public:
	int n;
	double    t;
	double *  Y;
	double * dY;

	int    MAX_STEPS = 10000;

	ODEderivFunc getDerivs;

	virtual void step( double dt ){
		getDerivs( t, n, Y, dY );
		for (int i=0; i<n; i++){
			Y[i] += dY[i] * dt;
		}
		t+=dt;
	}

	virtual int integrate( double dt, double Tmax ){
		int istep=0;
		while ( t<Tmax ){	
			step( dt ); 
			istep++;
			if ( istep> MAX_STEPS ){
				printf(" ERROR: MAX_STEPS ( %i ) achieved,  \n", istep );
				break;
			}
		}
		return istep;
	}

};

// =================================================
// ============ ODEintegrator_RKF45  ===============
// =================================================

class ODEintegrator_RKF45 : public ODEintegrator{

	// constants
	const double at  [6] = {         0,  1/4.0d,         3/8.0d,       12/13.0d,        1,  1/2.0d   };
	const double cs1 [6] = { 16/135.0d,       0,  6656/12825.0d, 28561/56430.0d, -9/50.0d,  2/55.0d  };
	const double cs2 [5] = { 25/216.0d,       0,   1408/2565.0d,   2197/4104.0d,  -1/5.0d            };
	const double 
	b10 =  1/4.0d,        
	b20 =  3/32.0d,      b21 = 9/32.0d,       
	b30 =  1932/2197.0d, b31 = -7200/2197.0d, b32 = 7296/2197.0d,
	b40 =  439/216.0d,   b41 = -8,            b42 = 3680/513.0d,   b43 = -845/4104.0d,
	b50 = -8/27.0d,      b51 = 2,             b52 = -3544/2565.0d, b53 = 1859/4104.0d,  b54 = -11/40.0d; 

	// axuliary state variables
	double *invMaxYerr;
	double *Ynew, *Yerr;
	double *dY0, *dY1, *dY2, *dY3, *dY4, *dY5;

	// adaptive timestep parameters
	double SAFETY    =  0.2d; 
	double PGROW     =  1.2d;
	double PSHRINK   =  0.7d;
	int    MAX_ADAPT =  10; 

	double dt_min, dt_max;
	double error;
	double dt_adapt;

	// ==================   step =======================

	virtual void step( double dt ){
		// predictor step ( get all derivatives )
		getDerivs( t + at[0]*dt, n, Y   , dY0 );  for (int i = 0; i<n; i++) {  Ynew[i] = Y[i] + dt *   b10*dY1[i];                                                        } 
		getDerivs( t + at[1]*dt, n, Ynew, dY1 );  for (int i = 0; i<n; i++) {  Ynew[i] = Y[i] + dt * ( b20*dY0[i] + b21*dY1[i]                                         ); } 
		getDerivs( t + at[2]*dt, n, Ynew, dY2 );  for (int i = 0; i<n; i++) {  Ynew[i] = Y[i] + dt * ( b30*dY0[i] + b31*dY1[i] + b32*dY2[i]                            ); } 
		getDerivs( t + at[3]*dt, n, Ynew, dY3 );  for (int i = 0; i<n; i++) {  Ynew[i] = Y[i] + dt * ( b40*dY0[i] + b41*dY1[i] + b42*dY2[i] + b43*dY3[i]               ); } 
		getDerivs( t + at[4]*dt, n, Ynew, dY4 );  for (int i = 0; i<n; i++) {  Ynew[i] = Y[i] + dt * ( b50*dY0[i] + b51*dY1[i] + b52*dY2[i] + b53*dY3[i] + b54*dY4[i]  ); } 
		getDerivs( t + at[5]*dt, n, Ynew, dY5 );  

		 // construct solutions and error estimator
		for (int i = 0; i < n; i++) {
			double dYdt1 = cs1[0]*dY0[i] + cs1[2]*dY2[i] + cs1[3]*dY3[i] + cs1[4]*dY4[i] + cs1[5]*dY5[i];
			double dYdt2 = cs2[0]*dY0[i] + cs2[2]*dY2[i] + cs2[3]*dY3[i] + cs2[4]*dY4[i];
			Ynew[i] = Y[i] + dt * dYdt1;
			Yerr[i] = dt * ( dYdt2 - dYdt1 );
		}
		t +=dt;
	}

	// ==================   adaptive step   =======================

	virtual void adaptive_step( ){
		for (int iadapt=0; iadapt<MAX_ADAPT; iadapt++ ){
			double tbak = t;
			step( dt_adapt ); 
			error = 0;	for (int i=0; i<n; i++ ){ error = fmax( error, abs( Yerr[i]*invMaxYerr[i] ) ); }
			if ( error < 1 ){
				for (int i=0; i<n; i++ ){ Y[i] = Ynew[i]; }
				if ( error < SAFETY ){	
					dt_adapt = fmin( dt_adapt * PGROW , dt_max );
				}
				break;
			}else{
				t = tbak;
			}
			dt_adapt = fmax( dt_adapt * PSHRINK, dt_min );
		}

	}

	// ============= integrate =======  

	virtual int integrate( double dt, double Tmax ){
		int istep=0;
		while ( t<Tmax ){	
			step( dt ); 
			istep++;
			if ( istep> MAX_STEPS ){
				printf(" ERROR in integrate : MAX_STEPS ( %i ) achieved,  \n", istep );
				break;
			}
		}
		return istep;
	}

	// ============= integrate_adaptive =======  

	virtual int integrate_adaptive(  double dt_start, double Tmax ){	
		int istep=0;
		dt_adapt = dt_start;
		while ( t < Tmax ){
			adaptive_step( );
			if ( istep> MAX_STEPS ){
				printf(" ERROR in integrate_adaptive : MAX_STEPS ( %i ) achieved, dt = %e t= %e error= %e \n", istep, dt_adapt, t, error );
				break;
			}
			istep++;
		}
		return istep;
	}  

};
