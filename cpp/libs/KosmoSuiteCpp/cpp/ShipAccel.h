

namespace ShipAccel{
	double GM;          // central body gravitation  GRAV_CONTS * Mass
	double accel;

	ODEintegrator_RKF45 solver;

	void forces( double t, int n_doubles, double * Ys, double * dYs ){
		Vec3d f,v;
		f.set(0.0);	
		gravity_force       ( ((Vec3d*)Ys)[0], GM, f  );
		v.set( ((Vec3d*)Ys)[1] );
		f.add_mul( v, accel / sqrt( v.norm2() + 1e-32 ) );  // acceleartion along velocity direction
		((Vec3d*)dYs)[0].set( v );
		((Vec3d*)dYs)[1].set( f );
	};

	void setup( double GM_, double accel_, double ivErr, double irErr, double * pos, double * vs ){
		solver.reallocate( 6 );
		GM = GM_;
		accel = accel_;
		solver.getDerivs = forces;
		solver.Y[ 0 ]    = pos [0];
		solver.Y[ 1 ]    = pos [1];
		solver.Y[ 2 ]    = pos [2];
		solver.Y[ 3 ]    = vs   [0];
		solver.Y[ 4 ]    = vs   [1];
		solver.Y[ 5 ]    = vs   [2];
		solver.invMaxYerr[ 0 ] = irErr;
		solver.invMaxYerr[ 1 ] = irErr;
		solver.invMaxYerr[ 2 ] = irErr;
		solver.invMaxYerr[ 3 ] = ivErr;
		solver.invMaxYerr[ 4 ] = ivErr;
		solver.invMaxYerr[ 5 ] = ivErr;
	}

	void run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
		//printf( " HERE  ShipAccel::run  HERE \n" );
		int i0 = 0;
		solver.t = 0;
		solver.dt_min = dt_min;
		solver.dt_max = dt_max;
		//printf( " nstep %i dt_start,dt_min,dt_max %e %e %e  \n",  nstep,   dt_start, dt_min, dt_max );
		for (int istep=0; istep<nstep; istep++ ){
			double tend = tsIn[ istep ];
			solver.integrate_adaptive( dt_start, tend );
			//printf( " DEBUG i %i  dt %e Y: ", istep, solver.dt_adapt ); for( int i=0; i<solver.n; i++ ){  printf( " %e ", solver.Y[i] ); }; printf( " \n" );
			tsOut[istep]  = solver.t;
			poss [ i0   ] = solver.Y[ 0 ];
			poss [ i0+1 ] = solver.Y[ 1 ];
			poss [ i0+2 ] = solver.Y[ 2 ];
			vs   [ i0   ] = solver.Y[ 3 ];
			vs   [ i0+1 ] = solver.Y[ 4 ];
			vs   [ i0+2 ] = solver.Y[ 5 ];
			i0+=3;
		}
	}
}

