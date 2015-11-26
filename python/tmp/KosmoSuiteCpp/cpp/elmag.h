

namespace ElMag{
	double charge;
	double mass;
	double QinvM;

	int ncoil;
	Vec3d * coil_ps;
	Vec3d * coil_dIs;

	ODEintegrator_RKF45 solver;

	void field_at( const Vec3d& where, Vec3d& B, int n, Vec3d * ps, Vec3d * dIs ){
		for (int i=0; i<n; i++){ biot_savart_element( where - ps[i], dIs[i], B ); }
	}

	void sample_field ( int m, Vec3d * wheres, Vec3d * Bs, int n, Vec3d * ps, Vec3d * dIs ){
		for (int i=0; i<m; i++){ Vec3d B; B.set(0.0); field_at( wheres[i], B, n, ps, dIs );  Bs[i].set(B); }
	}

	void forces( double t, int n_doubles, double * Ys, double * dYs ){
		Vec3d pos, v, f, B;
		B.set(0.0); f.set(0.0);
		pos.set( ((Vec3d*)Ys)[0] );
		v  .set( ((Vec3d*)Ys)[1] );
		field_at( pos, B, ncoil, coil_ps, coil_dIs );
		lorenz_force( QinvM, v, B, f );
		((Vec3d*)dYs)[0].set( v );
		((Vec3d*)dYs)[1].set( f );
	};

	void setup( double charge_, double mass_, double ivErr, double irErr, double * pos, double * vs, int ncoil_, Vec3d* coil_ps_, Vec3d* coil_dIs_ ){
		solver.reallocate( 6 );
		mass     = mass_;
		charge   = charge_;  
		QinvM    = charge/mass;
		ncoil    = ncoil_;
		coil_ps  = coil_ps_;
		coil_dIs = coil_dIs_;
		solver.getDerivs = forces;
		solver.Y[ 0 ]    = pos[0];
		solver.Y[ 1 ]    = pos[1];
		solver.Y[ 2 ]    = pos[2];
		solver.Y[ 3 ]    = vs [0];
		solver.Y[ 4 ]    = vs [1];
		solver.Y[ 5 ]    = vs [2];
		solver.invMaxYerr[ 0 ] = irErr;
		solver.invMaxYerr[ 1 ] = irErr;
		solver.invMaxYerr[ 2 ] = irErr;
		solver.invMaxYerr[ 3 ] = ivErr;
		solver.invMaxYerr[ 4 ] = ivErr;
		solver.invMaxYerr[ 5 ] = ivErr;
	}

	void run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
		int i0 = 0;
		solver.t = 0;
		solver.dt_min = dt_min;
		solver.dt_max = dt_max;
		//printf( " nstep %i dt_start,dt_min,dt_max %e %e %e  \n",  nstep, dt_start, dt_min, dt_max );
		for (int istep=0; istep<nstep; istep++ ){
			double tend = tsIn[ istep ];
			solver.integrate_adaptive( dt_start, tend );
			//printf( " DEBUG i %i t %e dt %e Y: ", istep, solver.t, solver.dt_adapt ); for( int i=0; i<solver.n; i++ ){  printf( " %e ", solver.Y[i] ); }; printf( " \n" );
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

