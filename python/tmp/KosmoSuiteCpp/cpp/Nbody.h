namespace NBody{
	int    n;
	double * mass     = NULL;
	double * invMass  = NULL;

	ODEintegrator_RKF45 solver;

	void forces( double t, int n_doubles, double * Ys, double * dYs ){
		//printf( "DEBUG forces 1 \n" );
		Vec3d *  rs = (  Vec3d * )Ys         ;
		Vec3d *  vs = (( Vec3d * )Ys  ) + n ;
		Vec3d * drs = (  Vec3d * )dYs        ;
		Vec3d * dvs = (( Vec3d * )dYs ) + n ;
		//printf( "DEBUG forces 2 \n" );
		for( int ia = 0; ia<n; ia++ ){
			Vec3d f;    f.set(0.0d);
			Vec3d  ra;  ra.set( rs[ia] );
			double ma = mass[ia];
			for( int ib = 0; ib<n; ib++ ){
				if( ia!=ib )gravity_force( rs[ib] - ra, GRAV_CONTS*ma*mass[ib], f );
			}
			drs[ ia ].set    ( vs[ia]               );
			dvs[ ia ].set_mul( f      , invMass[ia] );
		}
		//printf( " DEBUG dY: " ); for( int i=0; i<n_doubles; i++ ){  printf( " %e ", dYs[i] ); }; printf( " \n" );
		//printf( "DEBUG forces end \n" );
		//exit(0);
	};

	void reallocate( int n_ ){
		n = n_;
		solver.reallocate( 6*n );
		if(    mass != NULL ) delete  mass;
		if( invMass != NULL ) delete invMass;
		   mass = new double[ n ];
		invMass = new double[ n ];
	}

	void setup( int n_, double * mass_, double * poss, double * vs, double * errs ){
		//printf( " Y, dY %i %i \n", NBody::solver.Y, NBody::solver.dY );
		reallocate(n_);
		//printf( " solver allocated to %i %i %i \n", NBody::solver.n, sizeof(NBody::solver.Y), sizeof(NBody::solver.dY) );
		int n3 = 3*n;
		solver.getDerivs = forces;
		for (int i=0; i<n; i++ ){
			int i3 = 3*i;
			mass[i]             =   mass_ [i];
			invMass[i]          = 1/mass_ [i];
			solver.t = 0;
			solver.Y[    i3   ] = poss [i3  ];
			solver.Y[    i3+1 ] = poss [i3+1];
			solver.Y[    i3+2 ] = poss [i3+2];
			solver.Y[ n3+i3   ] = vs   [i3  ];
			solver.Y[ n3+i3+1 ] = vs   [i3+1];
			solver.Y[ n3+i3+2 ] = vs   [i3+2];
			double rErr = 1/errs[ (i<<1)     ];
			double vErr = 1/errs[ (i<<1) + 1 ];
			solver.invMaxYerr[ n3+i3   ] = rErr;
			solver.invMaxYerr[ n3+i3+1 ] = rErr;
			solver.invMaxYerr[ n3+i3+2 ] = rErr;
			solver.invMaxYerr[ n3+i3   ] = vErr;
			solver.invMaxYerr[ n3+i3+1 ] = vErr;
			solver.invMaxYerr[ n3+i3+2 ] = vErr;
			printf( " n,i,m  %i %i %e %e  r %e %e %e  v  %e %e %e \n",   NBody::n, i, mass [i], invMass [i],  poss[i3],poss[i3+1],poss[i3+2],    vs[i3],vs[i3+1],vs[i3+2]   );
		}
	}

	void run( double dt_start, double dt_min, double dt_max, int nstep, double * tsIn, double * tsOut, double * poss, double * vs ){
		int i0 = 0;
		int n3 = 3*n;
		solver.dt_min = dt_min;
		solver.dt_max = dt_max;
		for (int istep=0; istep<nstep; istep++ ){
			//printf( " istep nstep %i %i \n", istep, nstep );
			double tend = tsIn[ istep ];
			solver.integrate_adaptive( dt_start, tend );

			//printf( " .integrate_adaptive DONE \n" );
			tsOut[istep] = solver.t;
			//printf( " DEBUG solver.Y: " ); for( int i=0; i<solver.n; i++ ){  printf( " %e ", solver.Y[i] ); }; printf( " \n" );
			//printf( " store trajectory \n" );
			for (int ia=0; ia<n; ia++ ){
				int i3 = ia*3;
				//printf( " write to index %i %i   %i %i  %i\n", i0,i0+2, i3, i3+n3+2, solver.n );
				poss[ i0   ] = solver.Y[ i3    ];
				poss[ i0+1 ] = solver.Y[ i3 +1 ];
				poss[ i0+2 ] = solver.Y[ i3 +2 ]; i3+=n3;
				vs  [ i0   ] = solver.Y[ i3    ];
				vs  [ i0+1 ] = solver.Y[ i3 +1 ];
				vs  [ i0+2 ] = solver.Y[ i3 +2 ];
				i0 += 3;
			}
			//printf( " store DONE \n" );


		}
		//printf( " run DONE \n" );
	}
}

