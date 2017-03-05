

// according to http://nuclearweaponarchive.org/Nwfaq/Nfaq4-1.html#Nfaq4.1


namespace FissionPulse{

	// parameters

	double mass;       // inertial mass of the whole assembly
	double alphaDOF;   // alpha of degrees of freedom  3/2 for ideal gas  

	double total_crossection;  // [ m^2 ]
	double fission_share;      // [ 1   ]
	double back_scatter_share; // [ 1   ]

	double generation_rate;    // [ Hz ] inverse average generation time 
	double spontaneous_rate;   // [ Hz ] spontaneous fissions per second

	double E_fission;          // [ J ] energy per fission reaction
	double nmult;              // [ 1 ] neutrons per induced fission
	double nmult_spontal;      // [ 1 ] neutrons per spontal fission

	double time_trigger = INFINITY;
	double R_trigger    = 0;
	double Nn_spark     = 0;

	// axuliary variables
	double force;             //  [ N ] 
	double interaction_prob;  //  [ 1 ]

	ODEintegrator_RKF45 solver;

	void forces( double t, int n_doubles, double * Ys, double * dYs ){

		double R  = Ys[ 0 ];
		double v  = Ys[ 1 ];
		double Nf = Ys[ 2 ];
		double Nn = Ys[ 3 ];
		double Q  = Ys[ 4 ];

		//double U        = U0  +  Q  - mass*v*v;
		//double V        = 1.33333333333 * M_PI * R*R*R;
		//double pressure = U / ( alphaDOF * V );
		//double force    = 4*M_PI*R*R*pressure; 
		double W          = 0.5 * mass * v * v;
		double U          = Q  - W ;
		       U          = ( U > 0 ) ? U : 0;
		       force      = 3 * U / ( alphaDOF * R );

		//double interactions = Nn * exp( - total_crossection * Nf / (R*R) )    * generation_rate;

		       interaction_prob = total_crossection * Nf / (R*R);
		double interactions     = Nn * generation_rate * interaction_prob;
		double fissions         = fission_share    * interactions ; 
		double spontal          = spontaneous_rate * Nf;

		double dNf =                  -spontal -              fissions                                                                  ;
		double dNn =  nmult_spontal *  spontal +      nmult * fissions     + back_scatter_share * interactions  -  generation_rate *Nn  ;
		double dQ  =  E_fission     *( spontal +              fissions                                                                 );                                                     
		double dv  = force / mass;

 		double dR  = v;

		//printf( "  Ys: %e %e %e %e %e   aux:  %e %e %e %e %e %e   dYs:  %e %e %e %e %e \n",  R, v, Nf, Nn, Q,   W, U, force, interactions, fissions, spontal,    dR, dv, dNf, dNn, dQ );

		dYs[ 0 ] = 0;
		dYs[ 1 ] = 0;
		dYs[ 2 ] = 0;
		dYs[ 3 ] = 0;
		dYs[ 4 ] = 0;

		dYs[ 0 ] = dR;
		dYs[ 1 ] = dv;
		dYs[ 2 ] = dNf;
		dYs[ 3 ] = dNn;
		dYs[ 4 ] = dQ;

	};

	void set_params( 
		double mass_,				// inertial mass of the whole assembly
		double alphaDOF_,			// alpha of degrees of freedom  3/2 for ideal gas  
		double total_crossection_,  // [ m^2 ]
		double fission_share_,      // [ 1   ]
		double back_scatter_share_, // [ 1   ]
		double generation_rate_,    // [ Hz ] inverse average generation time 
		double spontaneous_rate_,   // [ Hz ] spontaneous fissions per second	
		double E_fission_,          // [ J ] energy per fission reaction
		double nmult_,              // [ 1 ] neutrons per induced fission
		double nmult_spontal_       // [ 1 ] neutrons per spontal fission
	){
		solver.reallocate( 5 );
		solver.getDerivs = forces;
		mass               = mass_;
		alphaDOF           = alphaDOF_;
		total_crossection  = total_crossection_; 
		fission_share      = fission_share_;     
		back_scatter_share = back_scatter_share_;
		generation_rate    = generation_rate_;
		spontaneous_rate   = spontaneous_rate_; 	
		E_fission          = E_fission_;          
		nmult              = nmult_;            
		nmult_spontal      = nmult_spontal_;         
	}

	void set_trigger( double time_trigger_, double R_trigger_, double Nn_spark_ ){  
		time_trigger = time_trigger_; R_trigger = R_trigger_; Nn_spark = Nn_spark_;
	}

	void set_initial( double R0, double v0, double Nf0, double Nn0, double Q0 ){
		solver.Y[ 0 ]    = R0;
		solver.Y[ 1 ]    = v0;
		solver.Y[ 2 ]    = Nf0;
		solver.Y[ 3 ]    = Nn0;
		solver.Y[ 4 ]    = Q0 + 0.5 * mass * v0 * v0;
		solver.t         = 0;
		//printf( " %e %e %e %e %e \n", R0, v0, Nf0, 0.0,  Q0 );
		//printf( " %e %e %e %e %e \n", solver.Y[ 0 ], solver.Y[ 1 ], solver.Y[ 2 ], solver.Y[ 3 ], solver.Y[ 4 ] );
	}

	void run_fixStep( double dt, int nstep, double * buff ){
		printf( " %i %f \n", nstep, dt );
		int ibuff = 0;
		for (int istep=0; istep<nstep; istep++ ){
			//if( (t<t_trigger)&&((t+dt)>t_trigger) ){ solver.Y[3] += Nn_spark; };
			solver.step( dt );
			for (int i=0; i<solver.n; i++ ){ solver.Y[i] = solver.Ynew[i]; }
			//printf( " %i %e %e %e %e %e \n", istep, solver.Y[ 0 ], solver.Y[ 1 ], solver.Y[ 2 ], solver.Y[ 3 ], solver.Y[ 4 ] );
			buff [ ibuff   ] = solver.Y[ 0 ];
			buff [ ibuff+1 ] = solver.Y[ 1 ];
			buff [ ibuff+2 ] = solver.Y[ 2 ];
			buff [ ibuff+3 ] = solver.Y[ 3 ];
			buff [ ibuff+4 ] = solver.Y[ 4 ];

			buff [ ibuff+5 ] = force;
			buff [ ibuff+6 ] = interaction_prob;
			ibuff += 7;
		}
	}
}

