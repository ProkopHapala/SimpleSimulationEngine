
class CompresibleMaterial{
	// parameters
	double density0; // [ kg/m3 ] density at zero pressure
	double cv;       // isochoric heat capacity
	//double cp;     // isobaric  heat capacity 
	double kappa;    // heat capacity ratio

	double adiabaticEOS_pressure( compression_ratio ){
		return  pow( compression_ratio, kappa );
	};
}


namespace ShockWaveSpherical{
	// parameters
	int nLayers;
	double p0;
	CompresibleMaterial * materials;

	// axuliary parameters
	double * Vls;     // equlibirum volume of layers
	//double * Mbs;     // mass of between layers
	double * Mps;     // mass of boundary points
	// variables
	// axuliary variables
	double * Areas;  // area     of boundary between layers

	ODEintegrator_RKF45 solver;

/*
	void evalAxuliaryVariables(){

	}
*/


//   m_eff   = delta_p/delta_v0  
//   delta_p = integral(r1,r2){ v(r) m(r) dr} 
//   delta_p = integral(r1,r2){ v0*(r-r1)/(r2-r1) 4*pi*r^2*rho dr} 
//   delta_p = ( ( rho*4*pi*rho*v0 )/( r2-r1 ) ) * (  (r2^4-r1^4)/4 - r2*(r2^3-r1^3)/3 ) ) 
	void effMass_L( double r1, double r2 ){
		double sqR1 = r1*r1;
		double sqR2 = r2*r2;
		return 12.5663706144*( ( sqR2*sqR2 - sqR1*sqR1 )*0.25 - r1*( sqR2*sqR2 - sqR1*sqR1 )*0.3333333333333333  )/( r2 - r1 );
	}
	void effMass_R( double r1, double r2 ){
		double sqR1 = r1*r1;
		double sqR2 = r2*r2;
		return -12.5663706144*( ( sqR2*sqR2 - sqR1*sqR1 )*0.25 - r2*( sqR2*sqR2 - sqR1*sqR1 )*0.3333333333333333  )/( r2 - r1 );
	}

	void forces( double t, int n_doubles, double * Ys, double * dYs ){
		double r1            = 0;
		double r2            = Ys[0];
		double Vsph          = 4.18879020479*r*r*r; 
		double compression   = Vls[0] / Vsph;
		double pressure      = p0          * materials[ 0 ].pressureEOS( compression );
		double rho           = compression * materials[ 0 ].density0; 
		for (int i=1; i<nLayes; i++ ){
			double r3          = Ys [i];
			double compression = Vls[i] / Vsph;
			double mb          = rho * effMass_L( r1, r2 ) + rho_ * effMass_R( r2, r3 );
			double area        = 12.5663706144*r*r;
			       r           = Ys[i];
			double Vsph_       = 4.18879020479*r*r*r;
			double volume      = Vsph_ - Vsph;  
			double pressure_   = p0 * materials[ i+1 ].pressureEOS( volume/Vls[i+1] );

			double accel    =  ( pressure_ - pressure ) * area / mp; 
			dYs[ i+1        ] = Ys[ i+1 + nLayers ];  // d r/dt = v 
			dYs[ i+1+nLayers] = accel;                // d v/dt = a

			pressure = pressure_;
			Vsph     = Vsph_;
		}
	};

	void evalAxuliaryParams(){
		double Vsph_old   = 0;
		for (int i=0; i<nLayes; i++ ){
			double r      = solver.Ys[0];
			double Vsph   = 4.18879020479*r*r*r; 
			Vls[i]        = ( Vsph - Vsoh_old );
		}
	}

	void setup( int nLayers, double ivErr, double irErr, double * Rs, double * Vs, double * materials_ ){
		materials = materials_;
		solver.reallocate( 2*nLayers );
		for( int i=0; i<nLayers; i++ ){
			
		}
		solver.getDerivs = forces;
		for( int i=0; i<nLayers; i++ ){
			materials = (CompresibleMaterial) materials_;
			solver.Y         [ i           ] = Rs[i]; 
			solver.Y         [ i + nLayers ] = Vs[i];
			solver.invMaxYerr[ i           ] = irErr;
			solver.invMaxYerr[ i + nLayers ] = ivErr;
		}
		evalAxuliary();
	}
 
}

