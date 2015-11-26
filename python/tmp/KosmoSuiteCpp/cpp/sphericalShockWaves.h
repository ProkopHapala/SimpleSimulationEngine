

class CompresibleMaterial{
	// parameters
	double density; // [ kg/m3 ] density at zero pressure
	

	double pressureEOS( density, temperature ){

	};

}




namespace NeutronTrace{
	ODEintegrator_RKF45 ;

	// parameters
	int nLayers;
	CompresibleMaterial * materials;
	// axuliary parameters
	double * Ms;     // mass     of boundary between layers
	// variables
	double * vs;     // velocity of boundary between layers
	double * Rs;     // position of boundary between layers
	// axuliary variables
	double * Areas;  // area     of boundary between layers

	void compute_area( ){
		for (int i=0; i<nLayers; i++){
			double R = Rs[i];
			Areas[i] = 12.5663706144*R*R;
		}
	}


	void forces( double t, int n_doubles, double * Ys, double * dYs ){

	};

	void move(){
		
		for(int i=0; i<nLayers i++){
			double e    = Rs[i];
			double m    = Ms[i];
			double area = 12.5663706144*r*r;
		}
		
		
	}
 
}

