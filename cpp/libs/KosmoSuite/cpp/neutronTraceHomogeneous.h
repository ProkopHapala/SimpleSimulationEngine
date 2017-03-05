

class NuclearMaterial{
	// parameters
	double spontalneous_fission; // [ Hz    ] number of spontaneous fissions per seconds
	double neutron_mult;         // [ 1     ] number of neutrons per fission
	double S_fission;            // [ m^2   ] fission crossection
	double S_capture;            // [ m^2   ] fission corssection
	double Ndens;                // [ 1/m^3 ] particle density

	// axuliary
	double S_total       ; 
	double fission_share ; 
	double capture_share ; 
	double Ndens         ;
	double mean_free_path;

	eval_axuliary( ){
		S_total        = S_fission + S_capture;
		fission_share  = S_fission / S_total;
		capture_share  = S_capture / S_total;
		Ndens          = 6.022140857e+23 * density / MolarMass; 
		mean_free_path = 1 / ( Ndens * S_total );
	}
}



namespace NeutronTrace{
	NuclearMaterial material;

	int nbind_E,nbins_R,nbins_tot;

	double * histogram;

	double getCollisonDist( double mean_free_path,  ){
		double rnf =   randf();
		double lp  = - log(rnf) * mean_free_path; 
		return lp;
	}

	void store( double R, double lnE ){
		int ibin_R  = (int)( R   / dRbin );
		int ibin_E  = (int)( lnE / dlnEbin );
		histogram[ ibin_E*nbin_E + ibin_R ] ++;    
	} 

	void traceNeutron( const Vec3d& pos0, const Vec3d& hray0, double velocity0 ){
		Vec3d hray,pos;
		hray.set(hray0);
		pos .set(pos0);
		double velocity;
		while( true ){
			double l = getCollisonDist( material.mean_free_path );
			pos.add_mul( hray, lp );
			
			velocity *= fslow; 
		}
	}
 
}

