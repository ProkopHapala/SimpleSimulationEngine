

class NuclearMaterial{
	// parameters
	double spontalneous_fission; // [ Hz    ] number of spontaneous fissions per seconds
	double neutron_mult;        // [ 1     ] number of neutrons per fission
	double S_fission;           // [ m^2   ] fission crossection
	double S_capture;           // [ m^2   ] fission corssection
	double Ndens;               // [ 1/m^3 ] particle density

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
	int nLayers;
	double * Rs;
	NuclearMaterial * materials;

	double * Rdens_fission;
	double * Tdens_fission;

	double * ls;    //  length of ray between outer shell intersections
	double * mfps;  //  dimension-less distance within each layer

	// segments
	double  *  dlds; // dimension les distances
	double  *    ds; // distances
	int     * imats; // material index
	


	double traceHomogeneousSegment( double length, const NuclearMaterial& mat, Vec3d& hRay, Vec3d& pos ){
		double rnf =   randf();
		double lp  = - log(rnf) * mean_free_path; 
		if ( lp < length ){

		}else{

		}
	}

	void traceNeutron( ){
		Vec3d hray,pos;
		double velocity;
		
		while( true ){
			double t0;
			double R2  = pos.norm2(); 
			double rt2 = rayPointDistance2( pos, hRay, SphereCenter, t0 );	
			int iin  = 0;
			int imin = 0;
			for( imin=nLayers-1; i>=0; i-- ){   // find intersection with layers boundary	
				double Ri   = Rs[imin]; 
				double Ri2  = Ri*R;
				if( Ri2 < R2 ) iin = i;
				double li2 = Ri2 - rt2;
				if( li2 > 0 ){
					double li = sqrt( dr2 );
					ls[imin] = li;
				}else{
					break;
				};
			}
			double mfpsum = 0;
			int nseg = 1;
			// build segments
			if( t0 > 0 ){ // toward ray
				double lin  = ls[iin+1] - t0;
				dlds[nseg] = (ls[i] - lin ) / materials[ iin ].mean_free_path; 
				for( int i=iin; i<nLayers; i++ ){	
					dlds[nseg] = (ls[i] - ls[i-1]) / materials [ i ].mean_free_path; 
					nseg++;
				}
			} else { // outward ray
				for( int i=iin; i>imin; i-- ){
				
				}
				for( int i=imin; i<nLayers; i++; ){

				}
			}
			// pass segments
			double rnf   =   randf();
			double dld   = - log(rnf); 
			double dlsum = 0; 
			for ( int iseg=0; iseg<nseg; iseg++ ){
				dlsum += dlds[i];
				if( dlsum > dld ){ // hit 
					double t = materials[ imats[ iseg ] ].mean_free_path * ( dlsum - dld ); 
					pos.add_mul( hRay, t1 + t2 +  );
					break;					
				}
			}
		}
	}
 
}

