

inline int find_val( int imin, int imax, double val, double * vals ){
	if( (vals[imin]>val)||(vals[imax]<val) ) return -1;
	while( true ){
		if( imin == imax){	return imin; };
		int i = imin + (imax-imin)>>1;
		if   ( vals[i]<val ){ imin = i; }
		else                { imax = i; }
	}	
}

inline int find_val_after( int imin, int n, double val, double * vals ){
	for( int i=imin; i<n; i++ ){
		if( vals[i]>val ) return i-1;
	}
	return -1;
}

class SpaceShip : public SpaceBody {
	public:

	int nPlan;
	int iPlan_last;
	double * plan_t;
	Vec3d  * plan_propAcc;

	bool timeSeqMode = false; // if reading from plan sequntially 

	inline void addPlanedThrustTo( int i, double dh, Vec3d& thrust ){
		thrust.add_mul( plan_propAcc[i  ], (1-dh) );
		thrust.add_mul( plan_propAcc[i+1],    dh  );
	}

	inline void addPlanedThrustTo( double t, Vec3d& thrust ){
		int i = find_val( 0, nPlan-1, t, plan_t );
		double dh = ( t - plan_t[i] )/( plan_t[i+1] - plan_t[i] );
		addPlanedThrustTo( i, dh, thrust );
	}

	virtual void evalForce( double t ) {
		SpaceBody::evalForce();
		int i;
		if   ( timeSeqMode ){  i = find_val_after( iPlan_last, nPlan,   t, plan_t ); }
		else                {  i = find_val      (          0, nPlan-1, t, plan_t ); }
		double dh = ( t - plan_t[i] )/( plan_t[i+1] - plan_t[i] );
		addPlanedThrustTo( i, dh, accel );
	}

};



