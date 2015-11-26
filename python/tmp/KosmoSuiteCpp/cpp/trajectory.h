inline int search_between( double t, int imin, int imax, double * ts ){
	int di = imax - imin;
	while ( di > 1 ){
		int i = (imin + imax) >> 1;
		if   ( t > ts[i] ){ imin=i; }
		else              { imax=i; }
		di = imax - imin;
	}
	return imin;
}

inline int search_behind( double t, int imin, double * ts, int n ){
	int di   = 1;
	int imax = imin+di; imax = (imax>n)?n:imax;
	while ( ts[imax] < t ){
		di   = di << 1;
		imax = imin+di; imax = (imax>n)?n:imax;
	}
	return search_between( t, imin, imax );
}

inline int search_before( double t, int imax, double * ts ){
	int di   = 1;
	int imin = imin-di; imin = (imin<0)?0:imin;
	while ( ts[imin] < t ){
		di   = di << 1;
		imin = imin-di; imin = (imin<0)?0:imin;
	}
	return search_between( t, imin, imax );
}


class Trajectory{
	public:
	int n;       // number of times
	int m;       // number of numbers per time 
	int ilast;
	double * ts;
	double * xs;

	void setup( int n, int m, ts ){

	}

	inline int search     ( double t ){	search_between( t, 0, n,  ts    ); }
	inline int search_next( double t ){	search_behind ( t, ilast, ts, n ); }
	inline int search_prev( double t ){	search_before ( t, ilast, ts    ); }

	inline double get_d( double t, int it ){
		double t0 = ts[ it   ];
		double t1 = ts[ it+1 ];
		return d  = (t-t0)/(t1-t0);
	}

	inline void interpolate_rv( double d, int it, Vec3d& r, Vec3d& v ){
		double m  = 1 - d;  
		Vec3d r0,r1,v0,v1;
		double xstmp = xs+(it*m);
		r0.set( ((Vec3d*)xstmp)[0] );	
		v0.set( ((Vec3d*)xstmp)[1] );
		xstmp += m;			
		r1.set( ((Vec3d*)xstmp)[0] );	
		v1.set( ((Vec3d*)xstmp)[1] );
		r.set_lincomb( m, r0, d, r1 );
		v.set_lincomb( m, v0, d, v1 );
	}

	inline void get_interpolated_rv( double t, Vec3d& r, Vec3d& v ){
		int it = search( t );
		int d  = get_d( t, it );	
		interpolate_rv( d, it, r, v )
	}

	inline void get_interpolated_rvs( double t, int m, double * tms, Vec3d * rs, Vec3d * vs ){
		int it = 0;
		for( int i=0; i<m; i++ ){
			double t = tms[i]; 
			int it   = search_behind( t, it, ts, n );
			int d    = get_d( t, it );	
			interpolate_rv( d, it, rs[i], vs[i] );
		}
	}

}
