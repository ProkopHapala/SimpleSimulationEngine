
// =========== sphere

inline double rayPointDistance2( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& point, double& t ){
	Vec3d pt;	
	pt.set_sub( point, ray0 );
	//printf( " %f %f %f   %f %f %f    %f %f %f \n", point.x, point.y, point.z,   ray0.x, ray0.y, ray0.z,   dPoint.x, dPoint.y, dPoint.z    );
	t  = pt.dot( hRay );
    pt.add_mul( hRay, -t );
	return pt.norm2();
}

inline double raySphere( const Vec3d& ray0, const Vec3d& hRay, double R, const Vec3d& center ){
	double t0;
	double rt2 = rayPointDistance2( ray0, hRay, center, t0 );
	double dr2 = R*R - rt2;
	//printf( " %f %f %f %f \n", t0, dr2, rt2, R*R );
	if( dr2 > 0 ){
		return t0 - sqrt( dr2 );
		//return t0 + sqrt( dr2 ); // this would be the second branch
	}else{
		return INFINITY;
	}
}

inline void sphereNormal( double t, const Vec3d& ray0, const Vec3d& hRay, const Vec3d& center, Vec3d& normal ){
	normal.set_sub( ray0, center );
	normal.add_mul( hRay, t );
}


