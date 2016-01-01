
#ifndef  raytrace_h
#define  raytrace_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>


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



// =========== Plane

inline double rayPlane( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& normal, const Vec3d& point ){
	double nh = normal.dot( hRay  );
	double np = normal.dot( point );
	double n0 = normal.dot( ray0  );
	//return ( n0 + np ) / nh;
	return ( np - n0 ) / nh;
}

// =========== Triangle

inline bool pointInTriangleEdges( const Vec3d& pa, const Vec3d& pb, const Vec3d& ab, const Vec3d& bc, const Vec3d& ca ){
	Vec3d nab,nbc,nca; 
	nab.set_cross( pa, ab );
	nbc.set_cross( pb, bc );
	nca.set_cross( pa, ca );
	double sgn1 = nca.dot( nab );
	double sgn2 = nca.dot( nbc );
	printf( " sgn1, sgn2 %f, %f \n", sgn1, sgn2 );
	return ( sgn1>0 )&&( sgn2>0 );
}


inline double rayTriangle(
	const Vec3d &ray0, const Vec3d &hRay,
	const Vec3d &a,    const Vec3d &b,    const Vec3d &c,
	bool& inside, Vec3d& p
){
	Vec3d ab,bc,ca, normal;
	ab.set_sub( b, a );
	bc.set_sub( c, b );
	ca.set_sub( a, c );

	normal.set_cross( ab, bc );
	
	double t = rayPlane( ray0, hRay, normal, a );
	p.set( ray0 ); p.add_mul( hRay, t );

	inside = pointInTriangleEdges( p-a, p-b, ab, bc, ca );

	return t;
} 


#endif







/*

//  http://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
bool rayTriangleIntersect(
	const Vec3d &orig, const Vec3d &dir,
	const Vec3d &a,    const Vec3d &b,    const Vec3d &c,
	double &t
){
	// compute plane's normal
	Vec3d ab = b - a;
	Vec3d ac = c - a;
	Vec3d cb = c - b;
	// no need to normalize
	Vec3d N; N.set_cross( ab, ac ); // N
	double area2 = N.norm();

	// Step 1: finding P

	// check if ray and plane are parallel ?
	double NdotRayDirection = N.dot( dir );
	//if (fabs(NdotRayDirection) < kEpsilon) // almost 0
	//return false; // they are parallel so they don't intersect !

	// compute d parameter using equation 2
	double d = N.dot( a );

	// compute t (equation 3)
	t = ( N.dot(orig) + d) / NdotRayDirection;

	
	//if (t < 0) return false;        // check if the triangle is in behind the ray

	
	Vec3d P = orig + t * dir;       // compute the intersection point using equation 1

	// Step 2: inside-outside test
	Vec3d C;                          // vector perpendicular to triangle's plane

	// edge 0
	Vec3d edge0 = v1 - v0;
	Vec3d vp0   = P - v0;
	C.set_cross( edge0, vp0 );
	if (N.dotProduct(C) < 0) return false; // P is on the right side

	// edge 1
	Vec3d edge1 = v2 - v1;
	Vec3d vp1 = P - v1;
	C = edge1.crossProduct(vp1);
	if (N.dotProduct(C) < 0) return false; // P is on the right side

	// edge 2
	Vec3d edge2 = v0 - v2;
	Vec3d vp2 = P - v2;
	C = edge2.crossProduct(vp2);
	if (N.dotProduct(C) < 0) return false; // P is on the right side;

	return true; // this ray hits the triangle
} 

*/
