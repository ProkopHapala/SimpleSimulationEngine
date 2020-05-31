
#ifndef  raytrace_h
#define  raytrace_h

/*

Quadric intersection:
    http://skuld.bmsc.washington.edu/people/merritt/graphics/quadrics.html

ToDo:
-Numerical solution of any intersection


*/


#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "Vec2.h"
#include "Vec3.h"

// =========== sphere

const double t_inf = 1e+300;

inline void setToRay( const Vec3d& ro, const Vec3d& rd, Vec3d& p ){
    Vec3d d; d.set_sub(p,ro); d.makeOrthoU(rd);
    //printf( "setToRay %g %g %g \n", d.x, d.y, d.z );
    p.sub(d);
}

inline double rayPointDistance2( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& point, double& t ){
	Vec3d pt;
	pt.set_sub( point, ray0 );
	//printf( " (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", point.x, point.y, point.z,   ray0.x, ray0.y, ray0.z,   pt.x, pt.y, pt.z    );
	t  = hRay.dot( pt );
	//printf( " %f (%3.3f,%3.3f,%3.3f) \n", t, pt.x, pt.y, pt.z    );
    pt.add_mul( hRay, -t );
    //printf( " %f %f (%3.3f,%3.3f,%3.3f) \n", t, pt.norm2(), pt.x, pt.y, pt.z    );
	return pt.norm2();
}

inline double linePointDistance2( const Vec3d& p0, const Vec3d& hdir, const Vec3d& point, double tmax ){
	Vec3d d;
	d.set_sub( point, p0 );
	double t  = hdir.dot( d );
	if      ( t>tmax ){ d.add_mul( hdir, -tmax ); }
	else if ( t>0    ){ d.add_mul( hdir, -t    ); }
	return d.norm2();
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

// =========== Line

inline double rayLineDist( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& p0, const Vec3d& hL ){
    // http://mathworld.wolfram.com/Line-LineDistance.html
    Vec3d rxl; rxl.set_cross(hRay, hL);
    double r = rxl.norm();
    return rxl.dot( p0-ray0 )/r;
}

inline double rayLine( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& p0, const Vec3d& hL, double& t1, double& t2 ){
    // https://math.stackexchange.com/questions/2083573/mutually-closest-point-of-two-lines-defined-by-vectors
    Vec3d d;     d.set_sub(p0,ray0);    // (p2-p1)
    Vec3d rxl; rxl.set_cross(hRay, hL); // normal vector to (d1,d2) plane
    rxl.normalize();
    double c = rxl.dot(d);
    d.add_mul(rxl,c); // (p2-p1) in d1 d2 plane
    // --- solve system of linear equtions    (p2-p1)_T = d1*t - d2*s;
    double crl  = hRay.dot(hL); // <d1|d2>
    double c1   = d.dot(hRay);
    double c2   = d.dot(hL);
    t1 =   (c1-c2*crl)/(1-crl*crl);
    t2 =  -(c2-c1*crl)/(1-crl*crl);
    return fabs(c); // distance
    //return (c1-c2*crl)/(1-crl*crl);
    //return
}

// =========== Plane

inline double rayPlane( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& normal, const Vec3d& point ){
	double nh = normal.dot( hRay  );
	double np = normal.dot( point );
	double n0 = normal.dot( ray0  );
	//printf( "rayPlane %g | %g %g %g\n", ( np - n0 ) / nh, nh, np, n0 );
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
	//printf( " sgn1, sgn2 %f, %f \n", sgn1, sgn2 );
	return ( sgn1>0 )&&( sgn2>0 );
}


inline double rayTriangle(
	const Vec3d &ray0, const Vec3d &hRay,
	const Vec3d &a,    const Vec3d &b,    const Vec3d &c,
	bool& inside, Vec3d& hitpos, Vec3d& normal
){
	Vec3d ab,bc,ca;
	ab.set_sub( b, a );
	bc.set_sub( c, b );
	ca.set_sub( a, c );

	normal.set_cross( ab, bc );

	double t = rayPlane( ray0, hRay, normal, a );
	hitpos.set( ray0 ); hitpos.add_mul( hRay, t );

	inside = pointInTriangleEdges( hitpos-a, hitpos-b, ab, bc, ca );

	return t;
}

inline bool rayInTriangle( const Vec3d& a_, const Vec3d& b_, const Vec3d& c_, const Vec3d& hX, const Vec3d& hY ){
/*
	Vec2d a; a.set( hX.dot(a_), hY.dot(a_)  );
	Vec2d b; b.set( hX.dot(b_), hY.dot(b_)  );
	if( 0 > ( a.x*(b.y-a.y) - a.y*(b.x-a.x) ) ) return false;
	Vec2d c; c.set( hX.dot(c_), hY.dot(c_)  );
	if( 0 > ( b.x*(c.y-b.y) - b.y*(c.x-b.x) ) ) return false;
	if( 0 > ( c.x*(a.y-c.y) - c.y*(a.x-c.x) ) ) return false;
*/

	Vec2d a; a.set( hX.dot(a_), hY.dot(a_)  );
	Vec2d b; b.set( hX.dot(b_), hY.dot(b_)  );
	double sgn = a.x*(b.y-a.y) - a.y*(b.x-a.x);
	Vec2d c; c.set( hX.dot(c_), hY.dot(c_)  );
	if( 0 > sgn*( b.x*(c.y-b.y) - b.y*(c.x-b.x) ) ) return false;
	if( 0 > sgn*( c.x*(a.y-c.y) - c.y*(a.x-c.x) ) ) return false;

    /*
    printf( "a  (%3.3f,%3.3f,%3.3f)\n", a_.x, a_.y, a_.z );
    printf( "b  (%3.3f,%3.3f,%3.3f)\n", b_.x, b_.y, b_.z );
    printf( "c  (%3.3f,%3.3f,%3.3f)\n", c_.x, c_.y, c_.z );
    printf( "hX (%3.3f,%3.3f,%3.3f)\n", hX.x, hX.y, hX.z );
    printf( "hY (%3.3f,%3.3f,%3.3f)\n", hY.x, hY.y, hY.z );
	Vec2d a; a.set( hX.dot(a_), hY.dot(a_) );
	Vec2d b; b.set( hX.dot(b_), hY.dot(b_) );
	Vec2d c; c.set( hX.dot(c_), hY.dot(c_) );
	double C1 = a.x*(b.y-a.y) - a.y*(b.x-a.x);
	double C2 = b.x*(c.y-b.y) - b.y*(c.x-b.x);
	double C3 = c.x*(a.y-c.y) - c.y*(a.x-c.x);
	printf( "C123 %f %f %f \n", C1, C2, C3 );
    if( (C1<0)||(C2<0)||(C3<0) ) return false;
	*/

    //printf("passed\n");
    return true;
}

inline double rayTriangle2(
	const Vec3d &ray0, const Vec3d &hRay, const Vec3d &hX, const Vec3d &hY,
	const Vec3d &a,    const Vec3d &b,    const Vec3d &c,
    Vec3d& normal
){
    //printf( "abc (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) \n", a.x,a.y,a.z,  b.x,b.y,b.z,  c.x,c.y,c.z );
    if( !rayInTriangle( a-ray0, b-ray0, c-ray0, hX, hY ) ) return t_inf;
    //printf("ray is in triangle\n");
	Vec3d ab,ac;
	ab.set_sub( b, a );
	ac.set_sub( c, a );
	normal.set_cross( ab, ac );
	//printf( "rayTriangle2 ab(%g,%g,%g) ac(%g,%g,%g) normal(%g,%g,%g)\n",  ab.x,ab.y,ab.z,  ac.x,ac.y,ac.z,  normal.x,normal.y,normal.z );
	double r2 = normal.norm2();
	if(r2<1.e-32){ return t_inf; }
	//printf( "rayTriangle2 ab (%g,%g,%g)\n", ab.x, ab.y, ab.z );
	//printf( "rayTriangle2 ac (%g,%g,%g)\n", ac.x, ac.y, ac.z );
    //printf( "rayTriangle2 normal (%g,%g,%g)\n", normal.x, normal.y, normal.z );
	return rayPlane( ray0, hRay, normal, a );
}

inline double rayTriangle2(
	const Vec3d &ray0, const Vec3d &hRay,
	const Vec3d &a,    const Vec3d &b,    const Vec3d &c,
    Vec3d& normal
){
    Vec3d hX,hY;
    hRay.getSomeOrtho(hX,hY);
	return rayTriangle2( ray0, hRay, hX, hY, a,b,c, normal );
}



// =========== Polygon

inline double rayPolygon(
	const Vec3d &ray0, const Vec3d &hRay, const Vec3d &hX, const Vec3d &hY,
	int n, const int * inds, const Vec3d * points,
    Vec3d& normal
){
    Vec3d tmp;
    Vec2d a,b;
	tmp.set_sub(points[inds[n-1]],ray0); a.set( hX.dot(tmp), hY.dot(tmp)  );
	tmp.set_sub(points[inds[0  ]],ray0); b.set( hX.dot(tmp), hY.dot(tmp)  );
	double sgn = a.x*(b.y-a.y) - a.y*(b.x-a.x);
	//printf( "0 %f  (%3.3f,%3.3f) (%3.3f,%3.3f)\n", sgn, a.x,a.y, b.x,b.y );
    for( int i=1; i<n; i++ ){
        a.set(b);
        tmp.set_sub(points[inds[i]],ray0); b.set( hX.dot(tmp), hY.dot(tmp)  );
        double sgn2 = ( a.x*(b.y-a.y) - a.y*(b.x-a.x) );
        //printf( "%i %i %f  (%3.3f,%3.3f) (%3.3f,%3.3f)\n", i, inds[i], sgn2, a.x,a.y, b.x,b.y );
        if( 0 > sgn*sgn2 ) return t_inf;
    }

    tmp.add(ray0);
	Vec3d ab,ac;
	ab.set_sub( points[inds[0]], tmp );
	ac.set_sub( points[inds[1]], tmp );
	normal.set_cross( ab, ac );

	return rayPlane( ray0, hRay, normal, tmp );
}



// =========== BoundingBox

//inline double rayCubeSide( double h, double ha, double hb, double amin, double amax, double bmin, double bmax){
//}

inline double rayBox( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& p1, const Vec3d& p2,  Vec3d& hitPos, Vec3d& normal ){
    Vec3d d1,d2;
    d1.set_sub( p1, ray0 );
    d2.set_sub( p2, ray0 );

    double tmin = t_inf ;
    //printf( " hRay (%3.3f,%3.3f,%3.3f) \n", hRay.x, hRay.y, hRay.z );

    // ---  along x
    if( (hRay.x*hRay.x) > 1e-32 ){

        double t,dir;
        if  ( hRay.x > 0 ){  t = d1.x/hRay.x;  dir = -1.0;  }
        else              {  t = d2.x/hRay.x;  dir =  1.0;  }
        //printf( " tx %f \n", t );
        //if( t < tmin ){
            double y = hRay.y * t;
            double z = hRay.z * t;
            //printf( " tx %f   y:  %3.3f %3.3f %3.3f    z: %3.3f %3.3f %3.3f \n", t, y, d1.y, d2.y,  z, d1.z, d2.z );
            if( ( y > d1.y ) && ( y < d2.y ) && ( z > d1.z ) && ( z < d2.z ) ){
                tmin = t;
                normal.set( dir,0,0 );
            }
        //}
    }else{
        if( (d1.x<0.0)||(d2.x>0.0) ){ // cannot hit
            return t_inf;
        }
    }

    // ---  along y
    if( (hRay.y*hRay.y) > 1e-32 ){
        double t,dir;
        if  ( hRay.y > 0 ){  t = d1.y/hRay.y;  dir = -1.0;  }
        else              {  t = d2.y/hRay.y;  dir =  1.0;  }
        //printf( " ty %f \n", t );
        if( t < tmin ){
            double x = hRay.x * t;
            double z = hRay.z * t;
            //printf( " tx %f   y:  %3.3f %3.3f %3.3f    z: %3.3f %3.3f %3.3f \n", t, x, d1.x, d2.x,  z, d1.z, d2.z );
            if( ( x > d1.x ) && ( x < d2.x ) && ( z > d1.z ) && ( z < d2.z ) ){
                tmin = t;
                normal.set( 0,dir,0 );
            }
        }
    }else{
        if( (d1.y<0.0)||(d2.y>0.0) ){ // cannot hit
            return t_inf ;
        }
    }

    // ---  along z
    if( (hRay.z*hRay.z) > 1e-32 ){
        double t,dir;
        if  ( hRay.z > 0 ){  t = d1.z/hRay.z;   dir = -1.0; }
        else              {  t = d2.z/hRay.z;   dir =  1.0; }
        //printf( " tz %f \n", t );
        if( t < tmin ){
            double x = hRay.x * t;
            double y = hRay.y * t;
            //printf( " tz %f   x:  %3.3f %3.3f %3.3f    y: %3.3f %3.3f %3.3f \n", t, x, d1.x, d2.x,  y, d1.y, d2.y );
            if( ( x > d1.x ) && ( x < d2.x ) && ( y > d1.y ) && ( y < d2.y ) ){
                tmin = t;
                normal.set( 0,0,dir );
            }
        }
    }else{
        if( (d1.z<0.0)||(d2.z>0.0) ){ // cannot hit
            return t_inf;
        }
    }

    hitPos.set(ray0);
    hitPos.add_mul( hRay, tmin );
    return tmin;
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
