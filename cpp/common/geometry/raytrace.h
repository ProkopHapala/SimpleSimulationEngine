
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

template<typename T>
inline T rayPointDistance2( const Vec3T<T>& ray0, const Vec3T<T>& hRay, const Vec3T<T>& point, T& t ){
	Vec3T<T> pt;
	pt.set_sub( point, ray0 );
	t  = hRay.dot( pt );;
    pt.add_mul( hRay, -t );
	return pt.norm2();
}

/*
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
*/

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

inline double capsulaIntersect( const Vec3d& ro, const Vec3d&  rd, const Vec3d&  pa, const Vec3d&  pb, double r ){
    // https://www.shadertoy.com/view/Xt3SzX
    Vec3d  ba = pb - pa;
    Vec3d  oa = ro - pa;
    double baba = ba.norm2();
    double oaoa = oa.norm2();
    double bard = ba.dot(rd);
    double baoa = ba.dot(oa);
    double rdoa = rd.dot(oa);
    double a = baba      - bard*bard;
    double b = baba*rdoa - baoa*bard;
    double c = baba*oaoa - baoa*baoa - r*r*baba;
    double h = b*b - a*c;
    if( h>=0 ){
        float t = (-b-sqrt(h))/a;
        float y = baoa + t*bard;
        // body
        if( y>0.0 && y<baba ) return t;
        // caps
        Vec3d oc = (y<=0.0) ? oa : ro - pb;
        b = rd.dot(oc);
        c = oc.dot(oc) - r*r;
        h = b*b - c;
        if( h>0.0 ) return -b - sqrt(h);
    }
    return 1e+300;
}

// =========== Plane

inline double rayPlane( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& normal, const Vec3d& point ){
	double nh = normal.dot( hRay  );
    if( nh*nh < 1.e-300 ){ return t_inf; }; // hRay and normal are perpendicular, they cannot intersect
	//double np = normal.dot( point );
	//double n0 = normal.dot( ray0  );
    //double t  = ( np - n0 ) / nh;
    // printf( "rayPlane t: %g | <n|h>: %g <n|p>: %g <n|o> %g\n", t, nh, np, n0 );
    Vec3d  d = point - ray0;
    double t = normal.dot( d ) / nh;
	//printf( "rayPlane t: %g | <n|h>: %g <n|d> %g \n", t, nh, d );
	//return ( n0 + np ) / nh;
	return t;
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

inline double rayTriangles( const Vec3d& ray0, const Vec3d& hRay, int n, const Vec3i* tris, const Vec3d* points, Vec3d& normal, int& imin ){
    Vec3d hX,hY;
    hRay.getSomeOrtho( hX, hY );
    double t_min = 1e+300;
    imin = -1;
    for(int i=0; i<n; i++ ){
        const Vec3i& it = tris[i];
        const Vec3d& A  = points[it.x];
        const Vec3d& B  = points[it.y];
        const Vec3d& C  = points[it.z];
        Vec3d normal_;
        bool inside_;
        double t = rayTriangle2( ray0, hRay, hX, hY, A, B, C, normal_ );
        inside_ = (t<0.9e+300 )&&(t>0);
        if( inside_ && ( t<t_min ) ){
            t_min  = t;
            normal = normal_;
            imin   = i;
        }
    }
    return t_min;
}

inline int pickParticle( const Vec3d& ray0, const Vec3d& hRay, double R, int n, Vec3d * ps, bool* ignore=0 ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<n; i++){
        if(ignore)if(ignore[i])continue;
        double ti = raySphere( ray0, hRay, R, ps[i] );
        //printf( "atom %i t %g %g hRay(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, tmin, hRay.x,hRay.y,hRay.z,  ps[i].x, ps[i].y, ps[i].z );
        //printf( "atom %i t %g %g %g ray0(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, R, tmin, ray0.x,ray0.y,ray0.z,  ps[i].x, ps[i].y, ps[i].z );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    return imin;
}

inline int pickPoinMinDist( const Vec3d &ray0, const Vec3d &hRay, int n, const Vec3d* points ){
    double r2min=1e+300;
    int imin=0;
    for(int i=0; i<n; i++){
        double t;
        double r2 = rayPointDistance2( ray0, hRay, points[i], t );
        if(r2<r2min){ imin=i; r2min=r2; }
    }
    return imin;
};

inline int pickBondCenter( int nb, const Vec2i* bonds, const Vec3d* ps, const Vec3d& ray0, const Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    //printf( "ray0(%g,%g,%g) hRay(%g,%g,%g) \n", ray0.x,ray0.y,ray0.z, hRay.x,hRay.y,hRay.z );
    for(int i=0; i<nb; i++){
        const Vec2i& b = bonds[i];
        Vec3d p = (ps[b.a]+ps[b.b])*0.5;
        double ti = raySphere( ray0, hRay, R, p );
        //printf( "ti[%i] %g  p(%g,%g,%g)\n", i, ti, p.x,p.y,p.z );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    //printf( ">>>> imin %i tmin %g \n", imin, tmin );
    return imin;
}

template<typename Func>
int rayPickBond( const Vec3d& ray0, const Vec3d& hRay, int nb, Func func, const double Rmax=0.5, const bool byCenter=false ){
    //printf( "rayPickBond() byCenter=%i Rmax=%g ray0(%g,%g,%g)\n", byCenter, Rmax, ray0.x,ray0.y,ray0.z );
    double R2max = Rmax*Rmax;
    double dist_min =  1e+300;
    int    imin = -1;
    for(int i=0; i<nb; i++){
        Vec3d pa,pb; func(i,pa,pb);
        if(byCenter){
            double ti = raySphere( ray0, hRay, Rmax, (pa+pb)*0.5 );
            if(ti<dist_min){ imin=i; dist_min=ti; }
        }else{
            double t1,t2;
            Vec3d hb = pb-pa; double l = hb.normalize();
            double dist = rayLine( ray0, hRay, pa, hb, t1, t2 );
            //printf( "rayPickBond()[%i] dist=%g rs(%g,%g)\n", i, dist, t1, t2 );
            if( (dist<Rmax) && (t2>0) && (t2<l) ){
                //printf( "rayPickBond() picked %i dist=%g dist_min=%g\n", i, dist, dist_min  );
                imin=i; dist_min=dist;
            }
        }
    }
    return imin;
}


/*
inline bool pointInRect( const Vec3d& p, const Mat3d& rot, Vec2d ){
    double x =
}
*/

// =========== BoundingBox

//inline double rayCubeSide( double h, double ha, double hb, double amin, double amax, double bmin, double bmax){
//}

inline double rayBox( const Vec3d& ray0, const Vec3d& hRay, const Vec3d& p1, const Vec3d& p2,  Vec3d& hitPos, Vec3d& normal ){
    // we assume that p1 and p2 are opposite corners of the box, and are already sorted (p1.x<p2.x, p1.y<p2.y, p1.z<p2.z)
    Vec3d d1,d2;
    d1.set_sub( p1, ray0 );
    d2.set_sub( p2, ray0 );

    double tmin = t_inf ;
    //printf( " hRay (%3.3f,%3.3f,%3.3f) \n", hRay.x, hRay.y, hRay.z );

    // ---  along x
    if( (hRay.x*hRay.x) > 1e-32 ){

        double t,dir;
        if  ( hRay.x > 0 ){  t = d1.x/hRay.x;  dir = -1.0;  } // we use assumption that p1.x<p2.x here
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
        if  ( hRay.y > 0 ){  t = d1.y/hRay.y;  dir = -1.0;  }  // we use assumption that p1.y<p2.y here
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
        if  ( hRay.z > 0 ){  t = d1.z/hRay.z;   dir = -1.0; }   // we use assumption that p1.y<p2.y here
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

/*
template<typename T>
inline int rayBox2( const Vec3TYPE<T>& o, const Vec3TYPE<T>& h, const Vec3TYPE<T>& p1, Vec3TYPE<T>& p2, T& tin, T& tout, Vec3TYPE<T>* normal ){
    Vec3TYPE<T> inv; inv.set_inv( hRay );
    Vec3d t1; t1.set_add_mul( o, p1, inv ); // hit to planes passing through p1
    Vec3d t2; t2.set_add_mul( o, p2, inv ); // hit to planes passing through p2

    if( t1.x > t2.x ){ T t=t1.x; t1.x=t2.x; t2.x=t; }
    if( t1.y > t2.y ){ T t=t1.y; t1.y=t2.y; t2.y=t; }
    if( t1.z > t2.z ){ T t=t1.z; t1.z=t2.z; t2.z=t; }

    // is closer hit to x-plane within y,z rectangle ?
    T y = o.y + h.y*t1.x;
    T z = o.z + h.z*t1.x;
    if( y )
    return tmin;
}
*/

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
