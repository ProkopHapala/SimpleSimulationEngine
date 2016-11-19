
// read also:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

#ifndef  geom3D_h
#define  geom3D_h

#include <vector>
#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"



inline double dist2_PointBox( Vec3d C1, Vec3d C2, Vec3d S){
    // from here : http://stackoverflow.com/questions/4578967/cube-sphere-intersection-test
    // assume C1 and C2 are element-wise sorted, if not, do that now
    double dist2 = 0.0d;
    if (S.x < C1.x){ dist2 += sq(S.x - C1.x); }else if(S.x > C2.x){ dist2 += sq(S.x - C2.x); };
    if (S.y < C1.y){ dist2 += sq(S.y - C1.y); }else if(S.y > C2.y){ dist2 += sq(S.y - C2.y); };
    if (S.z < C1.z){ dist2 += sq(S.z - C1.z); }else if(S.z > C2.z){ dist2 += sq(S.z - C2.z); };
    return dist2;
}


// ============ Plane3D

class Plane3D{
    public:
	union{
		struct{ double x,y,z,C;	        };
		struct{ Vec3d normal; double iso; };
		double array[4];
	};

	void set( const Vec3d& normal_, double iso_ ){ normal.set(normal_); iso=iso_; };
	void set( double x_, double y_, double z_, double C_){ x=x_; y=y_; z=z_; C=C_; };

	void fromPoints_noNorm( const Vec3d& a, const Vec3d& b, const Vec3d& c ){
        Vec3d ab,bc;
		ab.set_sub(b,a);
		bc.set_sub(c,b);
		normal.set_cross(ab,bc); // we don't need to renormalize
		normal.dot( a );
	};

	void fromPoints( const Vec3d& a, const Vec3d& b, const Vec3d& c ){
		fromPoints_noNorm( a, b, c );
		iso /= normal.normalize();
	};

	double dist( Vec3d point ){ return normal.dot( point ) - iso; };

};

// ============ Ellipsoide

class Ellipsoide{
    public:
	Vec3d pos;
	Vec3d span;
    Mat3d orientation;

    inline bool initOne(){
        pos.set(0.0,0.0,0.0);
        span.set(1.0,1.0,1.0);
        orientation.a.set(1.0,0.0,0.0);
        orientation.b.set(0.0,1.0,0.0);
        orientation.c.set(0.0,0.0,1.0);
    }

    inline bool fromVecs( const Mat3d& mat ){
        orientation.set(mat);
        span.a = orientation.a.normalize();
        span.b = orientation.b.normalize();
        span.c = orientation.c.normalize();
    }

    inline bool pointIn( const Vec3d& point ){
        // distance from bounding ellipsoide
        Vec3d d;
        d.set_sub( point, pos );
        d = orientation.dot( d );
        d.mul(span);
        return d.norm2()<1;
    };

    inline double ray( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal ){
        // distance from bounding volume
        // transform space
        //printf("here\n");
        Vec3d hRay_,ray0_,p,invspan;
        ray0_.set_sub(ray0, pos);
        invspan.set_inv(span);
        Mat3d m; m.set(orientation);
        m.mul(invspan);
        //m.mulT(p);
        //ray0_  = m.dotT(ray0_);
        //hRay_  = m.dotT(hRay);
        ray0_  = m.dot(ray0_);
        hRay_  = m.dot(hRay);
        //ray0_  = orientation.dotT(ray0_); ray0_.mul( p );
        //hRay_  = orientation.dotT(hRay);  hRay_.mul( p );
        //printf( " hRay_ (%3.3f,%3.3f,%3.3f) \n",  hRay_.x, hRay_.y, hRay_.z    );
        double lRay = hRay_.normalize(); // this sqrt() may be possible to optimize out
        // ray sphere
        double tdisk = -hRay_.dot( ray0_ );
        //printf( " tdisk %f ray0_ (%3.3f,%3.3f,%3.3f) \n", tdisk,  ray0_.x, ray0_.y, ray0_.z );
        p.set_add_mul( ray0_, hRay_, tdisk );
        double r2   = p.norm2();
        //printf( " r2 %f p (%3.3f,%3.3f,%3.3f) \n", r2,  p.x, p.y, p.z    );
        if( r2 > 1 ) return 1e+300;
        double thit = tdisk - sqrt( 1 - r2 );
        if(normal){
            // ======== NOT CORRECT; TO DO LATER
            // E = (x/Lx)^2 + (y/Ly)^2  + (z/Lz)^2
            // (dE/dx,dE/dy,dE/dz) = ( 2*x/Lx^2, 2*y/Ly^2, 2*x/Lz^2 )
            // normal
            Vec3d normal_;
            normal_.set_add_mul(ray0_, hRay_, thit );
            normal_.mul(invspan);
            normal_.mul(invspan);
            *normal = orientation.dotT( normal_ );
        }
        return thit/lRay;
    };

};

class Disk3D{  // used for acceleration of raytracing
    public:
    Vec3d    center;
    double   Rbound;
    Plane3D  plane;
};

class Polygon{
    public:
    //Disk3D * disk;
    std::vector<int> ipoints;
    std::vector<int> inormals;
    std::vector<int> iedges;
    //bool convex=false;
};







#endif


