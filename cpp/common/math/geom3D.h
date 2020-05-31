
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
#include "raytrace.h"


/*

ToDo:

classes should be differentiated by underlaying data-type, and casted between them
LineSegment{ Vec3d a ,b;                        }
LineSegment{ Vec3d p0,hdir; double l;           }
Capsula    { Vec3d p0,hdir; double l; double r; }
LineSegment{ Vec3d a,b;  }
Ray        { Vec3d p0,hdir;  }
AABB{ Vec3d a,b; }

Sphere{ Vec3d pos; double r;         }
Plane { Vec3d dir; double c;         }
Slab  { Vec3d dir; double cmin,cmax; }


*/
















inline Vec3d cog_of_points ( int n, Vec3d * points ){ Vec3d c;  c.set(0.0); for(int i=0;i<n; i++){ c.add(points[i]); }  c.mul(1.0d/n); return c; }

inline double Rbound2( const Vec3d& center, int n, Vec3d * points ){
    double r2max=0;
    for(int i=0;i<n; i++){
        Vec3d p; p.set_sub( points[i], center );
        double r2 = p.norm2();
        if(r2>r2max) r2max = r2;
    }
    return r2max;
}

class Sphere3d{ public:
    Vec3d   p;
    double  r;

    inline double ray( const Vec3d& ray0, const Vec3d& hRay ){
        return raySphere( ray0, hRay, r, p );
    }

    //inline void normal( Vec3d& normal ){
    //    sphereNormal( double t, ray0, hRay, p, normal ){
    //};

};

class Line3d{ public:
    Vec3d  a,b;

    void fromSubLine( const Vec3d& a_, const Vec3d& b_, double ca, double cb ){
        Vec3d d; d.set_sub(b_,a_);
        a.set_add_mul( a_, d, ca);
        b.set_add_mul( a_, d, cb);
    }

};

struct Rayt3d{
    Vec3d  p0,hdir;
    double t;
};

class LineInterval3d{ public:
    Vec3d  p0,hdir;
    double t0,t1;

    inline bool fromPlanes( Vec3d& dir1, double c1, Vec3d& dir2, double c2 ){
        double s = dir1.dot(dir2);
        if( s*s>0.999 ) return false;
        double s2=s*s;
        double denom = 1/(1-s2);
        p0.set_lincomb(  (c1-c2*s)*denom, dir1, (c2-c1*s)*denom,  dir2 );
        hdir.set_cross( dir1, dir2 );
        hdir.mul( sqrt(denom) );
        //printf( " s %g denom %g hdir (%g,%g,%g) p0 (%g,%g,%g) \n", s, denom, hdir.x, hdir.y, hdir.z,   p0.x, p0.y, p0.z );
        infiniteSpan();
        return true;
    };

    inline int trim( Vec3d& dir, double c ){
        double s = dir.dot(hdir);
        double d = c - dir.dot(p0);
        if( s*s<0.0001 ){
            //printf( "perpendiculer \n" );
            if( d > 0 ){ double t=t0; t0=t1; t1=t; }; // always outside
            return 0;
        }
        double t = d/s;
        //printf( " s %g c %g cp0 %g t %g | p0 (%g,%g,%g) dir (%g,%g,%g)  %g \n", s, c, dir.dot(p0), t,    p0.x,p0.y,p0.z,  dir.x,dir.y,dir.z , dir.dot(p0)  );
        if( s>0 ){ if(t>t0){t0=t; return -1; }; }
        else     { if(t<t1){t1=t; return  1; }; }
        //printf( "no-trim \n" );
        return 0;
        //return t1<t0;
    };

    inline void  infiniteSpan(){ t0 = -1e+300; t1 =  1e+300; }
    inline Vec3d endPoint0(){ return p0+hdir*t0; };
    inline Vec3d endPoint1(){ return p0+hdir*t1; };

};




// ### Quadric Intersection
// http://skuld.bmsc.washington.edu/people/merritt/graphics/quadrics.html
// F(x, y, z) = Ax2 + By2 + Cz2 + Dxy+ Exz + Fyz + Gx + Hy + Iz + J = 0
// Then substitute in ray equation R(t) = Ro + Rd and we get quadratic equation:
//   Aqt2 + Bqt + Cq = 0
//   Aq = Axd2 + Byd2 + Czd2 + Dxdyd + Exdzd + Fydzd
//   Bq = 2*Axoxd + 2*Byoyd + 2*Czozd + D(xoyd + yoxd) + E(xozd + zoxd) + F(yozd + ydzo) + Gxd + Hyd + Izd
//   Cq = Axo2 + Byo2 + Czo2 + Dxoyo + Exozo + Fyozo + Gxo + Hyo + Izo + J

inline bool quaric_intersec( Vec3d qp, Vec3d qb, Vec3d q, double R, Vec3d rd, Vec3d ro, double& t1, double& t2 ){
    double A = qp.bidot(rd,rd)   + qb.antidot(rd,rd);
    double B = qp.bidot(ro,rd)*2 + qb.antidot(ro,rd) + qb.antidot(rd,ro)  + q.dot(rd);
    double C = qp.bidot(ro,ro)   + qb.antidot(ro,ro)                      + q.dot(ro) + R;
    return quadratic_roots ( A, B, C, t1, t2 );
}

inline void rayQuadric( Mat3d rot, Vec3d pos, Vec3d sc, double R, Vec3d rd, Vec3d ro, double& t1, double& t2 ){
// Ax^2 + Bx^2 + Cx^2 = 0   // no parabola
// Ax^2 + Bx^2 + Cx^2 = 0
// ToDo : How to conveniently define General Unitary quadric (un-rotated, un-scaled)
};


class Capsula3D{ public:
    Vec3d  p,hdir;
    double r,l;

    inline double dist2( const Vec3d& pos ) const {
        Vec3d d; d.set_sub( pos, p );
        double x  = d.makeOrthoU(hdir);
        double y2 = d.norm2();
        double dist2 = y2;
        if (x <0    ){ dist2 += sq(x); }else if(x>l){ dist2 += sq(x-l); };
        return dist2;
    }

    inline double dist2_Cilinder( const Vec3d& pos ) const {
        Vec3d d; d.set_sub( pos, p );
        double x  = d.makeOrthoU(hdir);
        double y2 = d.norm2();
        double dist2 = 0.0;
        if (x <0    ){ dist2 += sq(x); }else if(x>l){ dist2 += sq(x-l); };
        if (y2>sq(r)){ dist2 += sq(sqrt(y2)-r); };
        //printf( "x %g y %g dist2 %g  |  l %g r %g \n", x, sqrt(y2), dist2,   l, r );
        return dist2;
    }

    void enclosePoints( int n, Vec3d* ps ){
        double lmin = 0.0;
        double lmax = l;
        double r2max=0.0;
        for( int i=0; i<n; i++ ){
            Vec3d d; d.set_sub( ps[i], p );
            double r2i,li;
            li =d.makeOrthoU(hdir);
            r2i=d.norm2();
            if(r2i>r2max)r2max=r2i;
            double dr2 = r2max-r2i;
            if      (li>lmax){ if( sq(li-lmax)>dr2 ){ lmax=li-sqrt(dr2);  } }
            else if (li<lmin){ if( sq(li-lmin)>dr2 ){ lmin=sqrt(dr2)-li;  } }
        }
        r = sqrt(r2max);
        l = lmax-lmin;
        p.set_add_mul( p, hdir, lmin );
	}

    void enclosePoints_Cylinder( int n, Vec3d* ps ){
        double lmin =0.0;
        double lmax =l;
        double r2max=0.0;
        for( int i=0; i<n; i++ ){
            Vec3d d; d.set_sub( ps[i], p );
            double r2i,li;
            li =d.makeOrthoU(hdir);
            r2i=d.norm2();
            if(r2i>r2max)r2max=r2i;
            if      (li>lmax){lmax=li;}
            else if (li<lmin){lmin=li;}
        }
        r = sqrt(r2max);
        l = lmax-lmin;
        p.set_add_mul( p, hdir, lmin );
	}

};

class Box{ public:
    Vec3d a,b;

    // ==== functions

    inline void setSymetric( const Vec3d& v ){ b={fabs(v.x),fabs(v.y),fabs(v.z)};  a.set_mul(b,-1); };
    inline void setSymetric( double f ){ f=fabs(f); a.set(-f,-f,-f); b.set(f,f,f); };
    inline void shift      ( const Vec3d& v ){ a.add(v); b.add(v);   };
    inline void scale      ( const Vec3d& v ){ a.mul(v); b.mul(v);   };
    inline void scale      ( double f ){ a.mul(f); b.mul(f);   };

    inline void setOrdered(const Vec3d& a_,const Vec3d& b_){
        if(a_.x<b_.x){ a.x=a_.x; b.x=b_.x; }else{  a.x=b_.x; b.x=a_.x; };
        if(a_.y<b_.y){ a.y=a_.y; b.y=b_.y; }else{  a.y=b_.y; b.y=a_.y; };
        if(a_.z<b_.z){ a.z=a_.z; b.z=b_.z; }else{  a.z=b_.z; b.z=a_.z; };
    };

    inline void order(){
        if(a.x>b.x){ double f=a.x; a.x=b.x; b.x=f; };
        if(a.y>b.y){ double f=a.y; a.y=b.y; b.y=f; };
        if(a.z>b.z){ double f=a.z; a.z=b.z; b.z=f; };
    };

    inline static bool pointIn( const Vec3d& p, const Vec3d& a, const Vec3d& b ){
        return ((p.x>a.x)&&(p.y>a.y)&&(p.z>a.z)&&
                (p.x<b.x)&&(p.y<b.y)&&(p.z<b.z));
    }

    inline static double dist2( const Vec3d& p, const Vec3d& a, const Vec3d& b ){
        // from here : http://stackoverflow.com/questions/4578967/cube-sphere-intersection-test
        // assume C1 and C2 are element-wise sorted, if not, do that now
        double dist2 = 0.0d;
        if (p.x < a.x){ dist2 += sq(p.x - a.x); }else if(p.x > b.x){ dist2 += sq(p.x - b.x); };
        if (p.y < a.y){ dist2 += sq(p.y - a.y); }else if(p.y > b.y){ dist2 += sq(p.y - b.y); };
        if (p.z < a.z){ dist2 += sq(p.z - a.z); }else if(p.z > b.z){ dist2 += sq(p.z - b.z); };
        return dist2;
    }

    void fromPoints( int n, Vec3d * points ){
        a.set(points[0]);
        b.set(points[0]);
        for(int i=1; i<n; i++){
            Vec3d& p = points[i];
            if(p.x<a.x){a.x=p.x;}else if(p.x>b.x){b.x=p.x;}
            if(p.y<a.y){a.y=p.y;}else if(p.y>b.y){b.y=p.y;}
            if(p.z<a.z){a.z=p.z;}else if(p.z>b.z){b.z=p.z;}
        }
    }

    inline bool pointIn( const Vec3d& p ) const { return pointIn( p,a, b ); }
    inline bool dist2  ( const Vec3d& p ) const { return pointIn( p,a, b ); }

    inline bool pointRot( const Vec3d& p ) const {
        return ((p.x>a.x)&&(p.y>a.y)&&(p.z>a.z)&&
                (p.x<b.x)&&(p.y<b.y)&&(p.z<b.z));
    }

    inline bool overlap( const Box& box ) const {
        if ( (box.a.x>b.x) || (box.b.x<a.x) ) return false;
        if ( (box.a.y>b.y) || (box.b.y<a.y) ) return false;
        if ( (box.a.z>b.z) || (box.b.z<a.z) ) return false;
        return true;
    }

    inline void combine( const Box& A, const Box& B ){
        a.x = (A.a.x<B.a.x) ? A.a.x : B.a.x;
        a.y = (A.a.y<B.a.y) ? A.a.y : B.a.y;
        a.z = (A.a.z<B.a.z) ? A.a.z : B.a.z;
        b.x = (A.b.x>B.b.x) ? A.b.x : B.b.x;
        b.y = (A.b.y>B.b.y) ? A.b.y : B.b.y;
        b.z = (A.b.z>B.b.z) ? A.b.z : B.b.z;
    }

    inline void enclose( const Box& B ){
        // is this "if" faster than branch-less SIMD operation with ternary operator ?
        if(B.a.x<a.x){ a.x = B.a.x; };
        if(B.a.y<a.y){ a.y = B.a.y; };
        if(B.a.z<a.z){ a.z = B.a.z; };
        if(B.b.x>b.x){ b.x = B.b.x; };
        if(B.b.y>b.y){ b.y = B.b.y; };
        if(B.b.z>b.z){ b.z = B.b.z; };
    }

    inline Vec3d  center    ()const{ Vec3d c; c.set_add(b,a); c.mul(0.5); return c; }
    inline Vec3d  dimensions()const{ Vec3d d; d.set_sub(b,a); return d;           }
    inline double volume    ()const{ Vec3d d; d.set_sub(b,a); return d.x*d.y*d.z; }
    inline double surfArea()const{
        Vec3d d; d.set_sub(b,a);
        return 2.0*(  d.x*d.y + d.x*d.z + d.y*d.x );
    }

    inline void spanAlongDir(const Vec3d hdir, Vec2d& span) const {
        span = (Vec2d){+1e+300,-1e+300};
        double xdir;
        xdir = hdir.dot( {a.x,a.y,a.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {a.x,a.y,b.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {a.x,b.y,a.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {a.x,b.y,b.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {b.x,a.y,a.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {b.x,a.y,b.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {b.x,b.y,a.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
        xdir = hdir.dot( {b.x,b.y,b.z} ); _setmin(span.x,xdir); _setmax(span.y,xdir);
    }

    inline Vec3d genRandomSample() const { return (Vec3d){randf(a.x,b.x),randf(a.y,b.y),randf(a.z,b.z)}; }

    inline double ray   ( const Vec3d& ray0, const Vec3d& hRay, Vec3d& hitPos, Vec3d& normal ) const { return rayBox( ray0, hRay, a, b, hitPos, normal ); }
    inline double rayRot( Vec3d ray0, Vec3d hRay, const Mat3d& rot, Vec3d& hitPos, Vec3d& normal )const {
        rot.dot_to(ray0,ray0);
        rot.dot_to(hRay,hRay);
        return rayBox( ray0, hRay, a, b, hitPos, normal );
        rot.dot_to(hitPos,hitPos);
        rot.dot_to(normal,normal);
    }

};


class Tetrahedron{ public:
    Vec3d a,b,c,d;

};



inline double normalAreaTriangle( Vec3d& nr, const Vec3d& a, const Vec3d& b, const Vec3d& c ){
    nr.set_cross(b-a,c-a);
    return nr.normalize() * 0.5;
}

class Triangle3D{
    public:
	union{
		struct{ Vec3d a,b,c; };
		Vec3d array[3];
	};
	inline double normalArea(Vec3d& nr)const{
        //nr.set_cross(b-a,c-a);
        //return nr.normalize() * 0.5;
        return normalAreaTriangle( nr, a,b,c );
	}
	inline bool rayIn( const Vec3d& ray0, const Vec3d& hX, const Vec3d& hY )const{
        return rayInTriangle( a-ray0, b-ray0, c-ray0, hX, hY );
	}
    inline double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal, bool& inside, const Vec3d& hX, const Vec3d& hY )const{
        return rayTriangle2( ray0, hRay, hX, hY, a,b,c, normal );
	}
	inline double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal, bool& inside )const{
        Vec3d hX,hY;
        hRay.getSomeOrtho(hX,hY);
        return rayTriangle2( ray0, hRay, hX, hY, a,b,c, normal );
    }
};


class Quad3d{ public:
    // two triangles abc,  acd (cad)
	union{
		struct{ Vec3d a,b,c,d; };
		struct{ Vec3d p00,p01,p10,p11; };
		struct{ Line3d l1,l2; };
		Vec3d array[4];
	};
	inline double normalArea(Vec3d& nr)const{
        return  normalAreaTriangle( nr, a,b,c );
              + normalAreaTriangle( nr, a,c,d );
	};
	inline bool rayIn( const Vec3d& ray0, const Vec3d& hX, const Vec3d& hY )const{
        Vec3d a_=a-ray0;
        Vec3d b_=b-ray0;
        Vec3d c_=c-ray0;
        Vec3d d_=d-ray0;
        return rayInTriangle( a_, b_, d_, hX, hY ) || rayInTriangle( a_, d_, c_, hX, hY );
	}
    inline double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal, const Vec3d& hX, const Vec3d& hY )const{
        //printf( "====== abcd (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) \n", a.x,a.y,a.z,  b.x,b.y,b.z,  c.x,c.y,c.z,  d.x,d.y,d.z );
        double t =         rayTriangle2( ray0, hRay, hX, hY, a,b,d, normal );
        //printf( "triangle 1 : %g \n", t );
        if(t>0.9*t_inf){ t=rayTriangle2( ray0, hRay, hX, hY, a,d,c, normal );}
        //printf( "triangle 2 : %g \n", t );
        return t;
	}
	inline double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal )const{
        Vec3d hX,hY;
        hRay.getSomeOrtho(hX,hY);
        return ray( ray0, hRay, normal, hX, hY );
    }
};


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

	void fromVec(const Vec3d& v){
        normal=v;
        iso = normal.normalize();
	}

    void fromPointAndNormal( const Vec3d& normal_, const Vec3d& p0 ){
        normal = normal_;
        C = normal.dot(p0);
	};


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

	double dist( Vec3d point )const{ return normal.dot( point ) - iso; };

};


inline bool dist_Box( const Vec3d& point,    const Vec3d& pos, const Mat3d& orientation, const Vec3d& span ){
    // distance from bounding ellipsoide
    Vec3d d;
    d.set_sub( point, pos );
    d = orientation.dot( d );
    double dm;
    dm =     ( fabs(d.a) - span.a     );
    dm = fmax( fabs(d.b) - span.b, dm );
    dm = fmax( fabs(d.c) - span.c, dm );
    return dm;
};

/*
inline double ray_Box( const Vec3d& ray0, const Vec3d& hRay,    const Vec3d& pos, const Mat3d& orientation, const Vec3d& span ){
    // WARRNING: unfinished
    Vec3d d ;
    d.set_sub(ray0, pos);
    d = orientation.dot( d );
    return ;
}
*/


// ============ Ellipsoide

inline double dist_Ellipsoide( const Vec3d& point, const Vec3d& pos, const Mat3d& orientation, const Vec3d& span ){
    // distance from bounding ellipsoide
    Vec3d d;
    d.set_sub( point, pos );
    d = orientation.dot( d );
    d.mul(span);
    return 1 - d.norm2();
};

inline Vec3d normal_Ellipsoide( const Vec3d& phit, const Vec3d& pos, const Mat3d& orientation, const Vec3d& invspan ){
    Vec3d normal_ = phit;
    //normal_.set_add_mul(ray0_, hRay_, thit );
    normal_.mul(invspan); normal_.mul(invspan);
    return orientation.dotT( normal_ );
}

inline double ray_Ellipsoide( const Vec3d& ray0, const Vec3d& hRay, Vec3d * normal, const Vec3d& pos, const Mat3d& orientation, const Vec3d& span ){
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
        normal_.mul(invspan); normal_.mul(invspan);
        *normal = orientation.dotT( normal_ );
    }
    return thit/lRay;
};

/*
inline double ray_Ellipsoide( const Vec3d& ray0, const Vec3d& hRay,   const Vec3d& pos, const Mat3d& orientation, const Vec3d& span ){
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
        return thit/lRay;
}
*/




class Ellipsoide{
    public:
	Vec3d pos;
	Vec3d span;
    Mat3d orientation;

    inline void initOne(){
        pos.set(0.0,0.0,0.0);
        span.set(1.0,1.0,1.0);
        orientation.a.set(1.0,0.0,0.0);
        orientation.b.set(0.0,1.0,0.0);
        orientation.c.set(0.0,0.0,1.0);
    }

    inline void fromVecs( const Mat3d& mat ){
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

    inline double ray( const Vec3d& ray0, const Vec3d& hRay ){
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
        return thit/lRay;
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
            normal_.mul(invspan); normal_.mul(invspan);
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

    int findEdgeIndex( int ip1, int ip2 ){
        int n = ipoints.size();
        int ip = ipoints[0];
        if      (ip==ip1){
            if(ipoints[1]==ip2){return 0;}else{ return n-1; }
        }else if(ip==ip2){
            if(ipoints[1]==ip1){return 0;}else{ return n-1; }
        }
        int i;
        for(i=1; i<n-1; i++){
            if( ipoints[i]==ip1 ) break;
        }
        if(ipoints[i-1]==ip2){return (i-1); }else{ return i; }
    }

    void insertPoint( int ip, int ito ){
        auto it = ipoints.begin();
        ipoints.insert(it+ito, ip );
    }

    void removePoint( int ip ){
        //auto it = find( ipoints.begin(), ipoints.end(), ip ); ipoints.erase( it );
        int i;
        for(i=0; i<ipoints.size(); i++ ){ if(ipoints[i]==ip){ break;} }
        if( i<ipoints.size() ) ipoints.erase( ipoints.begin()+i );

    }

    void printPoints( ){ for( int i=0; i<ipoints.size(); i++ ){ printf("%i ", ipoints[i]); } printf("\n"); }

};

class MeshEdge{ public:
    Vec2i verts;
    Vec2i faces;

    void setVerts( int ia, int ib ){ if(ia<ib){verts.set(ia,ib);}else{verts.set(ib,ia);} };
};

#endif


