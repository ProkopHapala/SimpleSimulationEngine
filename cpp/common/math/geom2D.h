
// read also:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

#ifndef  geom2D_h
#define  geom2D_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"


// TODO : Static vs dynamic data ... e.g. should be Line2d{ Vec2 *A,*B; } or Line2d{ Vec2 A,B; } ?

//         Basic fast functions


////////////////////////////////
//   FUNCTIONS : pointIn shapes
////////////////////////////////



template <class T>
inline double triangleArea( const Vec2T<T>& a, const Vec2T<T>& b, const Vec2T<T>& c ){
    //https://www.mathopenref.com/coordtrianglearea.html
    return ( a.x*(b.y-c.y) + b.x*(c.y-a.y) + c.x*(a.y-b.y) )*0.5;
}

template <class T>
inline void dtriangleArea( const Vec2T<T>& b, const Vec2T<T>& c, Vec2T<T>& da ){
    da.x =  (b.y-c.y)*0.5;
    da.y = -(b.x-c.x)*0.5;
}


template <class T>
inline bool pointInRect( const Vec2T<T>& p, T x0, T y0, T x1, T y1 ){ return ( p.x > x0 ) && ( p.x < x1 ) && ( p.y > y0 ) && ( p.y < y1 ); };

template <class T>
inline bool pointInRect( const Vec2T<T>& p, const Vec2T<T>& vmax, const Vec2T<T>& vmin ){ return ( p.x > vmin.x ) && ( p.x < vmax.x ) && ( p.y > vmin.y ) && ( p.y < vmax.y ); };

template <class T>
inline bool pointInCircle( const Vec2T<T>& p, const Vec2T<T>& center, T rmax ){
	Vec2T<T> d;
	d.set_sub( p, center );
	T r2 = d.norm2();
	return (r2 < ( rmax * rmax ) );
}

template <class T>
inline bool pointInConvexPolygon( const Vec2T<T>& p, T * buf, int n ){
    //printf( " point %f %f \n", p.x, p.y );
	for ( int i=0; i<n; i++ ){
		T cv = p.dot( *(Vec2d*)buf );
        //printf( " dir %f %f   cv %f %f   \n", (*(Vec2d*)buf).x, (*(Vec2d*)buf).y,    cv, buf[2] );
		if( cv < buf[2] ) return false;
		buf += 3;
	}
	//printf( " returning true \n" );
	return true;
}

template <class T>
inline bool pointInConvexParalelogram( const Vec2T<T>& p, T * buf, int n ){
	for ( int i=0; i<n; i++ ){
		T cv = p.dot( *buf );
		if( ( cv < buf[2] ) || ( cv > buf[3] ) ) return false;
		buf += 4;
	}
	return true;
}


/////////////////////////
//   CLASS :   Rect2d
/////////////////////////

class Rect2d{
	public:
	union{
		struct{ Vec2d a,b; };
		struct{ double x0,y0,x1,y1; };
		double array[4];
	};

	inline void set ( double x0_, double y0_, double x1_, double y1_ ){ x0=x0_; y0=y0_; x1=x1_; y1=y1_; };
	inline void set ( const Vec2d& a, const Vec2d b ) {
		x0 = a.x; x1 = b.x; y0 = a.y; y1 = b.y;
		if( x0 > x1 ){ SWAP( x0, x1, double ) }
		if( y0 > y1 ){ SWAP( y0, y1, double ) }
	}

	//inline      Rect2d   (){};
	//inline      Rect2d   ( double x0_, double y0_, double x1_, double y1_ ){ x0=x0_; y0=y0_; x1=x1_; y1=y1_; };
	//inline      Rect2d   ( const Vec2d& a, const Vec2d b                  ){ x0=x0_; y0=y0_; x1=x1_; y1=y1_; };

	inline bool pointIn( const Vec2d& p ) const { return pointInRect( p, x0, y0, x1, y1 );    		          }

	inline bool notOverlaps( const Rect2d& r ) const {
        return ( x0 > r.x1 ) || ( x1 < r.x0 ) || ( r.y0 > r.y1 ) || ( y1 < r.y0 );
	}

	inline void setEmpty(){ x0=1e+300; y0=1e+300; x1=-1e+300; y1=-1e+300; }
    inline uint8_t enclose( const Vec2d& p ){
        uint8_t mask=0;
        if(p.x<x0){ x0=p.x; mask |=1; }
        if(p.x>x1){ x1=p.x; mask |=4; }
        if(p.y<y0){ y0=p.y; mask |=2; }
        if(p.y>y1){ y1=p.y; mask |=8; }
		return mask;
	}
    inline void   margin(double R){ x0-=R; y0-=R; x1+=R; y1+=R; }
	inline Vec2d  cog() const     { return (Vec2d){(x0+x1)*0.5,(y0+y1)*0.5}; }
	inline double l2Diag() const  { return sq(x1-x0) + sq(y1-y0); }

};


/////////////////////////
//   CLASS :   Line2d
//////////////////////////

class Line2d{
	public:
	union{
		struct{ double a,b,c;        };
		struct{ Vec2d dir; double l; };
		double array[3];
	};

	inline void   set          ( double a_, double b_, double c_ ){ a=a_; b=b_; c=c_;                                   };
	inline void   set          ( const Vec2d& A, const Vec2d& B  ){ b = A.x - B.x; a = B.y - A.y; c = a*A.x + b*A.y;    };
	inline void   normalize    (                                 ){ double ir = sqrt( a*a + b*b ); a*=ir; b*=ir; c*=ir; };
	inline double dist_unitary ( const Vec2d& p                  ) const { return a*p.x + b*p.y - c;                    }; // Should normalize before

	inline void intersectionPoint( const Line2d& l, Vec2d& p ) const {
		double idet = 1/( a * l.b - b * l.a );
		p.x = ( c * l.b - b * l.c ) * idet;
		p.y = ( a * l.c - c * l.a ) * idet;
	}

    inline double intersection_t( const Vec2d& A, const Vec2d& dAB ) const{
        double idet = 1/( a * dAB.x + b * dAB.y );
        return ( c - a * A.x - b * A.y ) * idet;
	}

    inline unsigned char intersectionPoint( const Vec2d& A, const Vec2d& B, Vec2d& p ) const {
        Vec2d dAB;
        dAB.set_sub( B, A );
        double t  = intersection_t( A, dAB );
        unsigned char mask = 0;
        if ( t<0 ) { mask = mask | 1; };
        if ( t>1 ) { mask = mask | 2; };
        p.x = A.x + t * dAB.x;
        p.y = A.y + t * dAB.y;
        return mask;
	}

};

inline double line_side   ( const Vec2d& p, const Vec2d& a, const Vec2d& b ){ Line2d l; l.set( a, b ); return l.dist_unitary( p );  }

inline Vec2d dpLineSegment( const Vec2d& pos, const Vec2d& p1, const Vec2d& p2 ){
    Vec2d d  = p2-p1;
    Vec2d dp = pos-p1;
    double c = d.dot( dp );
    if( c>0 ){
        double r2 = d.norm2();
        if(c>r2){ dp=pos-p2;            }
        else    { dp.add_mul(d,-c/r2 ); }
    }
    return dp;
}

/*
// not usefull ... use Line2d::intersection_t( const Vec2d& A, const Vec2d& dAB ) instead
inline double intersection_t (  const Vec2d& l1a, const Vec2d& l1b, const Vec2d& l2a, const Vec2d& l2b ) {
    double dx1 = l1b.x - l1a.x;     double dy1 = l1b.y - l1a.y;
    double dx2 = l2b.x - l2a.x;     double dy2 = l2b.y - l2a.y;
	double x12 = ( l1a.x - l2a.x );
	double y12 = ( l1a.y - l2a.y );
	double invD = 1/( -dx2 * dy1 + dx1 * dy2 );
    return ( dx2 * y12 - dy2 * x12 ) * invD;
}
*/

inline void intersection_st (  const Vec2d& l1a, const Vec2d& l1b, const Vec2d& l2a, const Vec2d& l2b, double& s, double& t ) {
    double dx1 = l1b.x - l1a.x;     double dy1 = l1b.y - l1a.y;
    double dx2 = l2b.x - l2a.x;     double dy2 = l2b.y - l2a.y;
	double x12 = ( l1a.x - l2a.x );
	double y12 = ( l1a.y - l2a.y );
	double invD = 1/( -dx2 * dy1 + dx1 * dy2 );
    s = (-dy1 * x12 + dx1 * y12 ) * invD;
    t = ( dx2 * y12 - dy2 * x12 ) * invD;
}

inline unsigned char intersection_point (  const Vec2d& l1a, const Vec2d& l1b, const Vec2d& l2a, const Vec2d& l2b, Vec2d& p ) {
	double s,t;
	//intersection_st( l1a, l1b, l2a, l2b, s, t );
    double dx1 = l1b.x - l1a.x;     double dy1 = l1b.y - l1a.y;
    double dx2 = l2b.x - l2a.x;     double dy2 = l2b.y - l2a.y;
	double x12 = ( l1a.x - l2a.x );
	double y12 = ( l1a.y - l2a.y );
	double invD = 1/( -dx2 * dy1 + dx1 * dy2 );
    s = (-dy1 * x12 + dx1 * y12 ) * invD;
    t = ( dx2 * y12 - dy2 * x12 ) * invD;
	unsigned char mask = 0;
	if ( t<0 ) { mask = mask | 1; };
	if ( t>1 ) { mask = mask | 2; };
	if ( s<0 ) { mask = mask | 4; };
	if ( s>1 ) { mask = mask | 8; };
	//printf( " s t mask %f,%f %i \n", s, t, mask );
    p.x = l1a.x + (t * dx1 );
    p.y = l1a.y + (t * dy1 );
	return mask;
}


/////////////////////////
//   CLASS :   Segment2d
//////////////////////////

class Segment2d{
	public:
	int id;
	Vec2d *a,*b;
	inline      Segment2d( int id_ ){ id=id_;  };
	inline      Segment2d( int id_, Vec2d* a_, Vec2d* b_ ){ id=id_; a=a_; b=b_; };
	inline void setPoints( Vec2d* a_, Vec2d* b_             ){ a=a_; b=b_; };
	//inline void setPoints( const Vec2d& a_, const Vec2d& b_ ){ a=&a_; b=&b_; };

	inline bool inRect( const Vec2d& p ){
		double x0 = a->x; double x1 = b->x; double y0 = a->y; double y1 = b->y;
		if( x0 > x1 ){ SWAP( x0, x1, double ) }
		if( y0 > y1 ){ SWAP( y0, y1, double ) }
		//printf( " x0 x1 p.x %f %f %f       y0 y1 p.y %f %f %f \n", x0, x1, p.x, y0, y1, p.y );
		return pointInRect( p, x0, y0, x1, y1 );
	}

};

struct Ellipse2D{
    Vec2d pos;
    Vec2d dir;
    Vec2d sz;

    inline bool pointIn( Vec2d p ){
        p.sub(pos);
        double a = dir.dot     (p);
        double b = dir.dot_perp(p);
        return ( sq(a/sz.a) + sq(b/sz.b) ) > 1.0;
    }

    inline bool pointDistInterp( Vec2d p ){
        p.sub(pos);
        double a = dir.dot     (p);
        double b = dir.dot_perp(p);
        return 1.0 - sq(a/sz.a) + sq(b/sz.b);
    }

};


/////////////////////////
//   CLASS :   Triangle2d
//////////////////////////

class Triangle2d{ public:
	int id;
	Vec2d *a,*b,*c;
	inline      Triangle2d( int id_ ){ id=id_;  };
	inline      Triangle2d( int id_, Vec2d* a_, Vec2d* b_, Vec2d* c_ ){ id=id_; a=a_; b=b_; c=c_; };
	inline void setPoints ( Vec2d* a_, Vec2d* b_, Vec2d* c_ ){ a=a_; b=b_; c=c_; };

	inline bool pointIn( const Vec2d& p ){
		double inward = line_side( *c, *a, *b );
		return ( line_side( p, *a, *b )*inward > 0 ) && ( line_side( p, *b, *c )*inward > 0 ) && ( line_side( p, *c, *a )*inward > 0 );
	}

};

/*
// from here http://geomalgorithms.com/a01-_area.html#area2D_Triangle%28%29

inline double isLeft( const Vec2d& P0, const Vec2d& P1, const Vec2d& P2 ){
	Vec2d d10; d10.set_sub( P1, P0 );
	Vec2d d20; d20.set_sub( P2, P0 );
	return d10.cross( d20 );
	//return d10.x * d20.y - d20.x * d10.y;
    //return ( (P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y) );
}

inline double area2D_Triangle( const Vec2d& V0, const Vec2d& V1, const Vec2d& V2 ){
    return isLeft( V0, V1, V2 ) * 0.5;
}
*/

struct Ray2d{
    Vec2d p0;  //
    Vec2d dir; // normalized
    double l;  // length

    void getEnd(Vec2d& p1){
        p1.set_add_mul( p0, dir, l );
    }
    void fromPoints(const Vec2d& p0_, const Vec2d& p1){
        p0=p0_; dir.set_sub( p1, p0 ); l=dir.normalize();
    }
};

struct Circle2d{
    Vec2d  p0;
    double r;

    inline double pointDist(Vec2d p){
        p.sub(p0);
        return p.norm()-r;
    }

    inline double pointDist2(Vec2d p){
        p.sub(p0);
        return p.norm2()-r*r;
    }

    inline void fromCenterAndPoint(const Vec2d& p0_, const Vec2d& p1_){
        p0=p0_;
        r = (p0_-p1_).norm();
    }

    inline void from2points(const Vec2d& p0_, const Vec2d& p1_){
        p0 = (p0_+p1_)*0.5;
        r  = (p0_-p1_).norm()*0.5;
    }

    inline void from3points(const Vec2d& p1, const Vec2d& p2, const Vec2d& p3 ){

        // https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/

        Vec2d d12 = p1-p2;
        Vec2d d23 = p2-p3;
        Vec2d d31 = p3-p1;
        double r2_1 = p1.norm2();
        double r2_2 = p2.norm2();
        double r2_3 = p3.norm2();

        //denom = 0.5/( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) );
        //r     = sqrt( ( (x1-x2)**2 + (y1-y2)**2) * ((x1-x3)**2 + (y1-y3)**2) * ((x2-x3)**2 + (y2-y3)**2) ) * denom;
        //p0.x  = ( (x1**2 + y1**2)*(y2-y3)   +   (x2**2 + y2**2)*(y3-y1)   +   (x3**2+y3**2)*(y1-y2)      ) * denom;
        //p0.y  = ( (x1**2 + y1**2)*(x3-x2)   +   (x2**2 + y2**2)*(x1-x3)   +   (x3**2+y3**2)*(x2-x1)      ) * denom;

        double denom = 0.5/( p1.x*d23.y + p2.x*d31.y + p3.x*d12.y );
        //if(denom<0)denom=-denom;
        r     = sqrt( d12.norm2()  *   d23.norm2()  *   d31.norm2() ) *  denom;
        if(r<0)r=-r;
        p0.x  = (     r2_1*d23.y   +   r2_2*d31.y   +   r2_3*d12.y  ) *  denom;
        p0.y  = (     r2_1*d23.x   +   r2_2*d31.x   +   r2_3*d12.x  ) * -denom;

    }

    inline void fromCorner( const Vec2d& pc, const Vec2d& d1, const Vec2d& d2, double r_ ){
        r = r_;
        Vec2d d    = d1 + d2;
        double c   = d.dot(d1);
        //double l   = d.norm();
        double l2  = d.norm2();
        d.mul( r/sqrt( (l2-c*c) )  );
        //d.mul( r/(c*l) );
        p0 = pc + d;
    }
};

struct Arc2d{
    //Circle2d* circ;
    int icirc;
    //double angs[2];
    double ang0;
    double dang;

    inline double getEndAngle(){ return ang0+dang; };

    inline Vec2d getDir  ( double f )const                   { Vec2d d; d.fromAngle(ang0+(dang*f)); return d; };
    inline Vec2d getPoint( double f, const Circle2d* c )const{ Vec2d p= getDir(f); c+=icirc; p.mul(c->r); p.add(c->p0); return p; };

    inline void updateStart(const Vec2d& d ){ double a = d.toAngle(); dang+=ang0-a; ang0=a; };
    inline void updateEnd  (const Vec2d& d ){ double a = d.toAngle(); dang+=a-ang0;         };

    inline void from2vecs( const Vec2d d0, Vec2d d1 ){
        d1.udiv_cmplx(d0);
        ang0 = d0.toAngle();
        dang = d1.toAngle();
    }

    inline void fromCenter2points( const Vec2d& pc, const Vec2d& p1, const Vec2d& p2 ){
        //Vec2d d0,d1;
        //d0 = p1-pc;
        //d1 = p2-pc;
        //d1.udiv_cmplx(d0);
        //ang0 = d0.toAngle();
        //dang = d1.toAngle();
        //if(angs[0]>angs[1]) _swap(angs[0],angs[1]);
        from2vecs( p1-pc, p2-pc );
    }

    inline void fromCorner( Vec2d d1, Vec2d d2 ){
        Vec2d dc = d1+d2;
        d1.perp();
        d2.perp();
        d1.mul( -dc.dot(d1) );
        d2.mul( -dc.dot(d2) );
        from2vecs( d1, d2 );
        //d2.udiv_cmplx(d1);
        //ang0 = d1.toAngle();
        //dang = d2.toAngle();
        //angs[0] = atan2(-d1.x,d1.y);
        //angs[1] = atan2(-d2.x,d2.y);
    }

};

#endif


