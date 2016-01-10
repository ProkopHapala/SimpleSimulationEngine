
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

template <class TYPE>
inline bool pointInRect( const Vec2TYPE<TYPE>& p, TYPE x0, TYPE y0, TYPE x1, TYPE y1 ){ return ( p.x > x0 ) && ( p.x < x1 ) && ( p.y > y0 ) && ( p.y < y1 ); };

template <class TYPE>
inline bool pointInRect( const Vec2TYPE<TYPE>& p, const Vec2TYPE<TYPE>& vmax, const Vec2TYPE<TYPE>& vmin ){ return ( p.x > vmin.x ) && ( p.x < vmax.x ) && ( p.y > vmin.y ) && ( p.y < vmax.y ); };

template <class TYPE>
inline bool pointInCircle( const Vec2TYPE<TYPE>& p, const Vec2TYPE<TYPE>& center, TYPE rmax ){
	Vec2TYPE<TYPE> d;
	d.set_sub( p, center );
	TYPE r2 = d.norm2();
	return (r2 < ( rmax * rmax ) );
}

template <class TYPE>
inline bool pointInConvexPolygon( const Vec2TYPE<TYPE>& p, TYPE * buf, int n ){
    //printf( " point %f %f \n", p.x, p.y );
	for ( int i=0; i<n; i++ ){
		TYPE cv = p.dot( *(Vec2d*)buf );
        //printf( " dir %f %f   cv %f %f   \n", (*(Vec2d*)buf).x, (*(Vec2d*)buf).y,    cv, buf[2] );
		if( cv < buf[2] ) return false;
		buf += 3;
	}
	//printf( " returning true \n" );
	return true;
}

template <class TYPE>
inline bool pointInConvexParalelogram( const Vec2TYPE<TYPE>& p, TYPE * buf, int n ){
	for ( int i=0; i<n; i++ ){
		TYPE cv = p.dot( *buf );
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
		struct{ double x0,y0,x1,y1; };
		struct{ Vec2d a,b; };
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

	inline bool pointIn( const Vec2d& p ){ return pointInRect( p, x0, y0, x1, y1 );    		       }

};


/////////////////////////
//   CLASS :   Line2d
//////////////////////////

/*
inline double line_equation    ( const Vec2d& a, const Vec2d& b, double& A, double& B, double& C ){ B = a.x - b.x; A = b.y - a.y; C = A*a.x + B*a.y; }
inline double line_side        ( const Vec2d& p, double A, double B, double C   ){ return   A*p.x + B*p.y - C ;                    }
inline double dist_line_point  ( const Vec2d& p, double A, double B, double C   ){ return ( A*p.x + B*p.y + C )/sqrt( A*A + B*B ); }
inline double line_side        ( const Vec2d& p, const Vec2d& a, const Vec2d& b ){ double A,B,C; line_equation( a, b, A,B,C ); return line_side( p, A,B,C );  }
*/


inline double along_unitary( const Vec2d& a, const Vec2d& ab_hat, const Vec2d& p ){ Vec2d ap; ap.set_sub( p, a ); return ab_hat.dot(ap); }
inline double along        ( const Vec2d& a, const Vec2d& b,      const Vec2d& p ){
	Vec2d ab; ab.set_sub( b, a ); double rab = ab.norm2();
	Vec2d ap; ap.set_sub( p, a ); double ca  = ab.dot(ap);
	return ca / sqrt( rab );
}

class Line2d{
	public:
	union{
		struct{ double a,b,c;        };
		struct{ Vec2d dir; double l; };
		double array[3];
	};

	inline void   set       ( double a_, double b_, double c_ ){ a=a_; b=b_; c=c_;   };
	inline void   set       ( const Vec2d& A, const Vec2d& B  ){ b = A.x - B.x; a = B.y - A.y; c = a*A.x + b*A.y; }
	inline void   normalize (                                 ){ double ir = sqrt( a*a + b*b ); a*=ir; b*=ir; c*=ir; }
	inline double dist      ( const Vec2d& p                  ) const { return a*p.x + b*p.y - c;            };
};

inline double line_side     ( const Vec2d& p, const Vec2d& a, const Vec2d& b ){ Line2d l; l.set( a, b ); return l.dist( p );  }

inline void intersection( const Line2d& l1, const Line2d& l2, Vec2d& p ){
	double idet = 1/( l1.a * l2.b - l1.b * l2.a );
	p.x = ( l1.c * l2.b - l1.b * l2.c ) * idet;
	p.y = ( l1.a * l2.c - l1.c * l2.a ) * idet;
};

inline double intersection_t (  const Vec2d& l1a, const Vec2d& l1b, const Vec2d& l2a, const Vec2d& l2b ) {
    double dx1 = l1b.x - l1a.x;     double dy1 = l1b.y - l1a.y;
    double dx2 = l2b.x - l2a.x;     double dy2 = l2b.y - l2a.y;
	double x12 = ( l1a.x - l2a.x );
	double y12 = ( l1a.y - l2a.y );
	double invD = 1/( -dx2 * dy1 + dx1 * dy2 );
    return ( dx2 * y12 - dy2 * x12 ) * invD;
}

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


/////////////////////////
//   CLASS :   Triangle2d
//////////////////////////

class Triangle2d{
	public:
	int id;
	Vec2d *a,*b,*c;
	inline      Triangle2d( int id_ ){ id=id_;  };
	inline      Triangle2d( int id_, Vec2d* a_, Vec2d* b_, Vec2d* c_ ){ id=id_; a=a_; b=b_; c=c_; };
	inline void setPoints ( Vec2d* a_, Vec2d* b_, Vec2d* c_ ){ a=a_; b=b_; c=c_; };
	//inline void setPoints ( const Vec2d& a_, const Vec2d& b_, const Vec2d& c_ ){ a=&a_; b=&b_; c=&c_; };

	//inline bool pointIn( const Vec2d& p ){ return ( line_side( p, *a, *b ) < 0 ); }
	inline bool pointIn( const Vec2d& p ){
		double inward = line_side( *c, *a, *b );
		return ( line_side( p, *a, *b )*inward > 0 ) && ( line_side( p, *b, *c )*inward > 0 ) && ( line_side( p, *c, *a )*inward > 0 );
	}

};

/////////////////////////
//   CLASS :   Convex2d
//////////////////////////

class Convex2d{
	public:
	int id;
	int n;              // number of lines or corners ( it is the same )
	Vec2d  * corners;
	Line2d * lines;

	inline void update_lines( ){
		int nm1 = n - 1;
		for ( int i=0; i<nm1; i++ ){
			//lines[i].set( corners[i], corners[i+1] );
			lines[i].set( corners[i+1], corners[i] );
		}
		//lines[nm1].set( corners[nm1], corners[0] );
		lines[nm1].set( corners[0], corners[nm1] );
	}

	inline bool pointIn( const Vec2d& p ){
		return pointInConvexPolygon<double>( p, (double*)lines, n );
	}

    inline void boundingBox( Rect2d * rect ){
        //double a,b;
        double xmin,xmax,ymin,ymax;
        xmin=xmax=corners[0].x;
        ymin=ymax=corners[0].y;
        for ( int i=0; i<n; i++ ){
            double x = corners[0].x;
            double y = corners[0].y;
            if( x<xmin ) xmin = x;
            if( x>xmax ) xmax = x;
            if( y<ymin ) ymin = y;
            if( y>ymax ) ymax = y;
        }
        rect->set( xmin, xmax, ymin, ymax );
    }

	Convex2d( int n_ ){
	    id      = 0;
		n       = n_;
		corners = new Vec2d [ n ];
		lines   = new Line2d[ n ];
	}

	~Convex2d( ){
		delete corners;
		delete lines;
	}

};


#endif


