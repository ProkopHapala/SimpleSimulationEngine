
// read also:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

#ifndef  geom2D_h
#define  geom2D_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

class Vec2d;

/////////////////////////////
//  fast math functions
/////////////////////////////

inline double sin_taylor2( double a ){
	const double c3 = 1.0d/6;        
	const double c5 = 1.0d/120;
	double a2 = a*a;
	return    a * ( 1 - a2*( c3 - c5*a2 ) ); 
}

inline double cos_taylor2( double a ){
	const double c2 = 1.0d/2;        
	const double c4 = 1.0d/24;  
	double a2 = a*a;
	return    1 - a2*( c2 - c4*a2 );
}

inline void sincos_taylor2( double a, double& ca, double& sa ){
	const double c2 = 1.0d/2;       
	const double c3 = 1.0d/6;      
	const double c4 = 1.0d/24;     
	const double c5 = 1.0d/120; 
	double a2 = a*a;
	sa   = a * ( 1 - a2*( c3 - c5*a2 ) ) ; 
	ca   =       1 - a2*( c2 - c4*a2 )   ;
}

/////////////////////////
//   CLASS :   Vec2d
//////////////////////////

class Vec2d{
	public:
	union{
		struct{ double x,y; };
		struct{ double a,b; };
		double array[2];
	};

	//inline Vec2d(){};
	//inline Vec2d( double f              ) { x=f;   y=f;   };
    //inline Vec2d( double fx, double fy  ) { x=fx;  y=fy;  };
    //inline Vec2d( const  Vec2d& v       ) { x=v.x; y=v.y; };

	inline void set( double  f            ) { x=f;   y=f;   };
    inline void set( double fx, double fy ) { x=fx;  y=fy;  };
    inline void set( const Vec2d& v       ) { x=v.x; y=v.y; };

	// do not use assignement operator because it is not obious what is ment
	// also there is conflict with Mat3d: member ‘Vec2d Mat3d::<anonymous union>::<anonymous struct>::a’ with copy assignment operator not allowed in anonymous aggregate
	//inline Vec2d& operator =( const double f ) { x=f; y=f;       return *this; }; 
	//inline Vec2d& operator =( const Vec2d& v ) { x=v.x; y=v.y; return *this; }; 

    inline Vec2d& operator+=( double f ) { x+=f; y+=f; return *this; };
    inline Vec2d& operator*=( double f ) { x*=f; y*=f; return *this; };

    inline Vec2d& operator+=( const Vec2d&  v ) { x+=v.x; y+=v.y; return *this; };
    inline Vec2d& operator-=( const Vec2d&  v ) { x-=v.x; y-=v.y; return *this; };
    inline Vec2d& operator*=( const Vec2d&  v ) { x*=v.x; y*=v.y; return *this; };
    inline Vec2d& operator/=( const Vec2d&  v ) { x/=v.x; y/=v.y; return *this; };

	// This creates a new vectors? is it good?
    inline Vec2d operator+ ( double f   ) const { Vec2d vo; vo.x=x+f; vo.y=y+f; return vo; };
    inline Vec2d operator* ( double f   ) const { Vec2d vo; vo.x=x*f; vo.y=y*f; return vo; };

    inline Vec2d operator+ ( const Vec2d& vi ) const { Vec2d vo; vo.x=x+vi.x; vo.y=y+vi.y; return vo; };
    inline Vec2d operator- ( const Vec2d& vi ) const { Vec2d vo; vo.x=x-vi.x; vo.y=y-vi.y; return vo; };
    inline Vec2d operator* ( const Vec2d& vi ) const { Vec2d vo; vo.x=x*vi.x; vo.y=y*vi.y; return vo; };
    inline Vec2d operator/ ( const Vec2d& vi ) const { Vec2d vo; vo.x=x/vi.x; vo.y=y/vi.y; return vo; };

	inline void mul   ( const Vec2d& b                 ){ x*=b.x; y*=b.y;        };

	inline void set_add( const Vec2d& a, const Vec2d& b ){ x=a.x+b.x; y=a.y+b.y;    };
	inline void set_sub( const Vec2d& a, const Vec2d& b ){ x=a.x-b.x; y=a.y-b.y;    };
	inline void set_mul( const Vec2d& a, const Vec2d& b ){ x=a.x*b.x; y=a.y*b.y;    };
	inline void set_div( const Vec2d& a, const Vec2d& b ){ x=a.x/b.x; y=a.y/b.y;    };
    inline void set_fma( double       f, const Vec2d& b ){ x+=  f*b.x; y+=  f*b.y;  };

	inline double norm2( ) const { return        x*x + y*y  ; };
	inline double norm ( ) const { return  sqrt( x*x + y*y ); };
    inline double normalize() {
		double norm  = sqrt( x*x + y*y );
		double inorm = 1.0d/norm;
		x *= inorm;    y *= inorm;
		return norm;
    }

	inline void set_perp      ( const Vec2d& a ) {  x =a.y; y=-a.x; }
	inline void set_mulcomplex( const Vec2d& a, const Vec2d& b ){ double x_ = a.x*b.x - a.y*b.y;  y = a.x*b.y + a.y*b.x;  x=x_; }

	inline void set_angle   ( double angle                 ){ x = cos( angle ); y = sin( angle ); }
	inline void set_rotated ( double angle, const Vec2d& a ){ double ca = cos( angle ); double sa = sin( angle ); double x_ = a.x*ca - a.y*sa; y = a.x*sa + a.y*ca; x=x_; };
	inline void drot_taylor2( double angle                 ){ double ca, sa;  sincos_taylor2( angle, ca, sa );    double x_ =   x*ca -   y*sa; y =   x*sa +   y*ca; x=x_; };

};

inline double dot   ( const Vec2d& a, const Vec2d& b ){ return  a.x*b.x + a.y*b.y; }
inline double cross ( const Vec2d& a, const Vec2d& b ){ return  a.x*b.y - a.y*b.x; }
inline double dist2 ( const Vec2d& a, const Vec2d& b ){ double dx = b.x - a.x; double dy = b.y - a.y; return dx*dx + dy*dy; }
inline double dist  ( const Vec2d& a, const Vec2d& b ){ return sqrt( dist2( a, b ) );  }



/////////////////////////
//   CLASS :   Rect2d
/////////////////////////

inline bool pointInRect( const Vec2d& p, double x0, double y0, double x1, double y1 ){ return ( p.x > x0 ) && ( p.x < x1 ) && ( p.y > y0 ) && ( p.y < y1 ); };

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


inline double along_unitary( const Vec2d& a, const Vec2d& ab_hat, const Vec2d& p ){ Vec2d ap; ap.set_sub( p, a ); return dot(ap,ab_hat); }
inline double along        ( const Vec2d& a, const Vec2d& b,      const Vec2d& p ){
	Vec2d ab; ab.set_sub( b, a ); double rab = dot(ab,ab);
	Vec2d ap; ap.set_sub( p, a ); double ca  = dot(ap,ab);
	return ca / sqrt( rab );
}

class Lin2d{
	public:
	union{
		struct{ double a,b,c; };
		double array[3];
	};

	inline void   set       ( double a_, double b_, double c_ ){ a=a_; b=b_; c=c_;   };
	inline void   set       ( const Vec2d& A, const Vec2d& B  ){ b = A.x - B.x; a = B.y - A.y; c = a*A.x + b*A.y; }
	inline void   normalize ( ){ double ir = sqrt( a*a + b*b ); a*=ir; b*=ir; c*=ir; }
	inline double dist      ( const Vec2d& p ) const { return a*p.x + b*p.y - c;            };
};

inline double line_side     ( const Vec2d& p, const Vec2d& a, const Vec2d& b ){ Lin2d l; l.set( a, b ); return l.dist( p );  }   


inline void intersection( const Lin2d& l1, const Lin2d& l2, Vec2d& p ){  
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
    double dx1 = l1b.x - l1a.x;     double dy1 = l1b.y - l1a.y;
    double dx2 = l2b.x - l2a.x;     double dy2 = l2b.y - l2a.y;
	double x12 = ( l1a.x - l2a.x );
	double y12 = ( l1a.y - l2a.y );
	double invD = 1/( -dx2 * dy1 + dx1 * dy2 );
    double s = (-dy1 * x12 + dx1 * y12 ) * invD;
    double t = ( dx2 * y12 - dy2 * x12 ) * invD;
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


#endif


