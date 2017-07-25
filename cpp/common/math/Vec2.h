
#ifndef  Vec2_h
#define  Vec2_h

#include "fastmath.h"

//template <class TYPE,class VEC>
template <class TYPE>
class Vec2TYPE{
	using VEC = Vec2TYPE<TYPE>;
	public:
	union{
		struct{ TYPE x,y; };
		struct{ TYPE a,b; };
		TYPE array[2];
	};

	// ===== methods
	inline void set( TYPE f            ) { x=f;   y=f;   };
    inline void set( TYPE fx, TYPE fy  ) { x=fx;  y=fy;  };
    inline void set( const VEC& v      ) { x=v.x; y=v.y; };
	inline void set( TYPE* arr         ) { x=arr[0]; y=arr[1]; };

	//inline Vec2TYPE(){};
	//inline Vec2TYPE( TYPE f            ){ x=f;   y=f; };
	//inline Vec2TYPE( TYPE fx, TYPE fy  ){ x=fx;  y=fy; };
	//inline Vec2TYPE( const VEC& v      ){ x=v.x; y=v.y; };
	//inline Vec2TYPE( TYPE* arr         ){ x=arr[0]; y=arr[1]; };

    inline void get( TYPE& fx, TYPE& fy ) { fx=x;  fy=y;        };
	inline void get( TYPE* arr                    ) { arr[0]=x; arr[1]=y; };

    inline void add( TYPE f ) { x+=f; y+=f; };
    inline void mul( TYPE f ) { x*=f; y*=f; };

    inline void add( const VEC&  v ) { x+=v.x; y+=v.y; };
    inline void sub( const VEC&  v ) { x-=v.x; y-=v.y; };
    inline void mul( const VEC&  v ) { x*=v.x; y*=v.y; };
    inline void div( const VEC&  v ) { x/=v.x; y/=v.y; };

    inline void add( TYPE fx, TYPE fy) { x+=fx; y+=fy; };
    inline void sub( TYPE fx, TYPE fy) { x-=fx; y-=fy; };
    inline void mul( TYPE fx, TYPE fy ) { x*=fx; y*=fy; };
    inline void div( TYPE fx, TYPE fy ) { x/=fx; y/=fy; };

    inline void set_inv( const VEC&  v ) { x=1/v.x; y=1/v.y; };

	inline void set_add( const VEC& a, TYPE f ){ x=a.x+f; y=a.y+f; };
	inline void set_mul( const VEC& a, TYPE f ){ x=a.x*f; y=a.y*f; };
	inline void set_mul( const VEC& a, const VEC& b, TYPE f ){ x=a.x*b.x*f; y=a.y*b.y*f; };

	inline void set_add( const VEC& a, const VEC& b ){ x=a.x+b.x; y=a.y+b.y; };
	inline void set_sub( const VEC& a, const VEC& b ){ x=a.x-b.x; y=a.y-b.y; };
	inline void set_mul( const VEC& a, const VEC& b ){ x=a.x*b.x; y=a.y*b.y; };
	inline void set_div( const VEC& a, const VEC& b ){ x=a.x/b.x; y=a.y/b.y; };

	inline void add_mul( const VEC& a, TYPE f                ){ x+=a.x*f;     y+=a.y*f;     };
	inline void add_mul( const VEC& a, const VEC& b          ){ x+=a.x*b.x;   y+=a.y*b.y;   };
	inline void sub_mul( const VEC& a, const VEC& b          ){ x-=a.x*b.x;   y-=a.y*b.y;   };
	inline void add_mul( const VEC& a, const VEC& b, TYPE f  ){ x+=a.x*b.x*f; y+=a.y*b.y*f; };

	inline void set_add_mul( const VEC& a, const VEC& b, TYPE f ){ x= a.x + f*b.x;     y= a.y + f*b.y; };

	inline void set_lincomb( TYPE fa, const VEC& a, TYPE fb, const VEC& b ){ x = fa*a.x + fb*b.x;  y = fa*a.y + fb*b.y; };
	inline void add_lincomb( TYPE fa, const VEC& a, TYPE fb, const VEC& b ){ x+= fa*a.x + fb*b.x;  y+= fa*a.y + fb*b.y; };

	inline void set_lincomb( TYPE fa, TYPE fb, TYPE fc, const VEC& a, const VEC& b, const VEC& c ){ x = fa*a.x + fb*b.x + fc*c.x;  y = fa*a.y + fb*b.y + fc*c.y; };
	inline void add_lincomb( TYPE fa, TYPE fb, TYPE fc, const VEC& a, const VEC& b, const VEC& c ){ x+= fa*a.x + fb*b.x + fc*c.x;  y+= fa*a.y + fb*b.y + fc*c.y; };


    inline VEC operator+ ( TYPE f   ) const { VEC vo; vo.x=x+f; vo.y=y+f; return vo; };
    inline VEC operator* ( TYPE f   ) const { VEC vo; vo.x=x*f; vo.y=y*f; return vo; };

    inline VEC operator+ ( const VEC& vi ) const { VEC vo; vo.x=x+vi.x; vo.y=y+vi.y; return vo; };
    inline VEC operator- ( const VEC& vi ) const { VEC vo; vo.x=x-vi.x; vo.y=y-vi.y; return vo; };
    inline VEC operator* ( const VEC& vi ) const { VEC vo; vo.x=x*vi.x; vo.y=y*vi.y; return vo; };
    inline VEC operator/ ( const VEC& vi ) const { VEC vo; vo.x=x/vi.x; vo.y=y/vi.y; return vo; };

	inline TYPE dot      ( const VEC& a ) const { return x*a.x + y*a.y; };
	inline TYPE dot_perp ( const VEC& a ) const { return y*a.x - x*a.y; };
	inline TYPE norm2(              ) const { return x*x + y*y;     };

	inline TYPE norm ( ) const { return  sqrt( x*x + y*y ); };
    inline TYPE normalize() {
		TYPE norm  = sqrt( x*x + y*y );
		TYPE inVnorm = 1.0d/norm;
		x *= inVnorm;    y *= inVnorm;
		return norm;
    };

    inline TYPE dist2( const VEC& a) const { VEC d; TYPE dx = x-a.x; TYPE dy = y-a.y; return dx*dx + dy*dy; }
    inline TYPE dist ( const VEC& a) const { return sqrt( dist2(a) ); }

	inline void set_perp( const VEC& a )       { x=-a.y; y=a.x; }
	inline double cross ( const VEC& a ) const { return x*a.y - y*a.x; };

	inline void     mul_cmplx (               const VEC& b ){                            double x_ =    x*b.x -   y*b.y;         y =    y*b.x +   x*b.y;       x=x_;  }
	inline void pre_mul_cmplx ( const VEC& a               ){                            double x_ =  a.x*  x - a.y*  y;         y =  a.y*  x + a.x*  y;       x=x_;  }
	inline void set_mul_cmplx ( const VEC& a, const VEC& b ){                            double x_ =  a.x*b.x - a.y*b.y;         y =  a.y*b.x + a.x*b.y;       x=x_;  }
	inline void set_udiv_cmplx( const VEC& a, const VEC& b ){                            double x_ =  a.x*b.x + a.y*b.y;         y =  a.y*b.x - a.x*b.y;       x=x_;  }
	inline void set_div_cmplx ( const VEC& a, const VEC& b ){ double ir2 = 1/b.norm2();  double x_ = (a.x*b.x + a.y*b.y)*ir2;    y = (a.y*b.x - a.x*b.y)*ir2;  x=x_;  }

	inline void fromAngle        ( double phi ){	x = cos( phi );	y = sin( phi );	    }
	inline void fromAngle_taylor2( double phi ){	sincos_taylor2<TYPE>( phi, y, x );	}

	inline void rotate( double phi ){
		double bx = cos( phi );   		  double by = sin( phi );
		double x_ =    x*bx -   y*by;         y =    y*bx +   x*by;       x=x_;
	}

	inline void rotate_taylor2( double phi ){
		double bx,by;  sincos_taylor2<TYPE>( phi, by, bx );
		double x_ =    x*bx -   y*by;         y =    y*bx +   x*by;       x=x_;
	}

	inline double along_hat( const VEC& hat, const VEC& p ){ VEC ap; ap.set( p.x-x, p.y-y ); return hat.dot( ap ); }
	inline double along    ( const VEC& b,   const VEC& p ){
		VEC ab,ap;
		ab.set( b.x - x, b.y - y );
		ap.set( p.x - x, p.y - y );
		return ab.dot(ap) / ab.norm(ab);
	}

	inline TYPE angle( const VEC& a ){
		TYPE d = cdot ( a );
		TYPE c = cross( a );
		return atan2( d, c );
	}

};

using Vec2i = Vec2TYPE<int>;
using Vec2f = Vec2TYPE<float>;
using Vec2d = Vec2TYPE<double>;

inline uint64_t scalar_id  ( const Vec2i& v){ return (((uint64_t)v.a)<<32)|v.b; }
inline uint64_t symetric_id( const Vec2i& v){ if( v.a>v.b ){ return (((uint64_t)v.b)<<32)|v.a; }else{ return (((uint64_t)v.a)<<32)|v.b; }}

inline void convert( const Vec2f& from, Vec2d& to ){ to.x=from.x;        to.y=from.y;        };
inline void convert( const Vec2d& from, Vec2f& to ){ to.x=(float)from.x; to.y=(float)from.y; };

inline Vec2f toFloat( const Vec2d& from){ return (Vec2f){(float)from.x,(float)from.y}; }

#endif


