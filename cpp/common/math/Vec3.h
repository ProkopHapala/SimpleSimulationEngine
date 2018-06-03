
#ifndef  Vec3_h
#define  Vec3_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"

//template <class TYPE,class VEC>
template <class TYPE>
class Vec3TYPE{
	using VEC  = Vec3TYPE<TYPE>;
	using VEC2 = Vec2TYPE<TYPE>;
	public:
	union{
		struct{ TYPE x,y,z; };
		struct{ TYPE a,b,c; };
		TYPE array[3];
	};

	// ===== methods

	// Automatic conversion (works) but would be problematic
	//inline operator Vec3TYPE<float >()const{ return (Vec3TYPE<float >){(float)x,(float)y,(float)z}; }
	//inline operator Vec3TYPE<double>()const{ return (Vec3TYPE<double>){(double)x,(double)y,(double)z}; }
	//inline operator Vec3TYPE<int   >()const{ return (Vec3TYPE<int   >){(int)x,(int)y,(int)z}; }

	// Explicit conversion
    inline explicit operator Vec3TYPE<float >()const{ return (Vec3TYPE<float >){(float)x,(float)y,(float)z}; }
	inline explicit operator Vec3TYPE<double>()const{ return (Vec3TYPE<double>){(double)x,(double)y,(double)z}; }
	inline explicit operator Vec3TYPE<int   >()const{ return (Vec3TYPE<int   >){(int)x,(int)y,(int)z}; }

	//inline operator (const char*)()const{ return (; }

	//inline Vec3TYPE<double> toDouble()const{ return (Vec3TYPE<double>){ (double)x,(double)y,(double)z}; }
	//inline Vec3TYPE<float > toFloat ()const{ return (Vec3TYPE<float >){ (float)x, (double)y,(double)z}; }
	//inline Vec3TYPE<int >   toInt   ()const{ return (Vec3TYPE<int   >){ (int)x,      (int)y,   (int)z}; }

	// swizzles
	inline VEC2 xy() const { return {x,y}; };
	inline VEC2 xz() const { return {x,z}; };
	inline VEC2 yz() const { return {x,y}; };
    inline VEC2 yx() const { return {x,y}; };
	inline VEC2 zx() const { return {x,z}; };
	inline VEC2 zy() const { return {x,y}; };
    inline VEC xzy() const { return {x,z,y}; };
	inline VEC yxz() const { return {y,x,z}; };
	inline VEC yzx() const { return {y,z,x}; };
	inline VEC zxy() const { return {z,x,y}; };
	inline VEC zyx() const { return {z,y,x}; };

	inline void set( TYPE f                    ) { x=f;   y=f;   z=f;   };
    inline void set( TYPE fx, TYPE fy, TYPE fz ) { x=fx;  y=fy;  z=fz;  };
    inline void set( const VEC& v              ) { x=v.x; y=v.y; z=v.z; };
	inline void set( TYPE* arr                 ) { x=arr[0]; y=arr[1]; z=arr[2]; };

    inline void get( TYPE& fx, TYPE& fy, TYPE& fz ) { fx=x;  fy=y;  fz=z;           };
	inline void get( TYPE* arr                    ) { arr[0]=x; arr[1]=y; arr[2]=z; };

    inline void add( TYPE f ) { x+=f; y+=f; z+=f; };
    inline void mul( TYPE f ) { x*=f; y*=f; z*=f; };

    inline void add( const VEC&  v ) { x+=v.x; y+=v.y; z+=v.z; };
    inline void sub( const VEC&  v ) { x-=v.x; y-=v.y; z-=v.z; };
    inline void mul( const VEC&  v ) { x*=v.x; y*=v.y; z*=v.z; };
    inline void div( const VEC&  v ) { x/=v.x; y/=v.y; z/=v.z; };

    inline void set_inv( const VEC&  v ) { x=1/v.x; y=1/v.y; z=1/v.z; };
    inline VEC  get_inv()                { VEC o; o.x=1/x; o.y=1/y; o.z=1/z; return o; };

    inline void add( TYPE fx, TYPE fy, TYPE fz ) { x+=fx; y+=fy; z+=fz; };
    inline void sub( TYPE fx, TYPE fy, TYPE fz ) { x-=fx; y-=fy; z-=fz; };
    inline void mul( TYPE fx, TYPE fy, TYPE fz ) { x*=fx; y*=fy; z*=fz; };
    inline void div( TYPE fx, TYPE fy, TYPE fz ) { x/=fx; y/=fy; z/=fz; };

	inline void set_add( const VEC& a, TYPE f ){ x=a.x+f; y=a.y+f; z=a.z+f; };
	inline void set_mul( const VEC& a, TYPE f ){ x=a.x*f; y=a.y*f; z=a.z*f; };
	inline void set_mul( const VEC& a, const VEC& b, TYPE f ){ x=a.x*b.x*f; y=a.y*b.y*f; z=a.z*b.z*f; };

	inline void set_add( const VEC& a, const VEC& b ){ x=a.x+b.x; y=a.y+b.y; z=a.z+b.z; };
	inline void set_sub( const VEC& a, const VEC& b ){ x=a.x-b.x; y=a.y-b.y; z=a.z-b.z; };
	inline void set_mul( const VEC& a, const VEC& b ){ x=a.x*b.x; y=a.y*b.y; z=a.z*b.z; };
	inline void set_div( const VEC& a, const VEC& b ){ x=a.x/b.x; y=a.y/b.y; z=a.z/b.z; };

	inline void add_mul( const VEC& a, TYPE f                ){ x+=a.x*f;     y+=a.y*f;     z+=a.z*f;   };
	inline void add_mul( const VEC& a, const VEC& b          ){ x+=a.x*b.x;   y+=a.y*b.y;   z+=a.z*b.z; };
	inline void sub_mul( const VEC& a, const VEC& b          ){ x-=a.x*b.x;   y-=a.y*b.y;   z-=a.z*b.z; };
	inline void add_mul( const VEC& a, const VEC& b, TYPE f  ){ x+=a.x*b.x*f; y+=a.y*b.y*f; z+=a.z*b.z*f;   };


	inline void set_add_mul( const VEC& a, const VEC& b, TYPE f ){ x= a.x + f*b.x;     y= a.y + f*b.y;     z= a.z + f*b.z;  };


	inline void set_lincomb( TYPE fa, const VEC& a, TYPE fb, const VEC& b ){ x = fa*a.x + fb*b.x;  y = fa*a.y + fb*b.y;  z = fa*a.z + fb*b.z; };
	inline void add_lincomb( TYPE fa, const VEC& a, TYPE fb, const VEC& b ){ x+= fa*a.x + fb*b.x;  y+= fa*a.y + fb*b.y;  z+= fa*a.z + fb*b.z; };

	inline void set_lincomb( TYPE fa, TYPE fb, TYPE fc, const VEC& a, const VEC& b, const VEC& c ){ x = fa*a.x + fb*b.x + fc*c.x;  y = fa*a.y + fb*b.y + fc*c.y;  z = fa*a.z + fb*b.z + fc*c.z; };
	inline void add_lincomb( TYPE fa, TYPE fb, TYPE fc, const VEC& a, const VEC& b, const VEC& c ){ x+= fa*a.x + fb*b.x + fc*c.x;  y+= fa*a.y + fb*b.y + fc*c.y;  z+= fa*a.z + fb*b.z + fc*c.z; };

    inline void set_lincomb( const VEC& fs, const VEC& a, const VEC& b, const VEC& c ){ x = fs.a*a.x + fs.b*b.x + fs.c*c.x;  y = fs.a*a.y + fs.b*b.y + fs.c*c.y;  z = fs.a*a.z + fs.b*b.z + fs.c*c.z; };
	inline void add_lincomb( const VEC& fs, const VEC& a, const VEC& b, const VEC& c ){ x+= fs.a*a.x + fs.b*b.x + fs.c*c.x;  y+= fs.a*a.y + fs.b*b.y + fs.c*c.y;  z+= fs.a*a.z + fs.b*b.z + fs.c*c.z; };


    inline void set_cross( const VEC& a, const VEC& b ){ x =a.y*b.z-a.z*b.y; y =a.z*b.x-a.x*b.z; z =a.x*b.y-a.y*b.x; };
	inline void add_cross( const VEC& a, const VEC& b ){ x+=a.y*b.z-a.z*b.y; y+=a.z*b.x-a.x*b.z; z+=a.x*b.y-a.y*b.x; };

	TYPE makeOrthoU( const VEC& a ){ TYPE c = dot(a);          add_mul(a, -c); return c; }
	TYPE makeOrtho ( const VEC& a ){ TYPE c = dot(a)/a.norm(); add_mul(a, -c); return c; }

    inline VEC operator+ ( TYPE f   ) const { VEC vo; vo.x=x+f; vo.y=y+f; vo.z=z+f; return vo; };
    inline VEC operator* ( TYPE f   ) const { VEC vo; vo.x=x*f; vo.y=y*f; vo.z=z*f; return vo; };

    inline VEC operator+ ( const VEC& vi ) const { VEC vo; vo.x=x+vi.x; vo.y=y+vi.y; vo.z=z+vi.z; return vo; };
    inline VEC operator- ( const VEC& vi ) const { VEC vo; vo.x=x-vi.x; vo.y=y-vi.y; vo.z=z-vi.z; return vo; };
    inline VEC operator* ( const VEC& vi ) const { VEC vo; vo.x=x*vi.x; vo.y=y*vi.y; vo.z=z*vi.z; return vo; };
    inline VEC operator/ ( const VEC& vi ) const { VEC vo; vo.x=x/vi.x; vo.y=y/vi.y; vo.z=z/vi.z; return vo; };

	inline TYPE dot  ( const VEC& a ) const { return x*a.x + y*a.y + z*a.z;  };
	inline TYPE norm2(              ) const { return x*x + y*y + z*z;        };
	inline TYPE norm ( ) const { return  sqrt( x*x + y*y + z*z ); };
    inline TYPE normalize() {
		TYPE norm  = sqrt( x*x + y*y + z*z );
		TYPE inVnorm = 1.0d/norm;
		x *= inVnorm;    y *= inVnorm;    z *= inVnorm;
		return norm;
    };

    inline VEC getOrtho( VEC& up ) const {
        up.makeOrthoU(*this); up.normalize();
        VEC out; out.set_cross(*this,up);
        return out;
	}

	inline void getSomeOrtho( VEC& v1, VEC& v2 ) const {
		if(x<y){
//			x : y*vz - z*vy;
//			y : z*vx - x*vz;
//			z : x*vy - y*vx;
//			x : y*0 - z*0 ;
//			y : z*1 - x*0 ;
//			z : x*0 - y*1 ;
//			float vx = 0; float vy = z; float vz =-y;
			v1.x =  -y*y -z*z;
			v1.y =  x*y;
			v1.z =  x*z;
		}else{
//			x : y*0 - z*1;
//			y : z*0 - x*0;
//			z : x*1 - y*0;
//			float vx = -z; float vy = 0; float vz = x;
			v1.x =  y*x;
			v1.y =  -z*z -x*x;
			v1.z =  y*z;
		}
		v2.x = y*v1.z - z*v1.y;
		v2.y = z*v1.x - x*v1.z;
		v2.z = x*v1.y - y*v1.x;
	}

	// Rodrigues rotation formula: v' = cosa*v + sina*(uaxis X v) + (1-cosa)*(uaxis . v)*uaxis
	inline void rotate( TYPE angle, const VEC& axis  ){
		VEC uaxis;
		uaxis.set( axis * axis.norm() );
		TYPE ca   = cos(angle);
		TYPE sa   = sin(angle);
 		rotate_csa( ca, sa, uaxis );
	};

	inline void rotate_csa( TYPE ca, TYPE sa, const VEC& uaxis ){
		TYPE cu = (1-ca)*dot(uaxis);
		TYPE utx  = uaxis.y*z - uaxis.z*y;
		TYPE uty  = uaxis.z*x - uaxis.x*z;
		TYPE utz  = uaxis.x*y - uaxis.y*x;
		TYPE x_ = ca*x + sa*utx + cu*uaxis.x;
		TYPE y_ = ca*y + sa*uty + cu*uaxis.y;
		       z  = ca*z + sa*utz + cu*uaxis.z;
		x = x_; y = y_;
	};

	inline void rotateTo( const VEC& rot0, double coef ){
        //rot.add_mul( rot0, coef ); rot.normalize();
        VEC ax; ax.set_cross( *this, rot0 );
        double sa2 = ax.norm2();
        if( sa2 < coef*coef ){
            ax.mul( 1/sqrt(sa2) ); // this is safe if coef is large enough
            double ca = sqrt( 1-coef*coef );
            rotate_csa( ca, coef, ax );
        }else{
            set(rot0);
        }
    }

    inline void getInPlaneRotation( const VEC& rot0, const VEC& xhat, const VEC& yhat, double& ca, double& sa ){
        double x0 = rot0.dot(xhat);
        double y0 = rot0.dot(yhat);
        double x_ = dot(xhat);
        double y_ = dot(yhat);
        // http://mathworld.wolfram.com/ComplexDivision.html
        double renorm = 1.0/sqrt( (x0*x0 + y0*y0)*(x_*x_ + y_*y_) );
        ca = ( x0*x_ + y0*y_ ) * renorm;
        sa = ( y0*x_ - x0*y_ ) * renorm;
    }

	inline TYPE along_hat( const VEC& hat, const VEC& p ){ VEC ap; ap.set( p.x-x, p.y-y ); return hat.dot( ap ); }
	inline TYPE along    ( const VEC& b,   const VEC& p ){
		VEC ab,ap;
		ab.set( b.x - x, b.y - y, b.z - z );
		ap.set( p.x - x, p.y - y, b.z - z );
		return ab.dot(ap) / ab.norm(ab);
	}

    inline bool isLower  ( const VEC& vmax ) const { return (x<vmax.x)&&(y<vmax.y)&&(x<vmax.z); }
    inline bool isGreater( const VEC& vmin ) const { return (x>vmin.x)&&(y>vmin.y)&&(x>vmin.z); }
    inline bool isBetween( const VEC& vmin, const VEC& vmax ) const { return (x>vmin.x)&&(x<vmax.x)&&(y>vmin.y)&&(y<vmax.y)&&(z>vmin.z)&&(z<vmax.z); }

    inline void setIfLower  (const VEC& a){ if(a.x<x)x=a.x;if(a.y<y)y=a.y;if(a.z<z)z=a.z; }
    inline void setIfGreater(const VEC& a){ if(a.x>x)x=a.x;if(a.y>y)y=a.y;if(a.z>z)z=a.z; }
    //inline VEC min(VEC a){ return {fmin(x,a.x),fmin(y,a.y),fmin(z,a.z)}; };
    //inline VEC max(VEC a){ return {fmax(x,a.x),fmax(y,a.y),fmax(z,a.z)}; };
    //inline VEC set_min(VEC a,VEC b){ return {fmin(x,a.x),fmin(y,a.y),fmin(z,a.z)}; };
    //inline VEC set_max(VEC a,VEC b){ return {fmax(x,a.x),fmax(y,a.y),fmax(z,a.z)}; };

    inline TYPE dist2( const VEC& a ) const { VEC d; d.set( x-a.x, y-a.y, z-a.z ); return d.norm2(); }
    inline TYPE dist ( const VEC& a ) const { VEC d; d.set( x-a.x, y-a.y, z-a.z ); return d.norm (); }

    inline TYPE totprod(){ return x*y*z; };

    TYPE angleInPlane( const VEC& a, const VEC& b ){
        TYPE x = dot(a);
        TYPE y = dot(b);
        return atan2( y, x );
    }

    inline void fromLinearSolution( const VEC& va, const VEC& vb, const VEC& vc, const VEC& p ){
        // https://en.wikipedia.org/wiki/Cramer%27s_rule
        // 30 multiplications
        TYPE Dax = vb.y*vc.z - vb.z*vc.y;
        TYPE Day = vb.x*vc.z - vb.z*vc.x;
        TYPE Daz = vb.x*vc.y - vb.y*vc.x;
        TYPE idet = 1/( va.x*Dax - va.y*Day + va.z*Daz );
        x =  idet*( p.x*Dax - p.y*Day + p.z*Daz );
        y = -idet*( p.x*(va.y*vc.z - va.z*vc.y) - p.y*(va.x*vc.z - va.z*vc.x) + p.z*(va.x*vc.y - va.y*vc.x) );
        z =  idet*( p.x*(va.y*vb.z - va.z*vb.y) - p.y*(va.x*vb.z - va.z*vb.x) + p.z*(va.x*vb.y - va.y*vb.x) );
    }

};

template<typename VEC> inline VEC cross( VEC a, VEC b ){ return (VEC){ a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x }; }
template<typename VEC> inline VEC add  ( VEC a, VEC b ){ return (VEC){ a.x+b.x, a.z+b.z, a.z+b.z }; }

using Vec3i = Vec3TYPE<int>;
using Vec3f = Vec3TYPE<float>;
using Vec3d = Vec3TYPE<double>;

static constexpr Vec3d Vec3dZero = (Vec3d){0.0d,0.0d,0.0d};
static constexpr Vec3d Vec3dX    = (Vec3d){1.0d,0.0d,0.0d};
static constexpr Vec3d Vec3dY    = (Vec3d){0.0d,1.0d,0.0d};
static constexpr Vec3d Vec3dZ    = (Vec3d){0.0d,0.0d,1.0d};

static constexpr Vec3f Vec3fZero = (Vec3f){0.0f,0.0f,0.0f};
static constexpr Vec3f Vec3fX    = (Vec3f){1.0f,0.0f,0.0f};
static constexpr Vec3f Vec3fY    = (Vec3f){0.0f,1.0f,0.0f};
static constexpr Vec3f Vec3fZ    = (Vec3f){0.0f,0.0f,1.0f};

static constexpr Vec3i Vec3iZero = (Vec3i){0,0,0};
static constexpr Vec3i Vec3iX    = (Vec3i){1,0,0};
static constexpr Vec3i Vec3iY    = (Vec3i){0,1,0};
static constexpr Vec3i Vec3iZ    = (Vec3i){0,0,1};






inline uint64_t scalar_id  ( const Vec3i& v){ return ( v.x | (((uint64_t)v.y)<<16) | (((uint64_t)v.z)<<32) ); }
inline Vec3i    from_id    ( uint64_t id   ){
    Vec3i vi;
    vi.x=( id & 0xFFFF ); id>>16;
    vi.y=( id & 0xFFFF ); id>>16;
    vi.z=( id & 0xFFFF );
    return vi;
}


inline void print(Vec3d p){printf("(%.16g,%.16g,%.16g)", p.x,p.y,p.z);};
inline void print(Vec3f p){printf("(%.8g,%.8g,%.8g)", p.x,p.y,p.z);};
//inline void print(Vec3d p){printf("(%lg,%lg,%lg)", p.x,p.y,p.z);};
//inline void print(Vec3f p){printf("(%g,%g,%g)", p.x,p.y,p.z);};
inline void print(Vec3i p){printf("(%i,%i,%i)", p.x,p.y,p.z);};

inline void convert( const Vec3f& from, Vec3d& to ){ to.x=from.x;        to.y=from.y;        to.z=from.z; };
inline void convert( const Vec3d& from, Vec3f& to ){ to.x=(float)from.x; to.y=(float)from.y; to.z=(float)from.z; };

inline Vec3f toFloat( const Vec3d& from){ return (Vec3f){(float)from.x,(float)from.y,(float)from.z}; }

#endif



