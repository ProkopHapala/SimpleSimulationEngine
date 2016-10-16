
#ifndef  Mat3_h
#define  Mat3_h

#include <math.h>
#include <cstdlib>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"

//template <class TYPE, class VEC, class MAT>
//template <class TYPE, class VEC>
template <class TYPE>
class Mat3TYPE{
	using VEC = Vec3TYPE<TYPE>;
	using MAT = Mat3TYPE<TYPE>;
	public:
	union{
		struct{
			TYPE xx,xy,xz;
			TYPE yx,yy,yz;
			TYPE zx,zy,zz;
		};
		struct{
			TYPE ax,ay,az;
			TYPE bx,by,bz;
			TYPE cx,cy,cz;
		};
		struct{	VEC a,b,c; };
		TYPE array[9];
	};


// ====== initialization

	inline void setOne(        ){ xx=yy=zz=1; xy=xz=yx=yz=zx=zy=0; };
	inline void set   ( TYPE f ){ xx=yy=zz=f; xy=xz=yx=yz=zx=zy=0; };

	inline void set  ( const VEC& va, const VEC& vb, const VEC& vc ){ a.set(va); b.set(vb); c.set(vc); }
	inline void set  ( const MAT& M ){
		xx=M.xx; xy=M.xy; xz=M.xz;
		yx=M.yx; yy=M.yy; yz=M.yz;
		zx=M.zx; zy=M.zy; zz=M.zz;
	};


	inline void set_outer  ( const VEC& a, const VEC& b ){
		xx=a.x*b.x; xy=a.x*b.y; xz=a.x*b.z;
		yx=a.y*b.x; yy=a.y*b.y; yz=a.y*b.z;
		zx=a.z*b.x; zy=a.z*b.y; zz=a.z*b.z;
	};;

	inline VEC getColx(){ VEC out; out.x = xx; out.y = yx; out.z = zx; return out; };
    inline VEC getColy(){ VEC out; out.x = xy; out.y = yy; out.z = zy; return out; };
    inline VEC getColz(){ VEC out; out.x = xz; out.y = yz; out.z = zz; return out; };

	inline void  colx_to( VEC& out){ out.x = xx; out.y = yx; out.z = zx; };
    inline void  coly_to( VEC& out){ out.x = xy; out.y = yy; out.z = zy; };
    inline void  colz_to( VEC& out){ out.x = xz; out.y = yz; out.z = zz; };

	inline void  setColx( const VEC v ){ xx = v.x; yx = v.y; zx = v.z; };
	inline void  setColy( const VEC v ){ xy = v.x; yy = v.y; zy = v.z; };
	inline void  setColz( const VEC v ){ xz = v.x; yz = v.y; zz = v.z; };

	// Don't need this, because we use union: use representation a,b,c
	//inline VEC getRowx(){ VEC out; out.x = xx; out.y = xy; out.z = xz; return out; };
	//inline VEC getRowy(){ VEC out; out.x = yx; out.y = yy; out.z = yz; return out; };
	//inline VEC getRowz(){ VEC out; out.x = zx; out.y = zy; out.z = zz; return out; };
	//inline void rowx_to( VEC& out ){ out.x = xx; out.y = xy; out.z = xz; };
	//inline void rowy_to( VEC& out ){ out.x = yx; out.y = yy; out.z = yz; };
	//inline void rowz_to( VEC& out ){ out.x = zx; out.y = zy; out.z = zz; };
	//inline void  setRowx( const VEC& v ){ xx = v.x; xy = v.y; xz = v.z; };
	//inline void  setRowy( const VEC& v ){ yx = v.x; yy = v.y; yz = v.z; };
	//inline void  setRowz( const VEC& v ){ zx = v.x; zy = v.y; zz = v.z; };

// ====== transpose

	inline void T(){
		TYPE t1=yx; yx=xy; xy=t1;
		TYPE t2=zx; zx=xz; xz=t2;
		TYPE t3=zy; zy=yz; yz=t3;
	};

	inline void setT  ( const MAT& M ){
		xx=M.xx; xy=M.xy; xz=M.xz;
		yx=M.yx; yy=M.yy; yz=M.yz;
		zx=M.zx; zy=M.zy; zz=M.zz;
	};

	inline void setT  ( const VEC& va, const VEC& vb, const VEC& vc ){
		a.set( va.x, vb.x, vc.x );
		b.set( va.y, vb.y, vc.y );
		c.set( va.z, vb.z, vc.z );
	};

// ====== dot product with vector

	inline VEC dot( const VEC&  v ) const {
		VEC vout;
		vout.x = xx*v.x + xy*v.y + xz*v.z;
		vout.y = yx*v.x + yy*v.y + yz*v.z;
		vout.z = zx*v.x + zy*v.y + zz*v.z;
		return vout;
	}

	inline void dot_to( const VEC&  v, VEC&  vout ) const {
		vout.x = xx*v.x + xy*v.y + xz*v.z;
		vout.y = yx*v.x + yy*v.y + yz*v.z;
		vout.z = zx*v.x + zy*v.y + zz*v.z;
	};

	inline void dot_to_T( const VEC&  v, VEC&  vout ) const {
		vout.x = xx*v.x + yx*v.y + zx*v.z;
		vout.y = xy*v.x + yy*v.y + zy*v.z;
		vout.z = xz*v.x + yz*v.y + zz*v.z;
	};

// ====== matrix multiplication

	inline void set_mmul( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.xy*B.yx + A.xz*B.zx;
		xy = A.xx*B.xy + A.xy*B.yy + A.xz*B.zy;
		xz = A.xx*B.xz + A.xy*B.yz + A.xz*B.zz;
		yx = A.yx*B.xx + A.yy*B.yx + A.yz*B.zx;
		yy = A.yx*B.xy + A.yy*B.yy + A.yz*B.zy;
		yz = A.yx*B.xz + A.yy*B.yz + A.yz*B.zz;
		zx = A.zx*B.xx + A.zy*B.yx + A.zz*B.zx;
		zy = A.zx*B.xy + A.zy*B.yy + A.zz*B.zy;
		zz = A.zx*B.xz + A.zy*B.yz + A.zz*B.zz;
	};

	inline void set_mmul_NT( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.xy*B.xy + A.xz*B.xz;
		xy = A.xx*B.yx + A.xy*B.yy + A.xz*B.yz;
		xz = A.xx*B.zx + A.xy*B.zy + A.xz*B.zz;
		yx = A.yx*B.xx + A.yy*B.xy + A.yz*B.xz;
		yy = A.yx*B.yx + A.yy*B.yy + A.yz*B.yz;
		yz = A.yx*B.zx + A.yy*B.zy + A.yz*B.zz;
		zx = A.zx*B.xx + A.zy*B.xy + A.zz*B.xz;
		zy = A.zx*B.yx + A.zy*B.yy + A.zz*B.yz;
		zz = A.zx*B.zx + A.zy*B.zy + A.zz*B.zz;
	};

	inline void set_mmul_TN( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.yx*B.yx + A.zx*B.zx;
		xy = A.xx*B.xy + A.yx*B.yy + A.zx*B.zy;
		xz = A.xx*B.xz + A.yx*B.yz + A.zx*B.zz;
		yx = A.xy*B.xx + A.yy*B.yx + A.zy*B.zx;
		yy = A.xy*B.xy + A.yy*B.yy + A.zy*B.zy;
		yz = A.xy*B.xz + A.yy*B.yz + A.zy*B.zz;
		zx = A.xz*B.xx + A.yz*B.yx + A.zz*B.zx;
		zy = A.xz*B.xy + A.yz*B.yy + A.zz*B.zy;
		zz = A.xz*B.xz + A.yz*B.yz + A.zz*B.zz;
	};

	inline void set_mmul_TT( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.yx*B.xy + A.zx*B.xz;
		xy = A.xx*B.yx + A.yx*B.yy + A.zx*B.yz;
		xz = A.xx*B.zx + A.yx*B.zy + A.zx*B.zz;
		yx = A.xy*B.xx + A.yy*B.xy + A.zy*B.xz;
		yy = A.xy*B.yx + A.yy*B.yy + A.zy*B.yz;
		yz = A.xy*B.zx + A.yy*B.zy + A.zy*B.zz;
		zx = A.xz*B.xx + A.yz*B.xy + A.zz*B.xz;
		zy = A.xz*B.yx + A.yz*B.yy + A.zz*B.yz;
		zz = A.xz*B.zx + A.yz*B.zy + A.zz*B.zz;
	};

// ====== matrix solver

   inline TYPE determinant() {
        TYPE fCoxx = yy * zz - yz * zy;
        TYPE fCoyx = yz * zx - yx * zz;
        TYPE fCozx = yx * zy - yy * zx;
        TYPE fDet = xx * fCoxx + xy * fCoyx + xz * fCozx;
        return fDet;
    };

	inline void invert_to( MAT& Mout ) {
        TYPE idet = 1/determinant(); // we dont check det|M|=0
        Mout.xx = ( yy * zz - yz * zy ) * idet;
        Mout.xy = ( xz * zy - xy * zz ) * idet;
        Mout.xz = ( xy * yz - xz * yy ) * idet;
        Mout.yx = ( yz * zx - yx * zz ) * idet;
        Mout.yy = ( xx * zz - xz * zx ) * idet;
        Mout.yz = ( xz * yx - xx * yz ) * idet;
        Mout.zx = ( yx * zy - yy * zx ) * idet;
        Mout.zy = ( xy * zx - xx * zy ) * idet;
        Mout.zz = ( xx * yy - xy * yx ) * idet;
    };

    inline void adjoint_to( MAT& Mout ) {
        Mout.xx = yy * zz - yz * zy;
        Mout.xy = xz * zy - xy * zz;
        Mout.xz = xy * yz - xz * yy;
        Mout.yx = yz * zx - yx * zz;
        Mout.yy = xx * zz - xz * zx;
        Mout.yz = xz * yx - xx * yz;
        Mout.zx = yx * zy - yy * zx;
        Mout.zy = xy * zx - xx * zy;
        Mout.zz = xx * yy - xy * yx;
    };

// ======= Rotation

	inline void rotate( TYPE angle, const VEC& axis  ){
		Vec3d uaxis;
		uaxis.set( axis * axis.norm() );
		TYPE ca   = cos(angle);
		TYPE sa   = sin(angle);
 		rotate_csa( ca, sa, uaxis );
	};

	inline void rotate_csa( TYPE ca, TYPE sa, const VEC& uaxis ){
		a.rotate_csa( ca, sa, uaxis );
		b.rotate_csa( ca, sa, uaxis );
		c.rotate_csa( ca, sa, uaxis );
		//a.set(1);
		//b.set(2);
		//c.set(3);
	};

	// ==== generation


	// http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    // http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
    // http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.53.1357&rep=rep1&type=pdf
     //  RAND _ ROTATION   Author: Jim Arvo, 1991
     //  This routine maps three values (x[0], x[1], x[2]) in the range [0,1]
     //  into a 3x3 rotation matrix, M.  Uniformly distributed random variables
     //  x0, x1, and x2 create uniformly distributed random rotation matrices.
     //  To create small uniformly distributed "perturbations", supply
     //  samples in the following ranges
     //      x[0] in [ 0, d ]
     //      x[1] in [ 0, 1 ]
     //      x[2] in [ 0, d ]
     // where 0 < d < 1 controls the size of the perturbation.  Any of the
     // random variables may be stratified (or "jittered") for a slightly more
     // even distribution.
     //=========================================================================
	inline void fromRand( const VEC& vrand  ){
        TYPE theta = vrand.x * M_TWO_PI; // Rotation about the pole (Z).
        TYPE phi   = vrand.y * M_TWO_PI; // For direction of pole deflection.
        TYPE z     = vrand.z * 2.0;      // For magnitude of pole deflection.
        // Compute a vector V used for distributing points over the sphere
        // via the reflection I - V Transpose(V).  This formulation of V
        // will guarantee that if x[1] and x[2] are uniformly distributed,
        // the reflected points will be uniform on the sphere.  Note that V
        // has length sqrt(2) to eliminate the 2 in the Householder matrix.
        TYPE r  = sqrt( z );
        TYPE Vx = sin ( phi ) * r;
        TYPE Vy = cos ( phi ) * r;
        TYPE Vz = sqrt( 2.0 - z );
        // Compute the row vector S = Transpose(V) * R, where R is a simple
        // rotation by theta about the z-axis.  No need to compute Sz since
        // it's just Vz.
        TYPE st = sin( theta );
        TYPE ct = cos( theta );
        TYPE Sx = Vx * ct - Vy * st;
        TYPE Sy = Vx * st + Vy * ct;
        // Construct the rotation matrix  ( V Transpose(V) - I ) R, which
        // is equivalent to V S - R.
        xx = Vx * Sx - ct;   xy = Vx * Sy - st;   xz = Vx * Vz;
        yx = Vy * Sx + st;   yy = Vy * Sy - ct;   yz = Vy * Vz;
        zx = Vz * Sx;        zy = Vz * Sy;        zz = 1.0 - z;   // This equals Vz * Vz - 1.0
	}

};

/*
class Mat3i : public Mat3TYPE< int   , Vec3i, Mat3i >{};
class Mat3f : public Mat3TYPE< float , Vec3f, Mat3f >{};
class MAT : public Mat3TYPE< TYPE, VEC, MAT >{};
*/

using Mat3i = Mat3TYPE< int   >;
using Mat3f = Mat3TYPE< float >;
using Mat3d = Mat3TYPE< double>;

inline void convert( const Mat3f& from, Mat3d& to ){ convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); };
inline void convert( const Mat3d& from, Mat3f& to ){ convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); };

#endif

