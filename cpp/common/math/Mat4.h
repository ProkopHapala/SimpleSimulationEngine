
#ifndef  Mat4_h
#define  Mat4_h

#include <math.h>
#include <cstdlib>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

//template <class T, class VEC, class MAT>
//template <class T, class VEC>
template <class T>
class Mat4T{
	using VEC = Quat4T<T>;
	using MAT = Mat4T<T>;
	public:
	union{
		struct{
			T xx,xy,xz,xw;
			T yx,yy,yz,yw;
			T zx,zy,zz,zw;
			T wx,wy,wz,ww;
		};
		struct{
			T ax,ay,az,aw;
			T bx,by,bz,bw;
			T cx,cy,cz,cw;
			T dx,dy,dz,dw;
		};
		struct{	VEC a,b,c,d;    };
		struct{	VEC px,py,pz,s; };
		VEC  vecs [4];
		T    array[16];
		T    arr2d[4][4];
	};

// ====== initialization

	inline void setOne(     ){ xx=yy=zz=ww=1; xy=xz=xw=yx=yz=yw=zx=zy=zw=wx=wy=wz=0; };
	inline void set   ( T f ){ xx=yy=zz=ww=f; xy=xz=xw=yx=yz=yw=zx=zy=zw=wx=wy=wz=0; };

	inline void set_outer  ( const VEC& a, const VEC& b ){
		xx=a.x*b.x; xy=a.x*b.y; xz=a.x*b.z; xw=a.x*b.w;
		yx=a.y*b.x; yy=a.y*b.y; yz=a.y*b.z; yw=a.y*b.w;
		zx=a.z*b.x; zy=a.z*b.y; zz=a.z*b.z; zw=a.z*b.w;
		wx=a.w*b.x; wy=a.w*b.y; wz=a.w*b.z; ww=a.w*b.w;
	}

	inline void diag_add( T f ){ xx+=f; yy+=f; zz+=f; ww+=f; };

	inline VEC getColx(){ VEC out; out.x = xx; out.y = yx; out.z = zx; out.w = wx; return out; };
    inline VEC getColy(){ VEC out; out.x = xy; out.y = yy; out.z = zy; out.w = wy; return out; };
    inline VEC getColz(){ VEC out; out.x = xz; out.y = yz; out.z = zz; out.w = wz; return out; };
    inline VEC getColw(){ VEC out; out.x = xz; out.y = yz; out.z = zz; out.w = ww; return out; };

	inline void  colx_to( VEC& out){ out.x = xx; out.y = yx; out.z = zx; out.w = wx; };
    inline void  coly_to( VEC& out){ out.x = xy; out.y = yy; out.z = zy; out.w = wy; };
    inline void  colz_to( VEC& out){ out.x = xz; out.y = yz; out.z = zz; out.w = wz; };
    inline void  colw_to( VEC& out){ out.x = xw; out.y = yw; out.z = zw; out.w = ww; };

	inline void  setColx( const VEC v ){ xx = v.x; yx = v.y; zx = v.z; wx = v.w; };
	inline void  setColy( const VEC v ){ xy = v.x; yy = v.y; zy = v.z; wx = v.w; };
	inline void  setColz( const VEC v ){ xz = v.x; yz = v.y; zz = v.z; wx = v.w; };
	inline void  setColw( const VEC v ){ xw = v.x; yw = v.y; zw = v.z; wx = v.w; };

    inline void setTranspose( MAT& m ){ setColx(m.a); setColy(m.b); setColz(m.c); setColw(m.d); };
    inline MAT   transposed(){ MAT M; M.setColx(a); M.setColy(b); M.setColz(c); M.setColw(d); return M; };

// ====== transpose

	inline MAT operator* ( T f   ) const { MAT m; m.a.set_mul(a,f); m.b.set_mul(b,f); m.c.set_mul(c,f); return m; };

    inline void mul ( T f        ){ a.mul(f);    b.mul(f);    c.mul(f);    d.mul(f); };
    inline void mul ( const VEC& va ){ a.mul(va.x); b.mul(va.y); c.mul(va.z); d.mul(va.w); };

    inline void mulT ( const VEC& va ){
		ax*=va.x; ay*=va.y; az*=va.z; aw*=va.w;
		bx*=va.x; by*=va.y; bz*=va.z; bw*=va.w;
		cx*=va.x; cy*=va.y; cz*=va.z; cw*=va.w;
		dx*=va.x; dy*=va.y; dz*=va.z; dw*=va.w;
	};

// ====== dot product with vector

	inline void dot_to( const VEC&  v, VEC&  vout ) const {
		vout.x = xx*v.x + xy*v.y + xz*v.z + xw*v.w;
		vout.y = yx*v.x + yy*v.y + yz*v.z + yw*v.w;
		vout.z = zx*v.x + zy*v.y + zz*v.z + zw*v.w;
		vout.w = wx*v.x + wy*v.y + wz*v.z + ww*v.w;
	};

	inline void dot_to_T( const VEC&  v, VEC&  vout ) const {
		vout.x = xx*v.x + yx*v.y + zx*v.z + wx*v.w;
		vout.y = xy*v.x + yy*v.y + zy*v.z + wy*v.w;
		vout.z = xz*v.x + yz*v.y + zz*v.z + wz*v.w;
		vout.w = xw*v.x + yw*v.y + zw*v.z + ww*v.w;
	};

	inline VEC dot ( const VEC&  v ) const { VEC vout; dot_to  ( v, vout ); return vout; }
	inline VEC dotT( const VEC&  v ) const { VEC vout; dot_to_T( v, vout ); return vout; }

// ====== matrix multiplication

	inline void set_mmul( const MAT& A, const MAT& B ){ for(int i=0; i<4; i++){ for(int j=0; j<4; j++){
            arr2d[i][j] = A.arr2d[i][0]*B.arr2d[0][j] + A.arr2d[i][1]*B.arr2d[1][j] + A.arr2d[i][2]*B.arr2d[2][j] + A.arr2d[i][3]*B.arr2d[3][j];
        }}
	};

    inline void set_mmul_NT( const MAT& A, const MAT& B ){
        for(int i=0; i<4; i++){ for(int j=0; j<4; j++){
            arr2d[i][j] = A.arr2d[i][0]*B.arr2d[j][0] + A.arr2d[i][1]*B.arr2d[j][1] + A.arr2d[i][2]*B.arr2d[j][2] + A.arr2d[i][3]*B.arr2d[j][3];
        }}
	};

    inline void set_mmul_TN( const MAT& A, const MAT& B ){
        for(int i=0; i<4; i++){ for(int j=0; j<4; j++){
            arr2d[i][j] = A.arr2d[0][i]*B.arr2d[0][j] + A.arr2d[1][i]*B.arr2d[1][j] + A.arr2d[2][i]*B.arr2d[2][j] + A.arr2d[3][i]*B.arr2d[3][j];
        }}
	};

    inline void set_mmul_TT( const MAT& A, const MAT& B ){
        for(int i=0; i<4; i++){ for(int j=0; j<4; j++){
            arr2d[i][j] = A.arr2d[0][i]*B.arr2d[j][0] + A.arr2d[1][i]*B.arr2d[j][1] + A.arr2d[2][i]*B.arr2d[j][2] + A.arr2d[3][i]*B.arr2d[j][3];
        }}
	};

	void mmulL( const MAT& A ){ MAT& M=*this; set_mmul( M, A); }
	void mmulR( const MAT& A ){ MAT& M=*this; set_mmul( A, M); }

    void mmulLT( const MAT& A ){ MAT& M=*this; set_mmul_NT( M, A); }
	void mmulRT( const MAT& A ){ MAT& M=*this; set_mmul_TN( A, M); }

	//   http://www.songho.ca/opengl/gl_projectionmatrix.html
    void setPerspective( T xmin, T xmax, T ymin, T ymax, T zmin, T zmax ){
        //T invdx = xmax-xmin; T invdy = ymax-ymin; T invdz = zmax-zmin; // WARRNING : THIS IS WRONG
        T invdx = 1/(xmax-xmin); T invdy = 1/(ymax-ymin); T invdz = 1/(zmax-zmin);
        array[0 ]  = 2*zmin*invdx; array[1 ] = 0;            array[2 ] =  (xmax+xmin)*invdx;  array[3 ] = 0;
        array[4 ]  = 0;            array[5 ] = 2*zmin*invdy; array[6 ] =  (ymax+ymin)*invdy;  array[7 ] = 0;
        array[8 ]  = 0;            array[9 ] = 0;            array[10] = -(zmax+zmin)*invdz;  array[11] = -2*zmax*zmin*invdz;
        array[12]  = 0;            array[13] = 0;            array[14] = -1;                  array[15] = 0;
        /*
        array[0 ]  = 2*zmin*invdx; array[1 ] = 0;            array[2 ] =  (xmax+xmin)*invdx;  array[3 ] = 0;
        array[4 ]  = 0;            array[5 ] = 2*zmin*invdy; array[6 ] =  (ymax+ymin)*invdy;  array[7 ] = 0;
        array[8 ]  = 0;            array[9 ] = 0;            array[10] = -(zmax+zmin)*invdz;  array[11] = zmax*zmin*invdz;
        array[12]  = 0;            array[13] = 0;            array[14] = -1;                  array[15] = 0;
        */
    }

    //   https://www.khronos.org/registry/OpenGL-Refpages/gl2.1/xhtml/gluPerspective.xml
    void setPerspective( T ctgx, T ctgy, T zmin, T zmax ){
        T invdz = 1/(zmin-zmax);
        array[0 ]  = 2*ctgx; array[1 ] = 0;      array[2 ] =  0;                  array[3 ] = 0;
        array[4 ]  = 0;      array[5 ] = 2*ctgy; array[6 ] =  0;                  array[7 ] = 0;
        array[8 ]  = 0;      array[9 ] = 0;      array[10] = -(zmax+zmin)*invdz;  array[11] = -2*zmax*zmin*invdz;
        array[12]  = 0;      array[13] = 0;      array[14] = -1;                  array[15] = 0;
        // WARNING - Should be probably transposed !!!!!!
    }

    // http://www.songho.ca/opengl/gl_projectionmatrix.html
    void setOrthographic( T W, T H, T zmin, T zmax ){
        T invdz = 1/(zmax-zmin);
        /*
        array[0 ]  = 1/W; array[1 ] = 0;   array[2 ] =  0;       array[3 ] = 0;
        array[4 ]  = 0;   array[5 ] = 1/H; array[6 ] =  0;       array[7 ] = 0;
        array[8 ]  = 0;   array[9 ] = 0;   array[10] = 2*invdz;  array[11] = (zmax+zmin)*invdz;
        array[12]  = 0;   array[13] = 0;   array[14] =  0;       array[15] = 1;
        */
        array[0 ]  = 1/W; array[1 ] = 0;   array[2 ] =  0;                  array[3 ] = 0;
        array[4 ]  = 0;   array[5 ] = 1/H; array[6 ] =  0;                  array[7 ] = 0;
        array[8 ]  = 0;   array[9 ] = 0;   array[10] = 2*invdz;             array[11] = 0;
        array[12]  = 0;   array[13] = 0;   array[14] =  (zmax+zmin)*invdz;  array[15] = 1;
    }

    void setRot( Mat3T<T> M ){
        xx=M.xx; xy=M.xy; xz=M.xz;
		yx=M.yx; yy=M.yy; yz=M.yz;
		zx=M.zx; zy=M.zy; zz=M.zz;
    }

    inline void add_mul(const MAT& m,const VEC& v){
        for(int i=0; i<4; i++){ vecs[i].add_mul( m.vecs[i], v.array[i] ); }
    }

    inline void add_mul(const MAT& m,T f){
        for(int i=0; i<4; i++){ vecs[i].add_mul( m.vecs[i], f ); }
    }

    inline void makeOrthoU(const MAT& m ){
        for(int i=0; i<4; i++){ vecs[i].makeOrthoU( m.vecs[i] ); }
    }

    void normalize(){ for(int i=0;i<4;i++)vecs[i].normalize(); };

    inline T getOrthoForces(MAT& f){
        T c2sum=0;
        for(int i=1;i<4;i++){
            for(int j=0;j<i;j++){
                T c = -vecs[i].dot(vecs[j]);
                f.vecs[i].add_mul( vecs[j], c );
                f.vecs[j].add_mul( vecs[i], c );
                c2sum+=c*c;
                //printf( "getOrthoForces[%i,%i] c %g \n", i, j, c );
            }
        }
        return c2sum;
    }

    inline T orthoStep(const T sc){
        T err2 = 0;
        MAT f; f.set((T)0);
        for(int i=0;i<4;i++){ err2+=sq( vecs[i].normalize_taylor3() ); }
        err2+=getOrthoForces(f);
        for(int i=0;i<4;i++){ vecs[i].add_mul( f.vecs[i], sc ); }
        return err2;
    }

    inline int orthoRun(const T err2conv,const int nMaxIter ){
        T err2;
        MAT f;
        int itr;
        for(itr=0;itr<nMaxIter;itr++){
            err2=0;
            f.set((T)0);
            //for(int i=0;i<4;i++){ err2+=vecs[i].normalize_taylor3(); }
            for(int i=0;i<4;i++){ err2+=sq( vecs[i].normalize_hybrid( 0.04d+err2conv ) ); }
            err2+=getOrthoForces(f);
            //for(int i=0;i<4;i++){ vecs[i].add_mul( f.vecs[i], 0.5d ); }
            add_mul( f, 0.5d );
            //printf( "orthoRun[%i] %g\n", itr, err2 );
            if(err2<err2conv) break;
        }
        return itr;
    }

    //void setPos( Vec3T<T> p ){ xw=p.x; yw=p.y; xw=p.z; }
    void setPos( Vec3T<T> p ){ wx=p.x; wy=p.y; wz=p.z; }

// ====== matrix solver

	void print(){
        printf( " %f %f %f %f \n", ax, ay, az, aw );
        printf( " %f %f %f %f \n", bx, by, bz, bw );
        printf( " %f %f %f %f \n", cx, cy, cz, cw );
        printf( " %f %f %f %f \n", dx, dy, dz, dw );
    }

};

using Mat4i = Mat4T< int   >;
using Mat4f = Mat4T< float >;
using Mat4d = Mat4T< double>;

//inline void convert( const Mat3f& from, Mat3d& to ){ convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); };
//inline void convert( const Mat3d& from, Mat3f& to ){ convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); };

#endif

