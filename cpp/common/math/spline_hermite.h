
#ifndef  spline_hermite_h
#define  spline_hermite_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <Vec3.h>

//========================
//   Hermite splines
//========================
// https://en.wikipedia.org/wiki/Cubic_Hermite_spline
//
//       x3   x2   x   1
//  ----------------------
//  y0   2   -3        1
//  y1  -2   +3
// dy0   1   -2    1
// dy1   1   -1


//         x3   x2   x   1
//  ----------------------
//  p-1   -1   +2   -1       /2
//  p0    +3   -5       +2   /2
//  p1    -3   +4   +1       /2
//  p2    +1   -1            /2



// ============ optimized

namespace Spline_Hermite{

template<typename T>
struct Sampler{
    T    dx;
    int  ix;
    T y0,dy0,y01,p2,p3;

    inline void seek(T x){
        ix = (int)x;
        dx = x-ix;
    }
    inline void preval(T* buff){
        buff  += ix;
        T ym,y1,yp,dy1;
        ym  = buff[0];
        y0  = buff[1];
        y1  = buff[2];
        yp  = buff[3];
        //printf( "y[{1,2,3,4}] %g %g %g %g \n", ym,y0,y1,yp  );
        dy0 = (y1-ym)*0.5;
        dy1 = (yp-y0)*0.5;
        y01 = y0-y1;
        p2  = (-3*y01 -2*dy0 - dy1);
        p3  = ( 2*y01 +  dy0 + dy1)*dx;
    }
    inline void prepare(T x, T* buff){
        seek  (x);
        preval(buff);
    }
    //inline double y   (){ return y0 - dx*y01; }
    inline double y   (){ return y0 + dx*(dy0 + dx*(  p2 +   p3)); }
    inline double dy  (){ return          dy0 + dx*(2*p2 + 3*p3) ; }
    inline double ddyl(){ return                    2*p2 + 6*p3  ; }
};


const static double C[4][4] = {
{ 1, 0, 0, 0 },
{ 0, 1, 0, 0 },
{ 1, 1, 1, 1 },
{ 0, 1, 2, 3 }
};

const static double B[4][4] = {
{  1,  0,  0,  0 },
{  0,  1,  0,  0 },
{ -3, -2,  3, -1 },
{  2,  1, -2,  1 }
};

template <class T>
inline T val( T x,    T y0, T y1, T dy0, T dy1 ){
	T y01 = y0-y1;
	return      y0
		+x*(           dy0
		+x*( -3*y01 -2*dy0 - dy1
		+x*(  2*y01 +  dy0 + dy1 )));
}

template <class T>
inline T dval( T x,    T y0, T y1, T dy0, T dy1 ){
	T y01 = y0-y1;
	return                 dy0
		+x*( 2*( -3*y01 -2*dy0 - dy1 )
		+x*  3*(  2*y01 +  dy0 + dy1 ));
}

template <class T>
inline T ddval( T x, T y0, T y1, T dy0, T dy1 ){
	T y01 = y0-y1;
	return 2*( -3*y01 -2*dy0 - dy1 )
		+x*6*(  2*y01 +  dy0 + dy1 );
}

template <class T>
inline void valdval( T x, T& val, T& dval, T y0, T y1, T dy0, T dy1 ){
    T y01 = y0-y1;
    T p2  = (-3*y01 -2*dy0 - dy1)*x;
    T p3  = ( 2*y01 +  dy0 + dy1)*x*x;
    val   =  y0 + x*(dy0 +   p2 +   p3);
	dval  =          dy0 + 2*p2 + 3*p3;
}

template <class T>
inline void valdd( T x, T& val, T& dval, T& ddval, T y0, T y1, T dy0, T dy1 ){
    T y01 = y0-y1;
    T p2  = (-3*y01 -2*dy0 - dy1);
    T p3  = ( 2*y01 +  dy0 + dy1)*x;
    val   =  y0 + x*(dy0 + x*(  p2 +   p3));
	dval  =          dy0 + x*(2*p2 + 3*p3);
	ddval =                   2*p2 + 6*p3;
}

template <class T>
inline void basis( T x, T& c0, T& c1, T& d0, T& d1 ){
	T x2   = x*x;
	T K    =  x2*(x - 1);
	c0        =  2*K - x2 + 1;   //    2*x3 - 3*x2 + 1
	c1        = -2*K + x2    ;   //   -2*x3 + 3*x2
	d0        =    K - x2 + x;   //      x3 - 2*x2 + x
	d1        =    K         ;   //      x3 -   x2
}

template <class T>
inline void dbasis( T x, T& c0, T& c1, T& d0, T& d1 ){
	T K    =  3*x*(x - 1);
	c0        =  2*K        ;   //    6*x2 - 6*x
	c1        = -2*K        ;   //   -6*x2 + 6*x
	d0        =    K - x + 1;   //    3*x2 - 4*x + 1
	d1        =    K + x    ;   //    3*x2 - 2*x
}

template <class T>
inline void ddbasis( T x, T& c0, T& c1, T& d0, T& d1 ){
//               x3     x2    x  1
	T x6   =  6*x;
	c0        =  x6 + x6 -  6;   //    12*x - 6
	c1        =   6 - x6 - x6;   //   -12*x + 6
	d0        =  x6 -  4;        //     6*x - 4
	d1        =  x6 -  2;        //     6*x - 2
}

template <class T>
inline void curve_point( T u, const Vec3T<T>& p0, const Vec3T<T>& p1, const Vec3T<T>& t0, const Vec3T<T>& t1,	Vec3T<T>& p ){
	T c0,c1,d0,d1;
	basis<T>( u,  c0, c1, d0, d1 );
	p.set_mul( p0, c0 ); p.add_mul( p1, c1 ); p.add_mul( t0, d0 ); p.add_mul( t1, d1 );
}

template <class T>
inline void curve_tangent( T u, const Vec3T<T>& p0, const Vec3T<T>& p1, const Vec3T<T>& t0, const Vec3T<T>& t1,	Vec3T<T>& t ){
	T c0,c1,d0,d1;
	dbasis<T>( u,  c0, c1, d0, d1 );
	t.set_mul( p0, c0 ); t.add_mul( p1, c1 ); t.add_mul( t0, d0 ); t.add_mul( t1, d1 );
}

template <class T>
inline void curve_accel( T u, const Vec3T<T>& p0, const Vec3T<T>& p1, const Vec3T<T>& t0, const Vec3T<T>& t1,	Vec3T<T>& a ){
	T c0,c1,d0,d1;
	ddbasis<T>( u,  c0, c1, d0, d1 );
	a.set_mul( p0, c0 ); a.add_mul( p1, c1 ); a.add_mul( t0, d0 ); a.add_mul( t1, d1 );
}

template <class T>
inline T val2D( T x, T y,
	T f00, T f01, T f02, T f03,
	T f10, T f11, T f12, T f13,
	T f20, T f21, T f22, T f23,
	T f30, T f31, T f32, T f33
){
	T f0 = val<T>( x, f01, f02, 0.5*(f02-f00), 0.5*(f03-f01) );
	T f1 = val<T>( x, f11, f12, 0.5*(f12-f10), 0.5*(f13-f11) );
	T f2 = val<T>( x, f21, f22, 0.5*(f22-f20), 0.5*(f23-f21) );
	T f3 = val<T>( x, f31, f32, 0.5*(f32-f30), 0.5*(f33-f31) );
	return val<T>( y, f1, f2, 0.5*(f2-f0), 0.5*(f3-f1) );
}

template <class T>
inline T val( T x, const T * fs ){
	T f0 = fs[0];
	T f1 = fs[1];
	T f2 = fs[2];
	T f3 = fs[3];
	return val<T>( x, f1, f2, (f2-f0)*0.5, (f3-f1)*0.5 );
}

template <class T>
inline T dval( T x, const T * fs ){
	T f0 = fs[0];
	T f1 = fs[1];
	T f2 = fs[2];
	T f3 = fs[3];
	return dval<T>( x, f1, f2, (f2-f0)*0.5, (f3-f1)*0.5 );
}

template <class T>
inline T ddval( T x, const T * fs ){
	T f0 = fs[0];
	T f1 = fs[1];
	T f2 = fs[2];
	T f3 = fs[3];
	return ddval<T>( x, f1, f2, (f2-f0)*0.5, (f3-f1)*0.5 );
}

template <class T>
inline T value( T s, const T* ys ){
    //double  s  = r*invdr;
    int    is  = (int)s+1;
    T      dr  = s - is;
    return val<T>( dr, ys+is );
}

template <class T>
inline T deriv( T s, const T* ys ){
    //double  s  = r*invdr;
    int    is  = (int)s+1;
    T      dr  = s - is;
    return dval<T>( dr, ys+is );
}

template <class T>
inline T dderiv( T s, const T* ys ){
    //double  s  = r*invdr;
    int    is  = (int)s+1;
    T      dr  = s - is;
    return ddval<T>( dr, ys+is );
}


/*
template <class T>
inline void valAndDeriv( T s, const T* ys, double& val, double& dval ){
    //double  s  = r*invdr;
    int    is = (int)s;
    T      x  = s - is;

    T ym  =ys[i  ];
    T y0  =ys[i+1];
    T y1  =ys[i+2];
    T yp  =ys[i+3];

    T dy0 = (y1-ym)*0.5;
    T dy1 = (yp-y0)*0.5;
    T y01 = y0-y1;

    T p2  = (-3*y01 -2*dy0 - dy1)*x;
    T p3  = ( 2*y01 +  dy0 + dy1)*x*x;
    val a =  y0 + x*(dy0 +   p2 +   p3);
	dval  =          dy0 + 2*p2 + 3*p3;
}
*/

template <class T>
inline void valderiv( T s, const T* ys, T& val, T& dval ){
    int    is  = (int)s;
    T      x   =  s - is;
    ys+=is;
    T ym = ys[0];
	T y0 = ys[1];
	T y1 = ys[2];
	T yp = ys[3];
    Spline_Hermite::valdval( x, val, dval, y0, y1, (y1-ym)*0.5, (yp-y0)*0.5 ); // Overlap
}

template <class T>
inline void valdd( T s, const T* ys, T& val, T& dval, T& ddval ){
    int    is  = (int)s;
    T      x   =  s - is;
    ys+=is;
    T ym = ys[0];
	T y0 = ys[1];
	T y1 = ys[2];
	T yp = ys[3];
	//              valdd( x, val, dval, ddval, y0, y1, dy0,          dy1 );
    Spline_Hermite::valdd( x, val, dval, ddval, y0, y1, (y1-ym)*0.5, (yp-y0)*0.5 ); // Overlap
}

template <class T>
inline T val2D( T x, T y, const  T * f0s, const  T * f1s, const  T * f2s, const  T * f3s ){
	T f0 = val<T>( x, f0s );
	T f1 = val<T>( x, f1s );
	T f2 = val<T>( x, f2s );
	T f3 = val<T>( x, f3s );
	return val<T>( y, f1, f2, 0.5*(f2-f0), 0.5*(f3-f1) );
}


template <class T>
inline T dval2D( T x, T y, T& dfx, T& dfy,const  T * f0s,const  T * f1s,const  T * f2s,const  T * f3s ){
	T f00=f0s[0]; T f01=f0s[1]; T f02=f0s[2]; T f03=f0s[3];
	T f10=f1s[0]; T f11=f1s[1]; T f12=f1s[2]; T f13=f1s[3];
	T f20=f2s[0]; T f21=f2s[1]; T f22=f2s[2]; T f23=f2s[3];
	T f30=f3s[0]; T f31=f3s[1]; T f32=f3s[2]; T f33=f3s[3];
	T f0,f1,f2,f3;
	// by x
    f0  = val<T>( y, f10, f20, 0.5*(f20-f00), 0.5*(f30-f10) );
	f1  = val<T>( y, f11, f21, 0.5*(f21-f01), 0.5*(f31-f11) );
	f2  = val<T>( y, f12, f22, 0.5*(f22-f02), 0.5*(f32-f12) );
	f3  = val<T>( y, f13, f23, 0.5*(f23-f03), 0.5*(f33-f13) );
	dfx = dval<T>( x, f1, f2, 0.5*(f2-f0), 0.5*(f3-f1) );
	// by y
	f0  = val<T>( x, f01, f02, 0.5*(f02-f00), 0.5*(f03-f01) );
	f1  = val<T>( x, f11, f12, 0.5*(f12-f10), 0.5*(f13-f11) );
	f2  = val<T>( x, f21, f22, 0.5*(f22-f20), 0.5*(f23-f21) );
	f3  = val<T>( x, f31, f32, 0.5*(f32-f30), 0.5*(f33-f31) );
	f20 = 0.5*(f2-f0);
	f31 = 0.5*(f3-f1);
	dfy = dval<T>( y, f1, f2, f20, f31 );
	return val<T>( y, f1, f2, f20, f31 );
}


template <class T>
inline T dval2D( T x, T y, T& dfx, T& dfy,
	T f00, T f01, T f02, T f03,
	T f10, T f11, T f12, T f13,
	T f20, T f21, T f22, T f23,
	T f30, T f31, T f32, T f33
){
	T f0,f1,f2,f3;
	// by x
    f0  = val<T>( y, f10, f20, 0.5*(f20-f00), 0.5*(f30-f10) );
	f1  = val<T>( y, f11, f21, 0.5*(f21-f01), 0.5*(f31-f11) );
	f2  = val<T>( y, f12, f22, 0.5*(f22-f02), 0.5*(f32-f12) );
	f3  = val<T>( y, f13, f23, 0.5*(f23-f03), 0.5*(f33-f13) );
	dfx = dval<T>( x, f1, f2, 0.5*(f2-f0), 0.5*(f3-f1) );
	// by y
	f0  = val<T>( x, f01, f02, 0.5*(f02-f00), 0.5*(f03-f01) );
	f1  = val<T>( x, f11, f12, 0.5*(f12-f10), 0.5*(f13-f11) );
	f2  = val<T>( x, f21, f22, 0.5*(f22-f20), 0.5*(f23-f21) );
	f3  = val<T>( x, f31, f32, 0.5*(f32-f30), 0.5*(f33-f31) );
	f20 = 0.5*(f2-f0);
	f31 = 0.5*(f3-f1);
	dfy = dval<T>( y, f1, f2, f20, f31 );
	return val<T>( y, f1, f2, f20, f31 );
}




/*


Bicubic spline along line: (e.g. for terrain visibility raytracing)

according to wiki https://en.wikipedia.org/wiki/Bicubic_interpolation
F(x,y) = [y^3,y^2,y,1][B][P][B][x^3,x^2,x,1]      # 4x4 matrix product, quadratic form
where [B] is 4x4 matrix of cubic spline basis functions; and [P] 4x4 is matrix of control points

for direction vector s = [sx,sy]
x = x0 + sx*t
y = y0 + sy*t

F(t) = [(y0 + sy*t)^3,(y0 + sy*t)^2,(y0 + sy*t),1 ] [B][P][B]  [(x0 + sx*t)^3,(x0 + sx*t)^2,(x0 + sx*t),1 ]
which is 6 degree polynominal in t
F(t)  = c0 + c1*t + c2*t^2 + c3*t^3 + c4*t^4 + c5*t^5  +  c6*t^6
and its derivative is 5th degree polynominal
dF(t) = d0 + d1*t + d2*t^2 + d3*t^3 + d4*t^4 + d5*t^5
... coefitients

=> more effitient should be  line-search for maximum ( we even know second derivatives )

*/



};


#endif



