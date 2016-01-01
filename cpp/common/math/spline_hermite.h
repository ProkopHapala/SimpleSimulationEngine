
#ifndef  spline_hermite_h
#define  spline_hermite_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

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


// ============ optimized

namespace Spline_Hermite{

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



template <class TYPE>  
inline TYPE val( TYPE x,    TYPE y0, TYPE y1, TYPE dy0, TYPE dy1 ){  
	TYPE y01 = y0-y1;
	return      y0
		+x*(           dy0              
		+x*( -3*y01 -2*dy0 - dy1   
		+x*(  2*y01 +  dy0 + dy1 )));
}

template <class TYPE> 
inline TYPE dval( TYPE x,    TYPE y0, TYPE y1, TYPE dy0, TYPE dy1 ){  
	TYPE y01 = y0-y1;
	return                 dy0              
		+x*( 2*( -3*y01 -2*dy0 - dy1 )   
		+x*  3*(  2*y01 +  dy0 + dy1 ));
}

template <class TYPE> 
inline TYPE ddval( TYPE x, TYPE y0, TYPE y1, TYPE dy0, TYPE dy1 ){  
	TYPE y01 = y0-y1;
	return 2*( -3*y01 -2*dy0 - dy1 )   
		+x*6*(  2*y01 +  dy0 + dy1 );
}

template <class TYPE>  
inline void basis( TYPE x, TYPE& c0, TYPE& c1, TYPE& d0, TYPE& d1 ){
	TYPE x2   = x*x;
	TYPE K    =  x2*(x - 1); 
	c0        =  2*K - x2 + 1;   //    2*x3 - 3*x2 + 1
	c1        = -2*K + x2    ;   //   -2*x3 + 3*x2   
	d0        =    K - x2 + x;   //      x3 - 2*x2 + x 
	d1        =    K         ;   //      x3 -   x2
}

template <class TYPE>  
inline void dbasis( TYPE x, TYPE& c0, TYPE& c1, TYPE& d0, TYPE& d1 ){
	TYPE K    =  3*x*(x - 1); 
	c0        =  2*K        ;   //    6*x2 - 6*x
	c1        = -2*K        ;   //   -6*x2 + 6*x   
	d0        =    K - x + 1;   //    3*x2 - 4*x + 1 
	d1        =    K + x    ;   //    3*x2 - 2*x 
}

template <class TYPE>  
inline void ddbasis( TYPE x, TYPE& c0, TYPE& c1, TYPE& d0, TYPE& d1 ){
//               x3     x2    x  1
	TYPE x6   =  6*x; 
	c0        =  x6 + x6 -  6;   //    12*x - 6
	c1        =   6 - x6 - x6;   //   -12*x + 6   
	d0        =  x6 -  4;        //     6*x - 4 
	d1        =  x6 -  2;        //     6*x - 2   
}

template <class TYPE>  
inline void curve_point( TYPE u, const Vec3TYPE<TYPE>& p0, const Vec3TYPE<TYPE>& p1, const Vec3TYPE<TYPE>& t0, const Vec3TYPE<TYPE>& t1,	Vec3TYPE<TYPE>& p ){
	TYPE c0,c1,d0,d1;
	basis<TYPE>( u,  c0, c1, d0, d1 );
	p.set_mul( p0, c0 ); p.add_mul( p1, c1 ); p.add_mul( t0, d0 ); p.add_mul( t1, d1 );
}

template <class TYPE>  
inline void curve_tangent( TYPE u, const Vec3TYPE<TYPE>& p0, const Vec3TYPE<TYPE>& p1, const Vec3TYPE<TYPE>& t0, const Vec3TYPE<TYPE>& t1,	Vec3TYPE<TYPE>& t ){
	TYPE c0,c1,d0,d1;
	dbasis<TYPE>( u,  c0, c1, d0, d1 );
	t.set_mul( p0, c0 ); t.add_mul( p1, c1 ); t.add_mul( t0, d0 ); t.add_mul( t1, d1 );
}

template <class TYPE>  
inline void curve_accel( TYPE u, const Vec3TYPE<TYPE>& p0, const Vec3TYPE<TYPE>& p1, const Vec3TYPE<TYPE>& t0, const Vec3TYPE<TYPE>& t1,	Vec3TYPE<TYPE>& a ){
	TYPE c0,c1,d0,d1;
	ddbasis<TYPE>( u,  c0, c1, d0, d1 );
	a.set_mul( p0, c0 ); a.add_mul( p1, c1 ); a.add_mul( t0, d0 ); a.add_mul( t1, d1 );
}

};


#endif



