
#ifndef  fastmath_h
#define  fastmath_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include <stdint.h>
#include <stdint.h>

#include "integerOps.h"

#define GOLDEN_RATIO  1.61803398875
#define DEG2RAD  0.0174533
#define RAD2DEG  57.2958

constexpr double M_TWO_PI = M_PI * 2;

#define SWAP( a, b, TYPE ) { TYPE t = a; a = b; b = t; }

//#define sq(a) a*a
template <class TYPE> inline TYPE sq   (TYPE a){ return a*a; }
template <class TYPE> inline TYPE clip(TYPE x, TYPE xmin, TYPE xmax ){ if( x<xmin ) return xmin; if( x>xmax ) return xmax; return x; }


#define _max(a,b)      ((a>b)?a:b)
#define _min(a,b)      ((a<b)?a:b)
#define _abs(a)        ((a>0)?a:-a)
//#define _clamp(x,a,b)  max(a, min(b, x))
//#define _clamp(n, lower, upper) if (n < lower) n= lower; else if (n > upper) n= upper
#define _clamp(a, lower, upper) ((a>lower)?((a<upper)?a:upper):lower)

#define _minit( i, x, imin, xmin )  if( x<xmin ){ xmin=x; imin=i; }
#define _maxit( i, x, imax, xmax )  if( x>xmax ){ xmax=x; imax=i; }

#define _circ_inc( i, n )   i++; if(i>=n) i=0;
#define _circ_dec( i, n )   i--; if(i< 0) i=n-1;


typedef float  (*Func1f)( float  );
typedef double (*Func1d)( double );

typedef float  (*Func2f)( float,  float  );
typedef double (*Func2d)( double, double );

typedef float  (*Func3f)( float,  float,  float  );
typedef double (*Func3d)( double, double, double );

inline double x2grid( double x, double xstep, double invXstep, int& ix ){
    double x_=x*invXstep;
    ix=(int)x_;
    //printf( " %f %f %i \n", x, x_, ix );
    return x - ix*xstep;
}

inline double dangle(double da){
    if      (da> M_PI){ return da - 2*M_PI; }
    else if (da<-M_PI){ return da + 2*M_PI; }
    return da;
}

inline double clamp_abs( double x, double xmax ){
    if( x>0 ){ if(x>xmax) return xmax; }else{ if(x<-xmax) return -xmax; };
    return x;
}


#include "gonioApprox.h"

// https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions

inline double erf_4_plus(double x){
    double p = 1 + x*( 0.278393d + x*( 0.230389d + x*(0.000972d + x*0.078108d )));
    p=p*p; p=p*p;
    return 1 - 1/p;
}
inline double erf_4(double x){ if(x>0){ return erf_4_plus(x); }else{ return -erf_4_plus(-x); } }


inline double erf_6_plus(double x){
    double p = 1 + x*( 0.0705230784d + x*( 0.0422820123d + x*( 0.0092705272d + x*( 0.0001520143d + x*( 0.0002765672d + x*0.0000430638d ))))); 
    p=p*p; p=p*p; p=p*p; p=p*p;
    return 1 - 1/p;
}
inline double erf_6(double x){ if(x>0){ return erf_6_plus(x); }else{ return -erf_6_plus(-x); } }


// from http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
inline double fastPow(double a, double b) {
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
  return u.d;
}

// from http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// should be much more precise with large b
inline double fastPrecisePow(double a, double b) {
  // calculate approximation with fraction of the exponent
  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;

  // exponentiation by squaring with the exponent's integer part
  // double r = u.d makes everything much slower, not sure why
  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }

  return r * u.d;
}

/*
template <class FLOAT,class INT> INT fastFloor( FLOAT x ){
    if( x > 0 ){
        INT ix = static_cast <INT>(x);
        //dx = x - ix;
        return ix;
    }else{
        INT ix = static_cast <INT>(-x);
        return 1-ix;
    }
};
*/

inline int fastFloor( double x ){
     int ix = static_cast <int>(x);
     if( x < 0 ) ix--;
     return ix;
}


/*
template <class TYPE>
inline clamp( TYPE x, TYPE xmin, TYPE xmax ){
	if( x<xmin ) return xmin;
	if( x>xmax ) return xmax;
	return x;
}
*/

// ========= random ===========

const  float INV_RAND_MAX = 1.0f/RAND_MAX;
inline float randf(){ return INV_RAND_MAX*rand(); }
inline float randf( float min, float max ){ return randf()*( max - min ) + min; }

// there are some examples of hash functions
// https://en.wikipedia.org/wiki/Linear_congruential_generator
// https://en.wikipedia.org/wiki/Xorshift
// https://gist.github.com/badboy/6267743

// ========= Treshold functions ( Sigmoide, hevyside etc. ) ===========

template <class TYPE>
inline TYPE trashold_step( TYPE x, TYPE x1 ){
	if   (x<x1){ return 0.0; }
	else       { return 1.0; }
}

template <class TYPE>
inline TYPE trashold_lin( TYPE x, TYPE x1, TYPE x2 ){
	if      (x<x1){ return 0.0; }
	else if (x>x2){ return 1.0; }
	else    {       return (x-x1)/(x2-x1); };
}

template <class TYPE>
inline TYPE trashold_cub( TYPE x, TYPE x1, TYPE x2 ){
	if      (x<x1){ return 0.0; }
	else if (x>x2){ return 1.0; }
	else    {  double a =(x-x1)/(x2-x1); return a*a*( 3 - 2*a );  };
}

inline bool quadratic_roots( double a, double b, double c,  double& x1, double& x2 ){
    double D     = b*b - 4*a*c;
    if (D < 0) return false;
    double sqrtD = sqrt( D );
    double ia    = -0.5d/a;
    x1         = ( b - sqrtD )*ia;
    x2         = ( b + sqrtD )*ia;
    //printf( "  a,b,c, %f %f %f  x1,x2 %f %f \n", a,b,c, x1, x2 );
    return true;
}

#endif





