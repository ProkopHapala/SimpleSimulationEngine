
#ifndef  fastmath_h
#define  fastmath_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include <stdint.h>
#include <stdint.h>

#include "macroUtils.h"

#include "integerOps.h"

#define GOLDEN_RATIO  1.61803398875
#define DEG2RAD  0.0174533
#define RAD2DEG  57.2958

constexpr double M_TWO_PI = M_PI * 2;



//#define sq(a) a*a
template <class TYPE> inline TYPE sq   (TYPE a){ return a*a; }
template <class TYPE> inline TYPE clip(TYPE x, TYPE xmin, TYPE xmax ){ if( x<xmin ) return xmin; if( x>xmax ) return xmax; return x; }




typedef int (*Func1i)( int  );
typedef int (*Func2i)( int, int );

typedef float  (*Func1f)( float  );
typedef double (*Func1d)( double );

typedef float  (*Func2f)( float,  float  );
typedef double (*Func2d)( double, double );

typedef float  (*Func3f)( float,  float,  float  );
typedef double (*Func3d)( double, double, double );

typedef void   (*Func1d2)( double,                 double&, double& );
typedef void   (*Func2d2)( double, double,         double&, double& );
typedef void   (*Func3d2)( double, double, double, double&, double& );

typedef void   (*Func1d3)( double,                 double&, double&, double& );
typedef void   (*Func2d3)( double, double,         double&, double&, double& );
typedef void   (*Func3d3)( double, double, double, double&, double&, double& );




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


inline double clamp( double x, double xmin, double xmax ){
    if(x<xmin){ return xmin; }else{ if(x>xmax) return xmax; };
    return x;
}

inline double clamp_abs( double x, double xmax ){
    if( x>0 ){ if(x>xmax) return xmax; }else{ if(x<-xmax) return -xmax; };
    return x;
}

//inline int   fastFloor( float f         ){ int i=(int)f; if(f<0)i--; return i; }
//inline float fastFract( float f         ){ return f-fastFloor(f);              }
//inline float fastModf ( float f, int& i ){ i=fastFloor(f); return f-i;         }

inline int    fastFloor( double f         ){ int i=(int)f; if(f<0)i--; return i; }
inline double fastFract( double f         ){ return f-fastFloor(f);              }
inline double fastModf ( double f, int& i ){ i=fastFloor(f); return f-i;         }

#include "gonioApprox.h"

// https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions

inline double erf_4_plus(double x){
    double p = 1 + x*( 0.278393 + x*( 0.230389 + x*(0.000972 + x*0.078108 )));
    p=p*p; p=p*p;
    return 1 - 1/p;
}
inline double erf_4(double x){ if(x>0){ return erf_4_plus(x); }else{ return -erf_4_plus(-x); } }


inline double erf_6_plus(double x){
    double p = 1 + x*( 0.0705230784 + x*( 0.0422820123 + x*( 0.0092705272 + x*( 0.0001520143 + x*( 0.0002765672 + x*0.0000430638 )))));
    p=p*p; p=p*p; p=p*p; p=p*p;
    return 1 - 1/p;
}
inline double erf_6(double x){ if(x>0){ return erf_6_plus(x); }else{ return -erf_6_plus(-x); } }



// === Approx erf(x)/x
//   Problem - we need erf(w*x)/x which is function of two parameters x,w


inline double erfx_e6( double x_, double k, double& dy ){
    // approximation of erf(k*x)/x and its derivative with maximum error ~ 1e-6
    double x =x_*k;
    if( x>4.5 ){ double y=1/x_; dy=-y*y; return y; }
    double xx = x*x;
    double even  =  0.9850156202961753  +xx*(-0.02756061032579559  +xx*(-0.00188409579491924  +xx*(-0.003098629936170076 +xx*(-0.001348858853909826  +xx*(-3.98946569988845e-05 ) ) ) ) );
    double odd   = -0.13893350387140332 +xx*(-0.007664292475021448 +xx*( 0.003046826535877866 +xx*( 0.002879338499080343 +xx*( 0.0003260490382458129 +xx*( 1.97093650414204e-06 ) ) ) ) );
    double deven =                           -0.05512122065159118  +xx*(-0.00753638317967696  +xx*(-0.01859177961702045  +xx*(-0.01079087083127861   +xx*(-0.000398946569988845 ) ) ) )  ;
    double dodd  = -0.1389335038714033  +xx*(-0.02299287742506434  +xx*( 0.01523413267938933  +xx*( 0.0201553694935624   +xx*( 0.002934441344212316  +xx*(2.168030154556244e-05 ) ) ) ) );
    double  t = even    + x*odd;
    double dt = deven*x +  dodd;
    double t2 = t*t;
    double t4 = t2*t2;
    double dt8_dx = 8*dt*t*t2*t4;
    double y      = k/(t4*t4 + x);
    dy            = -y*y*(dt8_dx+1);
    // ToDo : We will need rather (dy/x) for Gauss:Coulomb() => we need to fit it onece more
    return y;
}

inline double erfx_e6( double x_, double k ){
    double x =x_*k;
    if( x>4.5 ){ return 1./x_; }
    double xx = x*x;
    double even =  0.9850156202961753  +xx*(-0.02756061032579559  +xx*(-0.00188409579491924  +xx*(-0.003098629936170076 +xx*(-0.001348858853909826  +xx*(-3.98946569988845e-05 ) ) ) ) );
    double odd  = -0.13893350387140332 +xx*(-0.007664292475021448 +xx*( 0.003046826535877866 +xx*( 0.002879338499080343 +xx*( 0.0003260490382458129 +xx*( 1.97093650414204e-06 ) ) ) ) );
    double  t = even + x*odd;
    t*=t; t*=t; t*=t; // ^8
    return k/( t + x );
}

inline double erfx_e9( double x_, double k ){
    double x =x_*k;
    if( x>4.5 ) return 1./x_;
    double xx = x*x;
    if(x<1./k){
        return 1.1283791662308296 +xx*(-0.3761262972953429 +xx*(0.1128363404233098 +xx*(-0.02685603827999912 +xx*(0.005192885862299865 +xx*(-0.0008053004722300972 +xx*(8.004020068129447e-05 ) ) ) ) ) );
    }
    double even =   0.9903386741213333  +xx*( 0.08180278811069948 +xx*( 0.219787883285348  +xx*( 0.0893543139653664  +xx*( 0.0071698531450102   +xx*( 8.644883946761633e-05 ) ) ) ) );
    double odd  =  -0.17511814497584813 +xx*(-0.2010794452848663  +xx*(-0.1692686167813105 +xx*(-0.03129254573733003 +xx*(-0.001037968593234627 +xx*(-3.164137211658646e-06 ) ) ) ) );
    double t = even + x*odd;
    t*=t; t*=t; t*=t; // ^8
    return k/( t + x );
}

inline double exp_p8( double x ){
    // optimized for decreasing exp(-x)
    if(x<-25) return 0;
    x *= 0.125;
    double xx = x*x;
    //double even = 1.0 +xx*(0.5000000000000000 +xx*(0.04166189077950237  +xx*(0.001321435070258156  ) ) );
    //double odd  = 1.0 +xx*(0.1666664718006032 +xx*(0.008304046626191663 +xx*(0.0001332637951696261 ) ) );
    //double p   = even + x*odd;
    double p = (1+x) +
               xx*( 0.5000000000000000   + 0.1666664718006032   *x +
               xx*( 0.04166189077950237  + 0.008304046626191663 *x +
               xx*( 0.001321435070258156 + 0.0001332637951696261*x ) ) );
    p*=p; p*=p; p*=p;
    return p;
}

inline double finiteExp( double x, double& fr,  double beta, double Rcut ){
    // y  = (1-C*x)^17 * (1-B*x)^2 
    // dy = (1-C x)^16   (1-B*x)    (-2 B - 17 C + 19 B C x)
    //printf( " x %g x/Rcut %g Rcut %g ", x, x/Rcut, Rcut );
    const int k=17;
    const double RN  = Rcut*0.5*k;
    const double cor = 1.15/RN; 
    const double C   = beta/k - cor;
    const double B   = 1/Rcut; 
    //printf( " x %g C %g cor %g beta/k %g ", x, C, cor, beta/k );
    double y1 = 1-x*C; 
    double y  = y1*y1;  // ^2
    y=y*y; // ^4
    y=y*y; // ^8
    y=y*y; // ^16
    double ycut = 1-x*B;
    y*=ycut;
    fr = y* ( -2*B - k*C + (k+2)*B*C*x );
    y*=ycut;
    //y=ycut;
    return y*y1;
}

inline double gauss_p8( double x ){ return exp_p8( -x*x );
    /*
    if(x_>5) return 0;
    double x  = x_*x_*0.125;
    double xx = x *x;
    double p = (1-x) +
               xx*( 0.5000000000000000   + -0.1666664718006032   *x +
               xx*( 0.04166189077950237  + -0.008304046626191663 *x +
               xx*( 0.001321435070258156 + -0.0001332637951696261*x ) ) );
    p*=p; p*=p; p*=p;
    return p;
    */
}

template <typename T>
inline T fastExp(T x, size_t n ){
    T e = 1 + x/(1<<n);
    for(int i=0; i<n; i++) e*=e;
    return e;
}

template <typename T>
inline T fastExp_n4(T x){
    T e = 1 + x*0.0625;
    e*=e; e*=e; e*=e; e*=e;
    return e;
}

template <typename T>
inline T fastExp_n4m(T x){
    T e = 1 + x*0.0625;
    if(e<0)e=0; // smooth landing at zero - cut of divergent part
    e*=e; e*=e; e*=e; e*=e;
    return e;
}

template <typename T>
inline T fastExp_n8(T x){
    T e = 1 + x*0.00390625;
    e*=e; e*=e; e*=e; e*=e;
    e*=e; e*=e; e*=e; e*=e;
    return e;
}

inline double pow3(double x) { return x*x*x; }
inline double pow4(double x) { x*=x;            return x*x;     }
inline double pow5(double x) { double x2=x*x;   return x2*x2*x; }
inline double pow6(double x) { x*=x;            return x*x*x;   }
inline double pow7(double x) { double x3=x*x*x; return x3*x3*x; }
inline double pow8(double x) { x*=x; x*=x;      return x*x;     }

inline double powN(double x, uint8_t n) {
    uint8_t mask=1;
    double xi     = x;
    double result = 1.0;
    while(mask<n){
        if(mask&n){ result*=xi; }
        xi*=xi;
        mask<<=1;
    }
    return result;
}

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


inline float fastInvSqrt( float x ){
	long i;
	float xh, y;
	xh  = x*0.5f;
	y   = x;
	i   = *(long*)&y;               // cast float->int
	i   = 0x5f3759df-(i>>1);
	y   = *(float*)&i;              // cast back
	y   = y * ( 1.5f - (xh*y*y) );   // 1st iteration
	y   = y * ( 1.5f - (xh*y*y) );   // 2nd iteration, this can be removed
	return y;
}



// sqrt(1+dx) taylor(dx=0)
// http://m.wolframalpha.com/input/?i=sqrt%281%2Bx%29+taylor+x%3D0
template <typename T> inline T sqrt1_taylor1( T dx ){ return 1 + dx*  0.5; }
template <typename T> inline T sqrt1_taylor2( T dx ){ return 1 + dx*( 0.5 + dx*  -0.125 ); }
template <typename T> inline T sqrt1_taylor3( T dx ){ return 1 + dx*( 0.5 + dx*( -0.125 + dx*  0.0625 ) ); }
template <typename T> inline T sqrt1_taylor4( T dx ){ return 1 + dx*( 0.5 + dx*( -0.125 + dx*( 0.0625 + dx *  -0.0390625 ) ) ); }
template <typename T> inline T sqrt1_taylor5( T dx ){ return 1 + dx*( 0.5 + dx*( -0.125 + dx*( 0.0625 + dx *( -0.0390625 + dx * 0.02734375 ) ) ) ); }

// 1/sqrt(1+dx) taylor(dx=0)
// http://m.wolframalpha.com/input/?i=1%2Fsqrt%281%2Bx%29+taylor+x%3D0
template <typename T> inline T invSqrt1_taylor1( T dx ){ return 1 + dx*  -0.5; }
template <typename T> inline T invSqrt1_taylor2( T dx ){ return 1 + dx*( -0.5 + dx*  0.375 ); }
template <typename T> inline T invSqrt1_taylor3( T dx ){ return 1 + dx*( -0.5 + dx*( 0.375 + dx*  -0.3125 ) ); }
template <typename T> inline T invSqrt1_taylor4( T dx ){ return 1 + dx*( -0.5 + dx*( 0.375 + dx*( -0.3125 + dx * 0.2734375 ) ) ); }
template <typename T> inline T invSqrt1_taylor5( T dx ){ return 1 + dx*( -0.5 + dx*( 0.375 + dx*( -0.3125 + dx *(0.2734375 + dx * -0.24609375 ) ) ) ); }

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

/*
inline int fastFloor( double x ){
     int ix = static_cast <int>(x);
     if( x < 0 ) ix--;
     return ix;
}
*/


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

inline double fhash_Wang( uint32_t h ){
    return (hash_Wang( h )&(0xffff))/((double)(0xffff));
}

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
    double ia    = -0.5/a;
    if( ia>0 ){
        x1       = ( b - sqrtD )*ia;
        x2       = ( b + sqrtD )*ia;
    }else{
        x1       = ( b + sqrtD )*ia;
        x2       = ( b - sqrtD )*ia;
    }
    //printf( "  a,b,c, %f %f %f  x1,x2 %f %f \n", a,b,c, x1, x2 );
    return true;
}

#endif





