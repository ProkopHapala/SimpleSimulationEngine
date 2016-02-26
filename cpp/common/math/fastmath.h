
#ifndef  fastmath_h
#define  fastmath_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include <stdint.h>
#include <stdint.h>

#define SWAP( a, b, TYPE ) { TYPE t = a; a = b; b = t; }

#define sq(a) a*a

#define _max(a,b)      (a>b)?a:b
#define _min(a,b)      (a<b)?a:b
#define _abs(a)        (a>0)?a:-a
#define _clamp(x,a,b)  max(a, min(b, x))

#define _circ_inc( i, n )   i++; (i>=n) ? 0   : i; 
#define _circ_dec( i, n )   i--; (i< 0) ? n-1 : i;


#include "gonioApprox.h"

typedef double (*Func1d)( double );

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

inline int rand_hash ( int r ){	return 1664525*r ^ 1013904223; }

inline int rand_hash2( int r ){
	r = 1664525*r ^ 1013904223;
	r = 1664525*r ^ 1013904223;
	return r;
}

inline unsigned int hash_Knuth( unsigned int i ){
	return ( i * 2654435761 >> 16 );
}

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

#endif





