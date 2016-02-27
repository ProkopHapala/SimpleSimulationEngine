
#ifndef  VecN_h
#define  VecN_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

namespace VecN{

	inline double dot(	int n,  double* a, double* b ){
		double sum = 0;
		for (int i=0; i<n; i++ ){		sum+= a[i]*b[i];	} 
		return sum;
	}

	inline void set( int n, double  f,            double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f;	      } }
	inline void add( int n, double  f, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f+b[i];    } }
	inline void mul( int n, double  f, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f*b[i];    } }

	inline void set( int n, double* a,            double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i];      } }
	inline void add( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]+b[i]; } }
	inline void sub( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]-b[i]; } }
	inline void mul( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]*b[i]; } }
	inline void div( int n, double* a, double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]/b[i]; } }
	inline void fma( int n, double* a, double* b, double f, double* out ){ for(int i=0; i<n; i++) { out[i]=a[i]+f*b[i]; }  }

	inline void random_vector ( int n, double xmin, double xmax, double * out ){
		double xrange = xmax - xmin;
		for (int i=0; i<n; i++ ){		out[i] = xmin + xrange*randf();	} 
	}

	inline void print_vector( int n, double * a ){
		for (int i=0; i<n; i++ ){	printf( "%f ", a[i] );	} 
		printf( "\n" );
	}

}

#endif
