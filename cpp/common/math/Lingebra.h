
#ifndef  Lingebra_h
#define  Lingebra_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "VecN.h"

namespace Lingebra{

	// ==== function declarations

	void mmul_ik_kj( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ik_jk( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ki_kj( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ki_jk( int ni, int nj, int nk, double** A, double** B, double** out );

	void random_matrix( int m, int n, double xmin, double xmax, double** out );
	void print_matrix( int m, int n, double ** A );

	void   makeQuadricFormMatrix( int m, int n, double * ks, double ** A, double ** Q );
	double evalQudraticForm( int n, double* x, double** Q );
	double evalQudraticFormDirs( int m, int n, double* x, double* k, double** A );

	void GaussElimination    ( int n, double ** A, double * c, int * index );
	void linSolve_gauss      ( int n, double ** A, double * b, int * index, double * x );
	void linSolve_CG         ( int n, double ** A, double * b, double * x );
	void linSolve_BCG        ( int n, int m, double ** A, double * b, double * x );
	void leastSquareFit_Gauss( int n, int m, double ** A, double * b, double * x );

	// ==== inline function 

	inline void set( int m, int n, double   f,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, f,       out[i] );  } }
	inline void add( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, f, B[i], out[i] );  } }
	inline void mul( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, f, B[i], out[i] );  } }

	inline void set( int m, int n, double** A,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, A[i],       out[i] );   } }
	inline void add( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, A[i], B[i], out[i] );   } }
	inline void sub( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::sub( n, A[i], B[i], out[i] );   } }
	inline void mul( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, A[i], B[i], out[i] );   } }
	inline void div( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::div( n, A[i], B[i], out[i] );   } }

	// creates double** from any continuous memory block folowing *p
	double ** from_continuous( int m, int n, double *p ){
		double ** A = new double*[m];
		for (int i=0; i<m; i++ ){ A[i] = &(p[i*n]); }
		return A;
	}

	double ** new_matrix( int m, int n ){
		double ** A = new double*[m];
		for (int i=0; i<m; i++ ){ A[i] = new double[n]; }
		return A;
	}

	double ** delete_matrix( int m, double** A ){
		for (int i=0; i<m; i++ ){ delete A[i]; }
		delete [] A;
	}

	void transpose( int m, int n, double** A, double** TA ){
		for (int i=0; i<m; i++ ){
			for (int j=0; j<n; j++ ){
				TA[i][j] = A[j][i];
			} 
		} 
	}

	void dot( int m, int n, double** A, double* x, double* out ){
		for (int i=0; i<m; i++ ){
			out[i] = VecN::dot( n, A[i], x );
		} 
	}

	void dotT( int m, int n, double** A, double* x, double* out ){
		for (int i=0; i<m; i++ ){
			double doti = 0;
			for (int j=0; j<n; j++ ){ doti += A[j][i] * x[j];	}
			out[i] = doti;
		} 
	}

};

#endif

