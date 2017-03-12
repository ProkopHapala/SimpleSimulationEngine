
#ifndef  Lingebra_h
#define  Lingebra_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "VecN.h"


//#define IJ(i,j) n*i+j

namespace Lingebra{

	// ==== function declarations

    double ** from_continuous( int m, int n, double *p );
	double ** new_matrix     ( int m, int n );
	double ** delete_matrix  ( int m,        double** A );
	void transpose           ( int m, int n, double** A, double** TA );
	void dot                 ( int m, int n, double** A, double* x, double* out );
	void dotT                ( int m, int n, double** A, double* x, double* out );

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

	int eig_Jacobi( int n, double* A, double* V, double* es, double tol, int nMaxIter  );
	double eig_Jacobi_init( int n, double* A, double* V, int* mjs, int& imax, int& jmax );
	void eig_Jacobi_step  ( int n, double* A, double* V, int* mjs, int& imax, int& jmax, double& vmax );

	// ==== inline function

	inline void set( int m, int n, double   f,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, f,       out[i] );  } }
	inline void add( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, f, B[i], out[i] );  } }
	inline void mul( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, f, B[i], out[i] );  } }

	inline void set( int m, int n, double** A,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, A[i],       out[i] );   } }
	inline void add( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, A[i], B[i], out[i] );   } }
	inline void sub( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::sub( n, A[i], B[i], out[i] );   } }
	inline void mul( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, A[i], B[i], out[i] );   } }
	inline void div( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::div( n, A[i], B[i], out[i] );   } }

};

#endif

