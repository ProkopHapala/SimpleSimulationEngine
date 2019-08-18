#include "Lingebra.h"

extern "C"{

    int eig_Jacobi( int n, double* A, double* V, double* es, double tol, int nMaxIter ){
        return Lingebra::eig_Jacobi( n, A, V, es, tol, nMaxIter );
    }

    double eig_Jacobi_init( int n, double* A, double* V, int* mjs, int* ijmax ){
        return Lingebra::eig_Jacobi_init( n, A, V, mjs, ijmax[0], ijmax[1] );
    }

    double eig_Jacobi_step( int n, double* A, double* V, int* mjs, int* ijmax, double vmax ){
        Lingebra::eig_Jacobi_step( n, A, V, mjs, ijmax[0], ijmax[1], vmax );
        return vmax;
    }

}

