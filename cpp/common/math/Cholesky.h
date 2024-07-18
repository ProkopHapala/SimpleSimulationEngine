
#ifndef  Cholesky_h
#define  Cholesky_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

//#include "fastmath.h"
//#include "VecN.h"

#define N_MAX_NEIGH 32
#define TOLERANCE 1e-16

//#define IJ(i,j) n*i+j

namespace Lingebra{

// =============== Cholesky solver (non-sparse)

void CholeskyDecomp_LDLT(double* A, double* L, double* D, int n) {
    for (int j = 0; j < n; j++) {
        double sum = A[j*n + j];
        for (int k = 0; k < j; k++) { sum -= L[j*n + k] * L[j*n + k] * D[k]; }
        D[j] = sum;
        for (int i = j+1; i < n; i++) {
            sum = A[i*n + j];
            for (int k = 0; k < j; k++) {  sum -= L[i*n + k] * L[j*n + k] * D[k]; }
            L[i*n + j] = sum / D[j];
        }
    }
}

void forward_substitution(double* L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        double sum = b[i];
        for (int j = 0; j < i; j++) {
            sum -= L[i*n + j] * y[j];
        }
        y[i] = sum;
    }
}

void forward_substitution_transposed(double* L, double* b, double* x, int n) {
    for (int i = n-1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i+1; j < n; j++) { sum -= L[j*n + i] * x[j]; }
        x[i] = sum;
    }
}

void solve_LDLT(double* L, double* D, double* b, double* x, int n) {
    double* z = (double*)malloc(n * sizeof(double));
    double* y = (double*)malloc(n * sizeof(double));

    forward_substitution(L, b, z, n);
    for (int i = 0; i < n; i++) { y[i] = z[i] / D[i]; } // Diagonal
    forward_substitution_transposed(L, y, x, n);

    free(z);
    free(y);
}


// =============== Sparse Cholesky solver

int find_or_add_to_neighlist(int* neighs, int i, int n, int max_nodes) {
    int* row = neighs + i * N_MAX_NEIGH;
    for (int j = 0; j < N_MAX_NEIGH; j++) {
        if      (row[j] == n) {  return j; } 
        else if (row[j] == 0) { row[j] = n;return j; }
    }
    return -1;  // list is full, couldn't add n
}

void CholeskyDecomp_LDLT_sparse(double* A, double* L, double* D, int* neighs, int n) {
    for (int j = 0; j < n; j++) {
        double sum1 = 0.0;
        int* neighs_j = neighs + j * N_MAX_NEIGH;
        for (int k = 0; k < N_MAX_NEIGH && neighs_j[k] != 0; k++) {
            if (neighs_j[k] < j) {
                int idx = neighs_j[k];
                sum1 += L[j*n + idx] * L[j*n + idx] * D[idx];
            }
        }
        D[j] = A[j*n + j] - sum1;
        for (int i = j+1; i < n; i++) {
            sum1 = 0.0;
            for (int k = 0; k < N_MAX_NEIGH && neighs_j[k] != 0; k++) {
                if (neighs_j[k] < j) {
                    int idx = neighs_j[k];
                    sum1 += L[i*n + idx] * L[j*n + idx] * D[idx];
                }
            }
            double Lij = (A[i*n + j] - sum1) / D[j];
            if (fabs(Lij) > TOLERANCE) {
                L[i*n + j] = Lij;
                find_or_add_to_neighlist(neighs, j, i, n);
                find_or_add_to_neighlist(neighs, i, j, n);
            }
        }
    }
}

void forward_substitution_sparse(double* L, double* b, int* neighs, double* y, int n) {
    for (int i = 0; i < n; i++) {
        double sum1 = 0.0;
        int* neighs_i = neighs + i * N_MAX_NEIGH;
        for (int k = 0; k < N_MAX_NEIGH && neighs_i[k] != 0; k++) {
            int j = neighs_i[k];
            if (j < i) { sum1 += L[i*n + j] * y[j]; }
        }
        
        y[i] = b[i] - sum1;
    }
}

void forward_substitution_transposed_sparse(double* L, double* b, int* neighs, double* x, int n) {
    for (int i = n-1; i >= 0; i--) {
        double sum1 = 0.0;
        int* neighs_i = neighs + i * N_MAX_NEIGH;
        for (int k = 0; k < N_MAX_NEIGH && neighs_i[k] != 0; k++) {
            int j = neighs_i[k];
            if (j > i) { sum1 += L[j*n + i] * x[j]; }
        }
        x[i] = b[i] - sum1;
    }
}

void solve_LDLT_sparse(double* L, double* D, double* b, int* neighs, double* x, int n) {
    double* z = (double*)malloc(n * sizeof(double));
    double* y = (double*)malloc(n * sizeof(double));
    
    forward_substitution_sparse(L, b, neighs, z, n);
    for (int i = 0; i < n; i++) { y[i] = z[i] / D[i]; } // Diagonal 
    forward_substitution_transposed_sparse(L, y, neighs, x, n);
    
    free(z);
    free(y);
}

};


#endif

