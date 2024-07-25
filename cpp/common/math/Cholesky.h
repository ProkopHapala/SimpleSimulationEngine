#pragma once

//#ifndef  Cholesky_h
//#define  Cholesky_h

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

template<typename T>
inline void CholeskyDecomp_LDLT(T* A, T* L, T* D, int n) {
    for (int j = 0; j < n; j++) {
        T sum = A[j*n + j];
        for (int k = 0; k < j; k++) { sum -= L[j*n + k] * L[j*n + k] * D[k]; }
        D[j] = sum;
        for (int i = j+1; i < n; i++) {
            sum = A[i*n + j];
            for (int k = 0; k < j; k++) {  sum -= L[i*n + k] * L[j*n + k] * D[k]; }
            L[i*n + j] = sum / D[j];
        }
    }
}

template<typename T>
inline void forward_substitution(T* L, T* b, T* y, int n) {
    for (int i = 0; i < n; i++) {
        T sum = b[i];
        for (int j = 0; j < i; j++) {
            sum -= L[i*n + j] * y[j];
        }
        y[i] = sum;
    }
}

template<typename T>
inline void forward_substitution_transposed(T* L, T* b, T* x, int n) {
    for (int i = n-1; i >= 0; i--) {
        T sum = b[i];
        for (int j = i+1; j < n; j++) { sum -= L[j*n + i] * x[j]; }
        x[i] = sum;
    }
}

template<typename T>
inline  void solve_LDLT(T* L, T* D, T* b, T* x, int n) {
    T* z = (T*)malloc(n * sizeof(T));
    T* y = (T*)malloc(n * sizeof(T));

    forward_substitution(L, b, z, n);
    for (int i = 0; i < n; i++) { y[i] = z[i] / D[i]; } // Diagonal
    forward_substitution_transposed(L, y, x, n);

    free(z);
    free(y);
}


// =============== Sparse Cholesky solver

inline int find_or_add_to_neighlist(int* neighs, int i, int n, int max_nodes) {
    int* row = neighs + i * N_MAX_NEIGH;
    for (int j = 0; j < N_MAX_NEIGH; j++) {
        if      (row[j] == n) {  return j; } 
        else if (row[j] == 0) { row[j] = n;return j; }
    }
    return -1;  // list is full, couldn't add n
}

template<typename T>
inline void CholeskyDecomp_LDLT_sparse(T* A, T* L, T* D, int* neighs, int n) {
    for (int j = 0; j < n; j++) {
        T sum = 0.0;
        const int* neighs_j = neighs + j * N_MAX_NEIGH;
        for (int k = 0; k < N_MAX_NEIGH && (neighs_j[k]>=0); k++) {
            if (neighs_j[k] < j) {
                int idx = neighs_j[k];
                sum += L[j*n + idx] * L[j*n + idx] * D[idx];
            }
        }
        D[j] = A[j*n + j] - sum;
        for (int i = j+1; i < n; i++) {
            sum = 0.0;
            for (int k = 0; k < N_MAX_NEIGH && (neighs_j[k]>=0); k++) {
                if (neighs_j[k] < j) {
                    int idx = neighs_j[k];
                    sum += L[i*n + idx] * L[j*n + idx] * D[idx];
                }
            }
            T Lij = (A[i*n + j] - sum) / D[j];
            if (fabs(Lij) > TOLERANCE) {
                L[i*n + j] = Lij;
                find_or_add_to_neighlist(neighs, j, i, n);
                find_or_add_to_neighlist(neighs, i, j, n);
            }
        }
    }
}

template<typename T>
inline void forward_substitution_sparse( int n, int m, const T* L, const T* b, T* x, const int* neighs ) {
    T sum[m];
    for (int i = 0; i < n; i++) {
        for(int s=0; s<m;s++){ sum[s]=0.0; }
        const int* neighs_i = neighs + i * N_MAX_NEIGH;
        for (int k = 0; k < N_MAX_NEIGH && (neighs_i[k]>=0); k++) {
            int j = neighs_i[k];
            if (j < i) { 
                for(int s=0; s<m;s++){ sum[s] += L[i*n+j] * x[j*m+s]; }
            }
        }
        for(int s=0; s<m;s++){  int ii=i*m+s; x[ii] = b[ii] - sum[s]; }
    }
}

template<typename T>
inline void forward_substitution_transposed_sparse( int n, int m, const T* L, const T* b, T* x, const int* neighs ) {
    T sum[m];
    for (int i = n-1; i >= 0; i--) {
        //T sum1 = 0.0;
        for(int s=0; s<m;s++){ sum[s]=0.0; }
        const int* neighs_i = neighs + i * N_MAX_NEIGH;
        for (int k = 0; k < N_MAX_NEIGH && (neighs_i[k]>=0); k++) {
            int j = neighs_i[k];
            if (j > i) { 
                //sum1 += L[j*n + i] * x[j]; 
                for(int s=0; s<m;s++){ sum[s] += L[j*n+i] * x[j*m+s]; }
            }
        }
        //x[i] = b[i] - sum1;
        for(int s=0; s<m;s++){  int ii=i*m+s; x[ii] = b[ii] - sum[s]; }
    }
}

template<typename T>
inline void solve_LDLT_sparse(int n, int m, const T* L, const T* D, T* b, T* x, int* neighs ) {
    T* z = new T[n*m];
    T* y = new T[n*m];
    forward_substitution_sparse( n,m,  L, b, z, neighs );
    for (int i = 0; i < n; i++){ y[i] = z[i] / D[i]; } // Diagonal 
    forward_substitution_transposed_sparse(n,m, L, y, x, neighs );
    delete [] z;
    delete [] y;
}

};


//#endif

