#pragma once

//#ifndef  Cholesky_h
//#define  Cholesky_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <unordered_set>

//#include "fastmath.h"
//#include "VecN.h"

#define TOLERANCE 1e-16

//#define IJ(i,j) n*i+j

namespace Lingebra{

// =============== Cholesky solver (non-sparse)

template<typename T>
inline void CholeskyDecomp_LDLT(T* A, T* L, T* D, int n){
    for(int j=0; j<n; j++){
        T sum = A[j*n+j];
        for(int k=0; k<j; k++){ sum -= L[j*n+k] * L[j*n+k] * D[k]; }
        D[j] = sum;
        for(int i=j+1; i<n; i++){
            sum=A[i*n+j];
            for(int k=0; k<j; k++){ sum -= L[i*n+k] * L[j*n+k] * D[k]; }
            L[i*n+j] = sum / D[j];
        }
    }
}


template<typename T>
inline void forward_substitution_m(T* L, T* b, T* y, int n, int m){
    constexpr double tol = 1e-32;
    T sum[m];
    for(int i=0; i<n; i++){
        printf("forward_substitution_T_m() [i=%i] \n", i);
        for (int s=0; s<m; s++){ sum[s]=0.0; };
        for (int j=0; j<i; j++){
            double l = L[i*n+j];
            if( fabs(l)>tol ){  printf("L[%2i,%2i]=%20.20f\n", i,j,l); };
            for (int s=0; s<m; s++){
                sum[s] += l*y[j*m+s];
            }
        }
        for(int s=0; s<m; s++){
            int ii = i*m+s;
            y[ii]  = b[ii] - sum[s];
            //if(s==0){ printf("fws(m=0)[%3i](x=%20.10f,b=%20.10f,sum=%20.10f)\n",i,x[ii],b[ii],sum[s]); }
            if(s==0){ printf("fws(m=0)[%3i]sum=%20.10f)\n",i,sum[s]); }
        }
        //y[i] = b[i] - sum;
    }
}

template<typename T>
inline void forward_substitution_T_m(T* L, T* b, T* y, int n, int m) {
    //constexpr double tol = 1e-32;
    T sum[m];
    for (int i=n-1; i>=0; i--){
        //printf("forward_substitution_T_m() [i=%i] \n", i);
        for (int s=0; s<m; s++){ sum[s]=0.0; };
        for (int j=i+1; j<n; j++){
            double l = L[j*n+i];
            //if( fabs(l)>tol ){  printf("L[%2i,%2i]=%20.20f\n", i,j,l); };
            for (int s=0; s<m; s++){
                sum[s] += l*y[j*m+s];
            }
        }
        for (int s=0; s<m; s++){
            int ii = i*m+s;
            y[ii]  = b[ii] - sum[s];
            //if(s==0){ printf("fws(m=0)[%3i](x=%20.10f,b=%20.10f,sum=%20.10f)\n",i,x[ii],b[ii],sum[s]); }
            //if(s==0){ printf("fwsT(m=0)[%3i]sum=%20.10f)\n",i,sum[s]); }
        }
        //y[i] = b[i] - sum;
    }
}

template<typename T>
inline void forward_substitution(T* L, T* b, T* y, int n) {
    for (int i=0; i<n; i++) {
        T sum = b[i];
        for (int j=0; j<i; j++) {
            sum -= L[i*n+j] * y[j];
        }
        y[i] = sum;
    }
}

template<typename T>
inline void forward_substitution_transposed(T* L, T* b, T* x, int n) {
    for (int i = n-1; i >= 0; i--) {
        T sum = b[i];
        for (int j=i+1; j<n; j++) { sum -= L[j*n+i] * x[j]; }
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

inline int find_or_add_neigh( int i, int ng, int* neighs, int nmax, int iVoid=-1, bool bPrint=false, bool bErrExit=true) {
    int* ngi = neighs+i*nmax;

    if(i==ng){ printf("ERROR in find_or_add_neigh() i=ng[%i=%i] \n", i, ng ); exit(0); }

    for (int k=0; k<nmax; k++) {
        if      ( ngi[k]== ng    ){   return k; } 
        else if ( ngi[k]== iVoid ){ 
            if(bPrint)printf("insert neigh[%i].insert(%i) at %i \n", i, ng, k );    
            ngi[k]=ng; 
            return k; 
        }
    }
    if(bErrExit){ printf("ERROR in find_or_add_neigh() cannot add ng=%i, neighs[n=%i] is full => exit() \n", ng, nmax); exit(0);  }
    return -1;  // list is full, couldn't add n
}

template<typename T>
inline void CholeskyDecomp_LDLT_sparse(T* A, T* L, T* D, int* neighs, int n, int nNeighMax, bool bErrExit=true ) {
    for (int j=0; j<n; j++){ L[j*n+j]=1; };  // Diagonal
    for (int j=0; j<n; j++){
        //printf(" ---- CholeskyDecomp_LDLT_sparse[%i] \n", j );
        T sum = 0.0;
        const int* neighs_j = neighs + j*nNeighMax;
        for (int k=0; k<nNeighMax; k++){
            int ng = neighs_j[k];
            if(ng<0) break;
            if (ng<j){
                T Lij = L[j*n+ng];
                T val = Lij*Lij * D[ng];
                sum  += val;
                //printf(  "sum[%i,%i]  %20.10f   %20.10f \n", j,k,  val,  sum,  Lij,  D[ng]  );
            }
        }
        D[j] = A[j*n+j] - sum;
        //printf(  "Ch[%i] D,A,sum  %20.10f   %20.10f   %20.10f \n", j,  D[j],  A[j*n+j],  sum  );
        for (int i=j+1; i<n; i++){
            sum = 0.0;
            for (int k=0; k<nNeighMax; k++){
                int ng = neighs_j[k];
                if(ng<0)break;
                if (ng<j){
                    sum += L[i*n+ng] * L[j*n+ng] * D[ng];
                }
            }
            if(i==j){ printf("ERROR in CholeskyDecomp_LDLT_sparse i=j[%i=%i] \n", i, j ); exit(0); }
            T Lij = (A[i*n+j]-sum)/D[j];
            if (fabs(Lij) > TOLERANCE) {
                L[i*n+j] = Lij;
                find_or_add_neigh( j, i, neighs, nNeighMax,-1,false,bErrExit);
                find_or_add_neigh( i, j, neighs, nNeighMax,-1,false,bErrExit);
            }
        }
    }
}

//#include <unordered_set>;

template<typename T>
void CholeskyDecomp_LDLT_sparse_set(T* A, T* L, T* D, int* neighs, int n, int nNeighMax, T tol=1.e-16 ) {
    // Initialize D and L
    for (int j = 0; j < n; j++) {
        D[j] = T(0);
        for (int i = 0; i < n; i++) {
            L[j * n + i] = (i == j) ? T(1) : T(0);
        }
    }
    std::unordered_set<int>  neigh_set[n];
    for (int j = 0; j < n; j++) {
        T sum1 = 0.0;

        // Convert the neighbor list array for column j to an unordered set
        for (int k = 0; k < nNeighMax; ++k) {
            int neighbor = neighs[j * nNeighMax + k];
            if (neighbor > -1 ) {  // Assuming 0 means no neighbor; use another sentinel if needed
                neigh_set[j].insert(neighbor);
            }
        }

        for (const int& k : neigh_set[j] ) {
            if (k < j) {
                sum1 += (L[j * n + k] * L[j * n + k]) * D[k];
            }
        }
        D[j] = A[j * n + j] - sum1;

        for (int i = j + 1; i < n; i++) {
            T sum2 = 0.0;
            for (const int& k : neigh_set[j] ) {
                if (k < j) {
                    sum2 += L[i * n + k] * L[j * n + k] * D[k];
                }
            }
            T Lij = (A[i * n + j] - sum2) / D[j];
            if (std::fabs(Lij) > tol) {
                L[i * n + j] = Lij;
                
                
                neigh_set[j].insert(i);
                neigh_set[i].insert(j);


            }
        }
    }
}

template<typename T>
inline void forward_substitution_sparse( int n, int m, const T* L, const T* b, T* x, const int* neighs, int nNeighMax ) {
    T sum[m];
    for (int i=0; i<n; i++){
        printf("forward_substitution_sparse() [i=%i] \n", i);
        for(int s=0;s<m;s++){ sum[s]=0.0; }
        const int* neighs_i = neighs + i*nNeighMax;
        for (int k=0; k<nNeighMax; k++){
            const int j = neighs_i[k];
            if(j<0)break;
            if (j<i){
                const T l = L[i*n+j];
                //if( fabs(l)>tol ){ 
                printf("L[%2i,%2i]=%20.20f     k=%i\n", i,j,l, k);
                //}; 
                for(int s=0;s<m;s++){ sum[s]+=l*x[j*m+s]; }
            }
        }
        for(int s=0;s<m;s++){ 
            const int ii=i*m+s; x[ii]=b[ii]-sum[s]; 
            //if(s==0){ printf("fws(m=0)[%3i](x=%20.10f,b=%20.10f,sum=%20.10f)\n",i,x[ii],b[ii],sum[s]); }
            //if(s==0){ printf("fws(m=0)[%3i]sum=%20.10f)\n",i,sum[s]); }
        }
    }
}

template<typename T>
inline void forward_substitution_transposed_sparse( int n, int m, const T* L, const T* b, T* x, const int* neighs, int nNeighMax ) {
    T sum[m];
    for (int i=n-1;i>=0;i--){
        for(int s=0;s<m;s++){ sum[s]=0.0; }
        const int* neighs_i = neighs + i*nNeighMax;
        for (int k=0; k<nNeighMax; k++){
            const int j = neighs_i[k];
            if(j<0)break;
            if(j>i){ 
                const T l = L[j*n+i];
                for(int s=0; s<m;s++){ sum[s] += l*x[j*m+s]; }
            }
        }
        for(int s=0;s<m;s++){  const int ii=i*m+s; x[ii]=b[ii]-sum[s]; }
    }
}

template<typename T>
inline void solve_LDLT_sparse(int n, int m, const T* L, const T* D, T* b, T* x, int* neighs, int nNeighMax ) {
    T* z = new T[n*m];
    T* y = new T[n*m];
    forward_substitution_sparse( n,m,  L, b, z, neighs, nNeighMax );
    for (int i = 0; i < n; i++){ y[i] = z[i] / D[i]; } // Diagonal 
    forward_substitution_transposed_sparse(n,m, L, y, x, neighs, nNeighMax );
    delete [] z;
    delete [] y;
}

};


//#endif

