
#ifndef  SparseMatrix2_h
#define  SparseMatrix2_h

#include "arrayAlgs.h"

/*

This sparse matrix has varying number of neighbors
That makes it efficient to store e.g. e.g. triangular matrixes (such as resulting from sparse Cholesky factorization)

The storage is densly packed, which makes it problematic for inserting new elements, => we need to re-create the matrix each time valency changes

*/

template<typename T>
class SparseMatrix2 { public:
    int  n;       // dimension (number of rows)
    int* nngs =0; // number of non-zero elements in each row
    int* i0s  =0; // index of row begginning in the folded array () 
    T*   vals =0; // folded non-zero values in the matrix 
    int* inds =0; // folded column-index (j) for each non-zero value

    void dot_dens_vector(const T* v, T* out) const {
        for (int i=0; i<n; i++){ // iteratie over rows
            int ni = nngs[i];    // number of non-zero values in the row
            int i0 = i0s [i];     // start index of the row in the folded arrays
            const int* indi  = inds + i0; // 
            const T*   Ai = vals + i0;
            T sum = 0;
            for (int k = 0; k < ni; k++) {
                int j = indi[k];
                sum += Ai[k] * v[j];
            }
            out[i] = sum;
        }
    }

    double get(int i, int j){
        int ni = nngs[i];  // number of non-zero values in the row
        int i0 = i0s[i];   // start index of the row in the folded arrays
        int k= binSearchBetween( j, ni, inds+i0 );
        if( inds[k] == j ){ return vals[i0s+k]; }else{ return 0; };
    }

    int fromDense( int n_, T* A, T tol ){
        n=n_;
        _realloc0(nngs, n, 0 );
        _realloc0(i0s,  n, -1 );
        int ntot = 0;
        for(int i=0; i<n; i++){
            int ni = countNonZero( n, A+i*n, tol );
            printf( "fromDense()[%i] ni=%i \n", i, ni );
            nngs[i]=ni;
            i0s [i]=ntot;
            ntot+=ni;
        }
        _realloc0(vals, ntot, 0.0 );
        _realloc0(inds, ntot, -1  );
        for(int i=0; i<n; i++){
            int i0 = i0s[i];
            if(i0<0){ printf("WTF? i=%i i0=%i \n", i, i0); exit(0); }
            int ni = exportNonZero( n, A+i*n, tol, vals+i0, inds+i0 );
            if(ni!=nngs[i]){ printf("ERROR in SparseMatrix::fromDense()[%i/%i] exportNonZero().count(%i) != countNonZero().count(%i) => exit \n", i,n, ni, nngs[i] );  exit(0); }
        }
        return ntot;
    }

    // int inds_to_file( FILE* fout ){
    //     return ntot;
    // }

};


#endif

