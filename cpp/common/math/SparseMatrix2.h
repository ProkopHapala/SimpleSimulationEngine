
#ifndef  SparseMatrix2_h
#define  SparseMatrix2_h

#include "arrayAlgs.h"

template<typename T>
class SparseMatrix2 { public:
    int  n;       // dimension (number of rows)
    int* lens =0; // number of non-zero elements in each row
    int* i0s  =0; // index of row begginning in the folded array () 
    T*   vals =0; // folded non-zero values in the matrix 
    int* inds =0; // folded column-index (j) for each non-zero value

    void dot_dens_vector(const T* v, T* out) const {
        for (int i=0; i<n; i++){ // iteratie over rows
            int ni = lens[i];    // number of non-zero values in the row
            int i0 = i0s[i];     // start index of the row in the folded arrays
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
        int ni = lens[i];  // number of non-zero values in the row
        int i0 = i0s[i];   // start index of the row in the folded arrays
        int k= binSearchBetween( j, ni, inds+i0 );
        if( inds[k] == j ){ return vals[i0s+k]; }else{ return 0; };
    }

};


#endif

