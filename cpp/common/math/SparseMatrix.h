
#ifndef  SparseMatrix_h
#define  SparseMatrix_h

#include "arrayAlgs.h"

//__attribute__((pure))

__attribute__((hot)) 
int binarySearch_ignor( int j, int end, const int* inds, int ignore=-1) {
    int start=0;
    while (start < end) {
        int mid = start + (end - start) / 2;
        if ( inds[mid] == ignore || inds[mid] > j) {
            end = mid;
        } else if (inds[mid] < j) {
            start = mid + 1;
        } else {
            return mid;
        }
    }
    return start;  // Return the position where j should be inserted
}


template<typename T>
void insertElement( T x, int i, int n, T* arr ) {
    for(int k=n-1; k>i; --k) { arr[k]=arr[k-1]; }
    arr[i] = x;
}

template<typename T>
class SparseMatrix { public:
    int  n;       // dimension (number of rows)
    int  m;       // maximum number of neighbors (i.e. non-zero elements) in each column
    int* nng  =0; // [n]   number of neighbors (non-zero values) in each row
    T*   vals =0; // [m*n] folded non-zero values in the matrix 
    int* inds =0; // [m*n] folded column-index (j) for each non-zero value, -1 for empty slots

    __attribute__((pure)) 
    __attribute__((hot)) 
    double get(int i, int j) const {
        // use binary search in sorted array to check if exist such non-zero value (neighbor), and return it, otherwise return 0;         
        const int k = i*m + binarySearch_ignor( j, m, inds+i*m );
        if( inds[k] == j ){ return vals[k]; }else{ return 0; };
    }

    void realloc(int n_,int m_){ 
        m=m_; n=n_;
        _realloc0(nng,n,0);
        _realloc0(inds,n*m,-1);
        _realloc0(vals,n*m,0.0);
    }

    __attribute__((hot)) 
    bool set(int i, int j, T x ){
        // use binary search in sorted array to check if exist such non-zero value (neighbor)
        // if it exist, change its value
        // if it does not exist, add new element. But to keep the neighbor list sorted must shift all the element right (like in insertion sort), if number of neighbors reached limit [m] return error
        const int k  = binarySearch_ignor( j, m, inds+i*m );
        const int ik = i*m+k;
        if( inds[ik]==j ){ 
            vals[ik]=x; 
            return false;
        }else{
            if( nng[i]>=m )[[unlikely]]{  printf("ERROR in SparseMatrix::set() cannot add new element row[i=%i].nneigh=%i reached the limit m=%i \n", i, nng[1], m); exit(0); }
            insertElement( j, k, m, inds+i*m );
            insertElement( x, k, m, vals+i*m );
            nng[i]++;
            return true; 
        };
    }

    __attribute__((hot)) 
    void dot_dens_vector(const T* v, T* out) const {
        // printf( "SparseMatrix::dot_dens_vector()\n" );
        for (int i=0; i<n; i++){ // iteratie over rows
            const int ni = nng[i];     // number of non-zero values in the row
            const int i0 = i*m;        // start index of the row in the folded arrays
            const int* indi = inds + i0;
            const T*   Ai   = vals + i0;
            T sum = 0;
            for (int k=0; k<ni; k++) {
                const int j  = indi[k];
                sum   += Ai[k] * v[j];
            }
            out[i] = sum;
        }
    }

    __attribute__((hot)) 
    void dot_dens_vector_m(int ns, const T* v, T* out) const {
        printf( "SparseMatrix::dot_dens_vector_m()\n" );
        for (int i=0; i<n; i++){ // iteratie over rows
            const int ni = nng[i];     // number of non-zero values in the row
            const int i0 = i*m;        // start index of the row in the folded arrays
            T sum[ns]; 
            for (int s=0; s<ns; s++) { sum[s]=0; };
            for (int k=0; k<ni; k++) {
                const int ik = i0+k; 
                const int j0 = inds[ik]*ns;
                const T  aj  = vals[ik];
                for (int s=0; s<ns; s++) {   
                    sum[s] += v[j0+s]*aj; 
                }
            }
            for (int s=0; s<ns; s++) { out[i*ns+s] = sum[s]; }
        }
    }

    bool checkRow( int i, const T* ref, T tol=1e-16, bool bPrint=true ) const {
        bool bErr=false;
        //for (int i=0; i<n; i++){ // iteratie over rows
        const int ni = nng[i];     // number of non-zero values in the row
        const int i0 = i*m;        // start index of the row in the folded arrays
        int ni_ = 0; for (int j=0; j<n; j++){ if( fabs(ref[j])>tol ){ ni_++; }; } // count number of non-zero elements
        if( ni_!=ni ){ bErr=true; if(bPrint){ printf( "checkRow[%i] ni=%i ref.ni=%i | tol=%g\n", i, ni, ni_, tol ); }; return true; }
        for (int k=0; k<ni; k++) {
            const int ik = i0+k; 
            const int j = inds[ik];
            if((j<0)||(j>=n)){ printf( "ERROR in SparseMatrix::checkRow()[i=%i,k=%i] j(%i) is out of range 0..n(%i) ni=%i => exit\n", i,k,j,n, ni ); exit(0); }
            //printf( "checkRow()[%i,%i] j=%i\n", i, k, j );
            const T x    = vals[ik];
            const T refj = ref [j];
            if( fabs( x-refj )>tol ){ bErr=true; if(bPrint){ printf( "checkRow[%i] vals[%i](%g) != ref[%i](%g) | tol=%g\n", i, k, x, j, refj, tol ); }; };
        }
        //}
        return bErr;
    }
    bool checkDens( const T* ref, T tol=1e-16, bool bPrint=true ) const {
        bool bErr=false;
        for (int i=0; i<n; i++){ 
            //printf( "#-- checkDens[%i]\n", i );
            bool b = checkRow( i, ref+i*n, tol, bPrint ); 
            if(b){ printf( "#-- checkDens[%i] does not match! \n", i ); }
            bErr |= b;

        };
        return bErr;
    }

};


#endif

