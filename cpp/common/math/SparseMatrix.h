
#ifndef  SparseMatrix_h
#define  SparseMatrix_h

#include "arrayAlgs.h"
#include "CGNE.h"

/// @brief Sparse matrix with constant maximum number of neighbors
///
/// This sparse matrix implementation is designed for efficient storage of:
/// - Nearest neighbor interactions
/// - Trusses
/// - Molecular structures
///
/// @note It is less efficient for storing triangular matrices (e.g., those resulting from sparse Cholesky factorization)
///
/// @details The assumption of constant maximum neighbor count provides several advantages:
/// - Efficient memory layout
/// - Easier parallelization
/// - Faster insertion of new elements
template<typename T>
class SparseMatrix { public:
    int  n;       // dimension (number of rows)
    int  m;       // maximum number of neighbors (i.e. non-zero elements) in each column
    int* nng  =0; // [n]   number of neighbors (non-zero values) in each row
    T*   vals =0; // [m*n] folded non-zero values in the matrix 
    int* inds =0; // [m*n] folded column-index (j) for each non-zero value, -1 for empty slots

    void realloc(int n_,int m_){ 
        m=m_; n=n_;
        _realloc0(nng ,n   ,  0  );
        _realloc0(inds,n*m , -1  );
        _realloc0(vals,n*m , (T)0 );
    }

    __attribute__((pure)) 
    __attribute__((hot)) 
    double get(int i, int j) const {
        // use binary search in sorted array to check if exist such non-zero value (neighbor), and return it, otherwise return 0;         
        const int k = i*m + binarySearch_ignor( j, m, inds+i*m );
        if( inds[k] == j ){ return vals[k]; }else{ return 0; };
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
        //printf( "SparseMatrix::dot_dens_vector()\n" );
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
        //printf( "SparseMatrix::dot_dens_vector_m()\n" );
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
        //if( ni_!=ni ){ bErr=true; if(bPrint){ printf( "checkRow[%i] ni=%i ref.ni=%i | tol=%g\n", i, ni, ni_, tol ); }; return true; }
        if( ni_>ni ){ bErr=true; if(bPrint){ printf( "checkRow[%i] ni=%i ref.ni=%i | tol=%g\n", i, ni, ni_, tol ); }; return true; }
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


void sparse_fsai( const SparseMatrix<double>& A, SparseMatrix<double>& G, double tol=1e-16, int niter=100 ) {
    int n = A.n;
    //std::vector<double> G(n*n, 0.0);  // Initialize G as a dense matrix for simplicity

    for (int i=0; i<n; ++i){ // loop over rows
        // Step 1: Determine the sparsity pattern for row i
        std::vector<int> pattern;
        int ni = A.nng[i];
        int i0 = A.m*i;
        const int* indi = A.inds + i0;
        for (int k = 0; k < ni; ++k) {
            if (indi[k] <= i) {    // only upper triangular (?)
                pattern.push_back(indi[k]);
            }
        }

        //if (pattern.size() > ni_max) {  pattern.resize(ni_max); }   // Ensure we don't exceed p non-zero elements
        //int pattern_size = pattern.size();

        // Step 2: Set up the least squares problem
        std::vector<double> A_i(ni*ni,0.0);
        std::vector<double> b_i(ni, 0.0);
        std::vector<double> g_i(ni, 0.0);

        // Build local linear system Ai*bi=gi  to fit one row of approximate matrix G   
        for (int k=0; k<ni; ++k){
            int row = pattern[k];
            //int row_ni = A.lens[row];
            //int row_i0 = A.i0s[row];
            int row_ni = A.nng[i];
            int row_i0 = A.m*row;
            const int*    row_indi   = A.inds + row_i0;
            const double* row_values = A.vals + row_i0;
            for (int l=0; l<ni; ++l){
                int col = pattern[l];
                for (int m = 0; m < row_ni; ++m){
                    if (row_indi[m] == col){
                        A_i[ k*ni + l ] = row_values[m];
                        break;
                    }
                }
            }
            if(row==i){  b_i[k]=1.0;  }  // basis vectors like b4 = {0,0,0,1,0,0}
        }

        CGNE( ni, A_i.data(), g_i.data(), b_i.data(), niter );

        for (int k=0; k<ni; ++k){
            int j = pattern[k];
            //G[i*n+j] = g_i[k];
            G.set( i,j, g_i[k] );
        }

    }  // for i     // loop over rows
}

#endif

