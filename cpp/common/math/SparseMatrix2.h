
#ifndef  SparseMatrix2_h
#define  SparseMatrix2_h

#include "arrayAlgs.h"

/// @brief Sparse matrix with varying number of neighbors per row
/// @details This implementation is efficient for storing triangular matrices, such as those resulting from sparse Cholesky factorization. 
/// The storage is densely packed, which optimizes memory usage but makes inserting new elements problematic. 
/// As a result, the matrix needs to be recreated each time the number of non-zero elements (valency) changes (e.g. when changing the number of neighbors in a truss or other interaction topology).

template<typename T>
class SparseMatrix2 { public:
    int  n,ntot;  // dimension (number of rows)
    int* nngs =0; // number of non-zero elements in each row
    int* i0s  =0; // index of row begginning in the folded array () 
    T*   vals =0; // folded non-zero values in the matrix 
    int* inds =0; // folded column-index (j) for each non-zero value

    __attribute__((hot)) 
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

    __attribute__((hot)) 
    T get(int i, int j){
        int ni = nngs[i];  // number of non-zero values in the row
        int i0 = i0s[i];   // start index of the row in the folded arrays
        int k= binSearchBetween( j, ni, inds+i0 );
        if( inds[k] == j ){ return vals[i0s+k]; }else{ return 0; };
    }




    __attribute__((hot)) 
    void fwd_subs_mi( int i, int ns, const T* b, T* x ){
        //printf("SparseMatrix2::fwd_subs_m() [i=%i] \n", i);
        T sum[ns];
        for(int s=0;s<ns;s++){ sum[s]=0.0; }
        const int i0 = i0s [i];
        const int ni = nngs[i];
        for (int k=0; k<ni; k++){
            const int ik = i0+k;
            const int j  = inds[ik]*ns;
            const T   l  = vals[ik];  
            for(int s=0;s<ns;s++){ sum[s]+=l*x[j+s]; }
        }
        for(int s=0;s<ns;s++){ 
            const int ii=i*ns+s; x[ii]=b[ii]-sum[s]; 
        }
    }
    __attribute__((hot)) 
    void fwd_subs_m_( int ns, const T* b, T* x ){
        //printf("SparseMatrix2::fwd_subs_m() [i=%i] n\", i);
        T sum[ns];
        for (int i=0; i<n; i++){
            fwd_subs_mi( i, ns, b, x );
        }
    }
    
    __attribute__((hot)) 
    void fwd_subs_T_mi( int i, int ns, const T* b, T* x ){
        //printf("SparseMatrix2::fwd_subs_m() [i=%i] \n", i);
        T sum[ns];
        for(int s=0;s<ns;s++){ sum[s]=0.0; }
        const int i0 = i0s [i];
        const int ni = nngs[i];
        for (int k=0; k<ni; k++){
            const int ik = i0+k;
            const int j  = inds[ik]*ns;
            const T   l  = vals[ik];  
            for(int s=0;s<ns;s++){ sum[s]+=l*x[j+s]; }
        }
        for(int s=0;s<ns;s++){ 
            const int ii=i*ns+s; x[ii]=b[ii]-sum[s]; 
        }
    }

    void fwd_subs_T_m_( int ns, const T* b, T* x ){
        //printf("SparseMatrix2::fwd_subs_m() [i=%i] ", i);
        T sum[ns];
        for (int i=n-1; i>=0; i--){
            fwd_subs_T_mi( ns, b, x );
        }
    }



    __attribute__((hot)) 
    void fwd_subs_m( int ns, const T* b, T* x ){
        printf("SparseMatrix2::fwd_subs_m() \n");
        T sum[ns];
        #pragma omp simd
        for (int i=0; i<n; i++){
            for(int s=0;s<ns;s++){ sum[s]=0.0; }
            const int i0 = i0s [i];
            const int ni = nngs[i];
            for (int k=0; k<ni; k++){
                const int ik = i0+k;
                const int j  = inds[ik]*ns;
                const T   l  = vals[ik];  
                for(int s=0;s<ns;s++){ sum[s]+=l*x[j+s]; }
            }
            for(int s=0;s<ns;s++){ 
                const int ii=i*ns+s; x[ii]=b[ii]-sum[s]; 
            }
        }
    }

    __attribute__((hot)) 
    void fwd_subs_T_m( int ns, const T* b, T* x ){
        //printf("SparseMatrix2::fwd_subs_m() [i=%i] ", i);
        T sum[ns];
        #pragma omp simd
        for (int i=n-1; i>=0; i--){
            for(int s=0;s<ns;s++){ sum[s]=0.0; }
            const int i0 = i0s [i];
            const int ni = nngs[i];
            for (int k=0; k<ni; k++){
                const int ik = i0+k;
                const int j  = inds[ik]*ns;
                const T   l  = vals[ik];  
                for(int s=0;s<ns;s++){ sum[s]+=l*x[j+s]; }
            }
            for(int s=0;s<ns;s++){ 
                const int ii=i*ns+s; x[ii]=b[ii]-sum[s]; 
            }
        }
    }

    T* reconstruct_dense( T* dens=0 ){
        if(dens==0) dens = new T[n*n];
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
                dens[i*n + j] = 0.0;
            }
            int i0 = i0s[i];
            int ni = nngs[i];
            for(int k=0; k<ni; k++){
                int j = inds[i0+k];
                dens[i*n + j] = vals[i0+k];
            }
        }
        return dens;
    }

    int fromDense( int n_, T* A, T tol, bool bRev=false ){
        n=n_;
        //printf( "SparseMatrix2::fromDense() n=%i \n", n );
        _realloc0(nngs, n,  0  );
        _realloc0(i0s,  n, -1 );
        ntot = 0;
        for(int i=0; i<n; i++){
            int ii = i; if(bRev)ii=n-i-1;
            int ni = countNonZero( n, A+ii*n, tol );
            //printf( "fromDense()[%i] ni=%i \n", i, ni );
            nngs[i]=ni;
            i0s [i]=ntot;
            ntot+=ni;
        }
        _realloc0(vals, ntot, (T)0 );
        _realloc0(inds, ntot, -1  );
        for(int i=0; i<n; i++){
            int ii = i; if(bRev)ii=n-i-1;
            int i0 = i0s[i];
            if(i0<0){ printf("ERROR in SparseMatrix::fromDense()[%i/%i] i0=%i \n", i, n, i0); exit(0); }
            int ni = exportNonZero( n, A+ii*n, tol, vals+i0, inds+i0 );
            if(ni!=nngs[i]){ printf("ERROR in SparseMatrix::fromDense()[%i/%i] exportNonZero().count(%i) != countNonZero().count(%i) => exit \n", i,n, ni, nngs[i] );  exit(0); }
        }
        return ntot;
    }

    int fromFwdSubT_( const T* L, bool bDo=true) {
        int ntot=0;
        for (int i=n-1; i>=0; i--){    
            int ni=0;
            for (int j=i+1; j<n; j++){
                const T l = L[j*n+i];
                if( fabs(l)>1e-12 ){ 
                    //printf("%6i ", j ); 
                    if(bDo){
                        inds[ntot+ni] = j;
                        vals[ntot+ni] = l;
                    }
                    ni++;
                }
            }
            nngs[i]=ni;
            i0s [i]=ntot;
            ntot+=ni;
        }
        return ntot;
    }
    int fromFwdSubT( int n_, const T* L) {
        n=n_;
        _realloc0(nngs, n, 0 );
        _realloc0(i0s,  n, -1 );
        ntot = fromFwdSubT_( L, false);
        _realloc0(vals, ntot, (T)0 );
        _realloc0(inds, ntot, -1  );
        fromFwdSubT_( L, true);
        return ntot;
    }

    int fprint_inds( const char* fname, const char* mode="w", const char* fmt="%6i " ){
        FILE* fout = fopen(fname,mode);
        for(int i=0; i<n; i++){
            int ni = nngs[i];
            int i0 = i0s [i];
            for(int k=0; k<ni; k++){
                fprintf( fout, fmt, inds[i0+k] );
            }
            fprintf( fout, "\n" );
        }
        fclose(fout);
        return 0;
    }

    int fprint_vals( const char* fname, const char* mode="w", const char* fmt="%20.20f " ){
        FILE* fout = fopen(fname,mode);
        for(int i=0; i<n; i++){
            int ni = nngs[i];
            int i0 = i0s [i];
            for(int k=0; k<ni; k++){
                fprintf( fout, fmt, vals[i0+k] );
            }
            fprintf( fout, "\n" );
        }
        fclose(fout);
        return 0;
    }

};


#endif

