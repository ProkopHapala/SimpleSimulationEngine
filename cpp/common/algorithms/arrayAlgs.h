
#ifndef  arrayAlgs_h
#define  arrayAlgs_h

// from here  http://stackoverflow.com/questions/3982348/implement-generic-swap-macro-in-c
//#define SWAP(x, y, TYPE) {TYPE tmp = x; x = y; y = tmp;}

#include "macroUtils.h"
#include "stdio.h"

/*
template<class TYPE>
inline TYPE * insertN( int i, int nTo, int nFrom, TYPE * to, TYPE * from ){
    TYPE * to_ = new TYPE[nTo+nFrom];
    for(int i=0; i< ; i++){   }
}

template<class TYPE>
inline int insertN( int i0, int len, TYPE * x, TYPE * x ){

}
*/

template< typename T, typename CondFunc >
inline int prune( int n, T* arr, CondFunc cond ){
    if(n==1){ if( cond( arr[0] ) ) {return 0;} else {return 1; } };
    int nlast = n-1;
    for(int i=0; i<nlast; i++){
        if( cond( arr[i] ) ){
            while( cond( arr[nlast] ) ){
                nlast--;
                if( i>=nlast ) return 0;
            }
            arr[i] = arr[nlast];
            nlast--;
        }
    }
    return nlast+1;
}


inline int objects2cells( int nobj, int ncell, int* obj2cell, int* cellNs, int* cellI0s, int* permut ){
    //printf("DEBUG 1 objects2cells \n");
    int nmax = 0;
    for(int k=0; k<ncell; k++ ){ cellNs[k]=0; }
    for(int i=0; i<nobj; i++ ){
        cellNs[ obj2cell[i] ]++;
    }
    //printf("DEBUG 2 objects2cells \n");
    int ntot=0;
    for(int k=0; k<ncell; k++ ){
        cellI0s[k]=ntot;
        int& ni = cellNs[k];
        ntot   += ni;
        if( ni>nmax ) nmax=ni;
        ni=0;
    }
    //printf("DEBUG 3 objects2cells \n");
    for(int i=0; i<nobj; i++){
        printf( "[%i] ", i );
        int k =  obj2cell[i];
        int j =  cellI0s[k] + cellNs[k];
        printf( " k %i j %i | nobj %i \n", k, j, nobj );
        permut[j] = i;
        cellNs[k]++;
    }
    //printf("DEBUG 4 objects2cells \n");
    return nmax;
}

struct I0n{
    int i0,n;
};

inline int objects2cells( int nobj, int ncell, int* obj2cell, I0n* cells, int* permut ){
    int nmax = 0;
    for(int i=0; i<nobj; i++ ){
        cells[ obj2cell[i] ].n++;
    }
    int ntot=0;
    for(int k=0; k<ncell; k++ ){
        cells[k].i0 = ntot;
        int& ni     = cells[k].n;
        ntot       += ni;
        if( ni>nmax ) nmax=ni;
        ni=0;
    }
    for(int i=0; i<nobj; i++){
        int k =  obj2cell[i];
        int j =  cells[k].i0 + cells[k].n;
        permut[j] = i;
        cells[k].n++;
    }
    return nmax;
}


template<typename T>
int insertSort( int n, int* permut, T* data ){
    //https://en.wikipedia.org/wiki/Insertion_sort
    int niter=0;
    int i=1;
    for(int i=1; i<n; i++){
        int ix = permut[i];
        const T& x = data[ix];
        int j=i-1;
        while( data[permut[j]] > x && (j>=0) ){
            permut[j+1] = permut[j];
            j=j-1; // backward iteration is not that great, but we already have it in cache
            niter++;
        }
        niter++;
        permut[j+1] = ix;
    }
    return niter;
}

template<typename T>
int insertSort_reverse( int n, int* permut, T* data ){
    //https://en.wikipedia.org/wiki/Insertion_sort
    int niter=0;
    int i=1;
    for(int i=1; i<n; i++){
        int ix = permut[i];
        const T& x = data[ix];
        int j=i-1;
        while( data[permut[j]] > x && (j>=0) ){
            permut[j+1] = permut[j];
            j=j-1; // backward iteration is not that great, but we already have it in cache
            niter++;
        }
        niter++;
        permut[j+1] = ix;
    }
    return niter;
}


//1 2 1 2.0000
template< class TYPE >
inline int binSearchBetween( TYPE x, int imax, TYPE * xs ){
	int di   = imax;
	int imin = 0;
	do{
        di  = di>>1;
        int  i  = imin + di;
		if( xs[i] < x ){ imin=i; di = (imax-i); }
		//printf("binSearchBetween %i %i %i %f %f \n", imin, i, di, xs[i], x );
	}while( di > 1 );
	//printf( " %f < %f < %f\n", xs[imin], x, xs[imin+1] );
	return imin;
}

template< class TYPE >
inline int binSearchFrom( TYPE x, int n, TYPE * xs ){
	int   i = 1;
	do{
        i  = i << 1;
		if(i>n){ if(xs[n]>x){ return binSearchBetween( x, n, xs ); }else{ return -1; } }
		//printf("binSearchForm    %i %f %f \n", i, xs[i], x);
	}while( xs[i] < x );
	return binSearchBetween( x, i, xs );
}

template< class TYPE >
inline void indexArray( int n, int * permut ){ for( int i=0; i<n; i++ ){ permut[i] = i; } }

template< class TYPE >
inline void permute( int * permut, TYPE * input, TYPE * output, int p, int q ){
	for( int i=p; i<q; i++ ){
		output[ i ] = input[ permut[i] ];
        //output[ permut[i] ] = input[ i ];
	}
}

template< class TYPE >
int quickSort_partition( TYPE * A, int * permut, int p, int q){
    //printf( " quickSort_partition %i %i \n", p, q );
    const TYPE& x = A[ permut[p] ];
    int i = p;
    int j;
    for( j = p+1; j<q; j++ ){
        //printf( " %i %i %i %i \n", permut[p], permut[j],  A[ permut[p] ], A[ permut[j] ] );
        if( A[ permut[j] ] < x ){
            i++;
            SWAP( permut[i], permut[j], int );
        }
    }
    SWAP( permut[i], permut[p], int );
    return i;
}

template< class TYPE >
void quickSort( TYPE * A, int * permut, int p, int q){
    //printf( " quickSort %i %i \n", p, q );
    if( p < q){
        int r = quickSort_partition<TYPE>( A, permut, p, q );
        quickSort<TYPE>( A, permut, p  , r );
        quickSort<TYPE>( A, permut, r+1, q );
    }
}

template< class TYPE >
int quickSort_partition_reverse( TYPE * A, int * permut, int p, int q){
    //printf( " quickSort_partition %i %i \n", p, q );
    const TYPE& x = A[ permut[p] ];
    int i = p;
    int j;
    for( j = p+1; j > q; j++ ){
        //printf( " %i %i %i %i \n", permut[p], permut[j],  A[ permut[p] ], A[ permut[j] ] );
        if( A[ permut[j] ] < x ){
            i++;
            SWAP( permut[i], permut[j], int );
        }
    }
    SWAP( permut[i], permut[p], int );
    return i;
}

template< class TYPE >
void quickSort_reverse( TYPE * A, int * permut, int p, int q){
    //printf( " quickSort %i %i \n", p, q );
    if( p > q){
        int r = quickSort_partition<TYPE>( A, permut, p, q );
        quickSort<TYPE>( A, permut, p  , r );
        quickSort<TYPE>( A, permut, r+1, q );
    }
}








template< class TYPE >
int quickSort_partition_inplace( TYPE * A, int p, int q){
    //printf( " substep %i %i \n", p, q );
    TYPE x = A[p];
    int i = p;
    int j;
    for( j = p+1; j<q; j++ ){
        //printf( " %i %i %i %i \n", p, j,  A[ p ], A[ j ] );
        if( A[j] < x ){
            i++;
            SWAP( A[i], A[j], TYPE );
        }
    }
    SWAP( A[i], A[p], TYPE );
    return i;
}

template< class TYPE >
void quickSort_inplace( TYPE * A, int p, int q){
    if(p<q){
        int r = quickSort_partition_inplace<TYPE>( A, p, q );
        quickSort_inplace<TYPE>( A, p  , r );
        quickSort_inplace<TYPE>( A, r+1, q );
    }
}

template<typename T>
void sort_permut( int n, int* permut, T* data, bool bReverse=false, bool almostSorted=false ){
    /*
    if(bReverse){
        if(almostSorted){ printf( "insertSort_reverse \n" ); insertSort_reverse( n, permut, data    );  }
        else            { printf( "quickSort_reverse  \n" ); quickSort_reverse ( data, permut, 0, n ); };
    }else{
        if(almostSorted){ printf( "insertSort \n" );          insertSort( n, permut, data   );  }
        else            { printf( "quickSort  \n" );          quickSort( data, permut, 0, n );  };
    }
    */
    if(bReverse){
        if(almostSorted){ insertSort_reverse( n, permut, data    ); }
        else            { quickSort_reverse ( data, permut, 0, n ); };
    }else{
        if(almostSorted){ insertSort( n, permut, data   );  }
        else            { quickSort( data, permut, 0, n );  };
    }
}


#endif
