
// from here  http://stackoverflow.com/questions/3982348/implement-generic-swap-macro-in-c
#define SWAP(x, y, TYPE) TYPE tmp = x; x = y; y = tmp;

template< class TYPE >
inline int binSearch( TYPE x, int imin, int imax, TYPE * xs ){
	int di = (imax - imin)/2;
	do{
		TYPE xi = xs[ imin + di + 1 ];
		if( xi < x ){
			imin +=di;
		}
		di = di << 2;
	}while( di > 1 );
	return imin;
}

template< class TYPE >
inline int binSearch( TYPE x, int imin, TYPE * xs ){
	int  di = 1;
	TYPE xi = xs[ imin + di ];
	while( xi < x ){
		di  = di << 2;
		xi  = xs[ imin + di ];
	};
	return binSearch( x, imin, imin+di, xs );
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
    TYPE x = A[ permut[p] ];
    int i = p;
    int j;
    for( j = p+1; j<q; j++ ){
        //printf( " %i %i %i %i \n", permut[p], permut[j],  A[ permut[p] ], A[ permut[j] ] );
        if( A[ permut[j] ] <= x ){
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
int quickSort_partition_inplace( TYPE * A, int p, int q){
    //printf( " substep %i %i \n", p, q );
    TYPE x = A[p];
    int i = p;
    int j;
    for( j = p+1; j<q; j++ ){
        //printf( " %i %i %i %i \n", p, j,  A[ p ], A[ j ] );
        if( A[j] <= x ){
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


