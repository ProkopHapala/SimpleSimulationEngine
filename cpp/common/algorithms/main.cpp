
#include "testUtils.h"
#include "arrayAlgs.h"

const int n = 10;
//static int A[n] = {6,10,13,5,8,3,2,25,4,11};
//static int out[n];

#define TYPE int

static TYPE  A[n] = { 6, 10, 13, 5, 8, 3, 2, 25, 4, 11 };
static TYPE out[n];
static int permut[n];

int main(){

/*
	printf("A:      " ); printArray( n, &A[0]  );
    quickSort_inplace<int>( A, 0, n );
    printf("A:      " ); printArray( n, &A[0]  );
*/

	indexArray<int>( n, permut );

    printf("A:      " ); printArray( n, &A[0]  );
	printf("permut: " ); printArray( n, permut );

	quickSort<TYPE>( A, permut,      0, n );
	permute  <TYPE>( permut, A, out, 0, n );

	printf("out:    " ); printArray( n, out    );
	printf("permut: " ); printArray( n, permut );

}
