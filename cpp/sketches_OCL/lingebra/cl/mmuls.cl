

//  ==========   mmul_elem       

__kernel void mmul_elem(
    const int N,
    __global float* A,
    __global float* B,
    __global float* C
){
    int k;
    int i = get_global_id(0);
    int j = get_global_id(1);
    float tmp;
    if ((i < N) && (j < N)) {
        tmp = 0.0;
        for (k = 0; k < N; k++) tmp += A[i*N+k] * B[k*N+j];
        C[i*N+j] = tmp;
    }
};

//  ==========   mmul_row      

__kernel void mmul_row(
    const int N,
    __global float* A,
    __global float* B,
    __global float* C
){
    int k, j;
    int i = get_global_id(0);
    float tmp;
    if (i < N) {
        for (j = 0; j < N; j++) {
            tmp = 0.0;
            for (k = 0; k < N; k++) tmp += A[i*N+k] * B[k*N+j];
            C[i*N+j] = tmp;
        }
    }
};

//  ==========   mmul_row_priv     

__kernel void mmul_row_priv(
    const int N,
    __global float* A,
    __global float* B,
    __global float* C
){
    int k, j;
    int i = get_global_id(0);
    float Awrk[1024];
    float tmp;
    if (i < N) {
        for (k = 0; k < N; k++) Awrk[k] = A[i*N+k];
        for (j = 0; j < N; j++) {
            tmp = 0.0f;
            for (k = 0; k < N; k++) tmp += Awrk[k] * B[k*N+j];
            C[i*N+j] = tmp;
        }
    }
};

//  ==========   mmul_row_priv_block

__kernel void mmul_row_priv_block(
    const int N,
    __global float* A,
    __global float* B,
    __global float* C,
    __local float* Bwrk
){
    int k, j;
    int i    = get_global_id(0);
    int iloc = get_local_id(0);
    int nloc = get_local_size(0);
    float Awrk[1024];
    float tmp;
    if (i < N) {
        for (k = 0; k < N; k++) Awrk[k] = A[i*N+k];
        for (j = 0; j < N; j++) {
            barrier(CLK_LOCAL_MEM_FENCE);
            for (k = iloc; k < N; k += nloc) Bwrk[k] = B[k*N+j];
            barrier(CLK_LOCAL_MEM_FENCE);
            tmp = 0.0f;
            for (k = 0; k < N; k++) tmp += Awrk[k] * Bwrk[k];
            C[i*N+j] = tmp;
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }
}

// from http://gpgpu-computing4.blogspot.cz/2009/10/matrix-multiplication-3-opencl.html
#define BLOCK_SIZE 16

__kernel void mmul_local(
          __global float* Aik, 
          __global float* Bkj, 
          __global float* Cij, 
          __const int ni, 
          __const int nj,
          __const int nk
){
    //   WARRNING : interchange of  i  and  j  dimension  lower the performance >2x on my nV GT275 GPU    
    int gj = get_global_id(0);    int gi = get_global_id(1); 
    int bj = get_group_id(0);     int bi = get_group_id(1);  // Block index
    int tj = get_local_id(0);     int ti = get_local_id(1);  // Thread index
    int oj = bi*BLOCK_SIZE;       int oi = bj*BLOCK_SIZE; 
    float Csub =0; 
    __local float As   [BLOCK_SIZE][BLOCK_SIZE];
    __local float Bs   [BLOCK_SIZE][BLOCK_SIZE];
    for (int ok = 0; ok < nk; ok += BLOCK_SIZE )   {
        As[ti][tj] = Aik[ nk*(gi   ) + tj + ok ];   // A[i][k]
        Bs[ti][tj] = Bkj[ nj*(ti+ok) + gj ];        // B[k][j]
        barrier(CLK_LOCAL_MEM_FENCE);
        #pragma unroll
        for (int k = 0; k < BLOCK_SIZE; ++k) Csub += As[ti][k] * Bs[k][tj];
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    Cij[ nj * ( gi ) + gj ] = Csub;
}

//-------------------------------------------------------------
//
//  PROGRAM: Blocked Matrix Multipliplication kernel
//
//  PURPOSE: Computes an element of the proudct matrix
//
//              C = A * B
//
//           Using the well known blocked algorithm.  
//
//           To derive this algorithm, start with the naive
//           triply nested loop algorithm with a dot product 
//           for each element of C.  Decompose each loop 
//           into blocks of size blcksz.  This gives you 6
//           nested loops with three loops over blocks
//           and three loops over indices inside the blocks.
// 
//           Rearrange the loops to put the 3 loops over blocks 
//           at the outermost loops of the loop nest.  You'll
//           see that the three "inner" loops are just the 
//           regular matrix product between blocks.
//
//           The algorithms is simple.  Keeping all the indices
//           straight is not.  We will use the following 
//           conventions:
//
//             i,j,k            ... indices of full, global matrices 
//             Iblk, Jblk, Kblk ... indices of matrix blocks
//             iloc, jloc, kloc ... indices inside blocks
//                 
//  HISTORY: Written by Tim Mattson, November 2013 
//           Updated by Simon McIntosh-Smith, August 2014 
//
//  LICENSE: This work is licensed under the Creative Commons
//           Attribution 4.0 International License.
//           To view a copy of this license, visit
//           http://creativecommons.org/licenses/by/4.0/
//           or send a letter to:
//              Creative Commons,
//              444 Castro Street, Suite 900,
//              Mountain View, California, 94041, USA.
//
//-------------------------------------------------------------

// It turns out that the compiler generates much better code if we "hardwire" this block size.
// 16 works well for an NVIDIA  GPU, 32 works well for a CPU

#define blksz 16

__kernel void mmul_block(
    const unsigned int             N,
    __global const float* restrict A,
    __global const float* restrict B,
    __global       float* restrict C,
    __local        float* restrict Awrk,
    __local        float* restrict Bwrk
){
    int kloc, Kblk;
    float Ctmp=0.0f;

    const int i    = get_global_id(0); const int j    = get_global_id(1); // This work-item will compute element C(i,j)
    const int Iblk = get_group_id(0);  const int Jblk = get_group_id(1);  // Element C(i,j) is in block C(Iblk,Jblk)
    const int iloc = get_local_id(0);  const int jloc = get_local_id(1);  // C(i,j) is element C(iloc, jloc) of block C(Iblk, Jblk)
    const int Num_BLK = N/blksz;                                          // The number of blocks are the same in each dimension

    // Setup the upper-left-corner (base address) for the A and
    // B blocks plus the increments to advance base addresses as
    // we loop over blocks
            
    const int Ainc  = blksz;
    const int Binc  = blksz*N;
          int Abase = Iblk*N*blksz;  
          int Bbase = Jblk*blksz;
    
    // C(Iblk,Jblk) = (sum over Kblk) A(Iblk,Kblk)*B(Kblk,Jblk)
    for (Kblk = 0;  Kblk<Num_BLK;  Kblk++){
       // Load A(Iblk,Kblk) and B(Kblk,Jblk) into local memory.
       // Each work-item loads a single element of the two blocks
       // which are shared with the entire work-group.
       Awrk[jloc*blksz+iloc] = A[Abase+jloc*N+iloc];
       Bwrk[jloc*blksz+iloc] = B[Bbase+jloc*N+iloc];
       barrier(CLK_LOCAL_MEM_FENCE);
       // Compute dot products over local blocks to find
       // the contribution to C(i,j) from this block
       #pragma unroll
       for (kloc=0; kloc<blksz; kloc++)  Ctmp += Awrk[jloc*blksz+kloc] * Bwrk[kloc*blksz+iloc];
       barrier(CLK_LOCAL_MEM_FENCE);
       Abase += Ainc;
       Bbase += Binc;
    }
 
    // update global C matrix 
    C[j*N+i] = Ctmp;

}



/*
__kernel void  mmul_block(
                const unsigned int             N,
                __global const float* restrict A,
                __global const float* restrict B,
                __global       float* restrict C,
                __local        float* restrict Awrk,
                __local        float* restrict Bwrk)
{
    int kloc, Kblk;
    float Ctmp=0.0f;

    //  This work-item will compute element C(i,j)
    const int i = get_global_id(0);
    const int j = get_global_id(1);

    // Element C(i,j) is in block C(Iblk,Jblk)
    const int Iblk = get_group_id(0);
    const int Jblk = get_group_id(1);

    // C(i,j) is element C(iloc, jloc) of block C(Iblk, Jblk)
    const int iloc = get_local_id(0);
    const int jloc = get_local_id(1);

    // The number of blocks are the same in each dimension
    const int Num_BLK = N/blksz;

    // Setup the upper-left-corner (base address) for the A and
    // B blocks plus the increments to advance base addresses as
    // we loop over blocks
          int Abase = Iblk*N*blksz;    
    const int Ainc  = blksz;

          int Bbase = Jblk*blksz;
    const int Binc  = blksz*N;


    // C(Iblk,Jblk) = (sum over Kblk) A(Iblk,Kblk)*B(Kblk,Jblk)
    for (Kblk = 0;  Kblk<Num_BLK;  Kblk++)
    {
       // Load A(Iblk,Kblk) and B(Kblk,Jblk) into local memory.
       // Each work-item loads a single element of the two blocks
       // which are shared with the entire work-group.

       Awrk[jloc*blksz+iloc] = A[Abase+jloc*N+iloc];
       Bwrk[jloc*blksz+iloc] = B[Bbase+jloc*N+iloc];

       barrier(CLK_LOCAL_MEM_FENCE);

       // Compute dot products over local blocks to find
       // the contribution to C(i,j) from this block
       #pragma unroll
       for (kloc=0; kloc<blksz; kloc++)
          Ctmp += Awrk[jloc*blksz+kloc] * Bwrk[kloc*blksz+iloc];

       barrier(CLK_LOCAL_MEM_FENCE);
       Abase += Ainc;
       Bbase += Binc;
    }
 
    // update global C matrix 
    C[j*N+i] = Ctmp;

}
*/





