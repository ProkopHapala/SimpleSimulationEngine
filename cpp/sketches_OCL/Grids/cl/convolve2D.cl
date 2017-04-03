
// http://www.cmsoft.com.br/opencl-tutorial/case-study-high-performance-convolution-using-opencl-__local-memory/


__kernel void blur2D_naive(
    __global float* I, 
    __global float* O
){
    const int ix = get_global_id (0)+1;
    const int iy = get_global_id (1)+1;
    const int nx = get_global_size(0)+2;
    
    int i = iy * nx + ix;
    
    // 3.9 Mticks
    O[i] =( I[i-nx-1] + I[i-nx] + I[i-nx+1] +
            I[i   -1] + I[i   ] + I[i   +1] +
            I[i+nx-1] + I[i+nx] + I[i+nx+1] ) * 0.11111111111;
}


__kernel void blur2D_Gauss(
    __global float* I, 
    __global float* O
){
    const int ix = get_global_id (0)+1;
    const int iy = get_global_id (1)+1;
    const int nx = get_global_size(0)+2;
    
    int i = iy * nx + ix;
    
    // 11.0 Mticks
    //O[i] =( 0.0625*I[i-nx-1] + 0.125*I[i-nx] + 0.0625*I[i-nx+1] +
    //        0.125 *I[i   -1] + 0.25 *I[i   ] + 0.125 *I[i   +1] +
    //        0.0625*I[i+nx-1] + 0.125*I[i+nx] + 0.0625*I[i+nx+1] );
    
    // 4.9 Mticks
    O[i] = 0.25  *  I[i] 
         + 0.0625*( I[i-nx-1] + I[i+nx-1] + I[i+nx+1] + I[i-nx+1] )
         + 0.125 *( I[i-nx]   + I[i-1]    + I[i+1]    + I[i+nx]   );
}


#define NBx 18
#define NBy 18


// seems to be slower than naive method
__kernel void blur2D_local(
    __global float* I, 
    __global float* O
){
    __local float L[NBx*NBy];
    const int2 iG  = (int2)(get_global_id  (0)+1 , get_global_id  (1)+1 );
    const int2 nG  = (int2)(get_global_size(0)+2 , get_global_size(1)+2 );
    const int2 iL  = (int2)(get_local_id   (0)+1 , get_local_id   (1)+1 );
    const int2 nL  = (int2)(get_local_size (0)+2 , get_local_size (1)+2 );
    const int2 iGR = (int2)(get_group_id   (0)   , get_group_id   (1)   );
    
    int ig = iG.y*nG.x + iG.x;
    int il = iL.y*nL.x + iL.x;
    
    L[il] = I[ig];  // each thread copy one pixel
    
    // some threads has to copy one more for boundary (halo)
    switch( get_local_id(1) ){
        case 4: 
        switch( get_local_id(0) ){
            case 0: L[        0      ] = I[ nG.x* get_group_id(1)*get_local_size(1)          + get_group_id(0)*get_local_size(0)         ]; break;
            case 1: L[         NBx-1 ] = I[ nG.x* get_group_id(1)*get_local_size(1)          + get_group_id(0)*get_local_size(0)+(NBx-1) ]; break;
            case 2: L[ (NBy-1)*NBx   ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+(NBy-1)) + get_group_id(0)*get_local_size(0)         ]; break;
            case 3: L[ NBy*    NBx-1 ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+(NBy-1)) + get_group_id(0)*get_local_size(0)+(NBx-1) ]; break;
        }
        case 0: L[               iL.x    ] = I[ nG.x* get_group_id(1)*get_local_size(1)                   + iG.x                                        ]; break;
        case 1: L[ NBx*(NBy-1) + iL.x    ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+(NBy-1)         ) + iG.x                                        ]; break;
        case 2: L[ NBx*iL.x              ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+get_local_id(0) ) +  get_group_id(0)*get_local_size(0)          ]; break;
        case 3: L[ NBx*iL.x    + (NBx-1) ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+get_local_id(0) ) + (get_group_id(0)*get_local_size(0)+(NBx-1)) ]; break;
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    const float renorm = 1.0/9.0;
    O[ig] =( L[il-NBx-1] + L[il-NBx] + L[il-NBx+1] +
             L[il    -1] + L[il    ] + L[il    +1] +
             L[il+NBx-1] + L[il+NBx] + L[il+NBx+1] ) / 9.0;
}



// more efficient algorithm - scanline
//
//   l1 = ***********************
//   l2 = *+++++++++++++++++++++*
//   l3 = ***********************
//
//   we swap __local memory buffers l1,l2,l3  only l3 is read from __global 

// seems to be slower than naive method
__kernel void blur2D_scanline(
    int nx, int ny,
    __global float* I, 
    __global float* O
){
    __local float L[32*3];
    __local float *Lm = L;
    __local float *L0 = L+32;
    __local float *Lp = L+64;
    
    const int il = get_local_id(0);
    
    Lm[il] = I[      get_global_id(0) ];
    L0[il] = I[ nx + get_global_id(0) ];
    //Lp[iL] = I[ (iy+1)*nx +get_global_id(0)]
    barrier(CLK_LOCAL_MEM_FENCE);    
    for(int iy=1; iy<(ny-1); iy++ ){
        Lp[il] = I[ (iy+1)*nx +get_global_id(0)];
        barrier(CLK_LOCAL_MEM_FENCE);   
        O[iy*nx +get_global_id(0)] = 
            ( Lm[il-NBx-1] + Lm[il-NBx] + Lm[il-NBx+1] +
              L0[il    -1] + L0[il    ] + L0[il    +1] +
              Lp[il+NBx-1] + Lp[il+NBx] + Lp[il+NBx+1] ) * 0.11111111111;
        barrier(CLK_LOCAL_MEM_FENCE);
        __local float *Ltmp=Lm; Lm=L0; L0=Lp; Lp=Ltmp;
        
    }

}


__kernel void blur2D_scanline_priv(
    int nx, int ny,
    __global float* I, 
    __global float* O
){ 
    int ig    = get_global_id(0);
    float3 Lm = (float3)( I[ig-1], I[ig], I[ig+1] );  ig += nx;
    float3 L0 = (float3)( I[ig-1], I[ig], I[ig+1] ); 
    for(int iy=1; iy<(ny-1); iy++ ){
        ig += nx;
        float3 Lp= (float3)( I[ig-1], I[ig], I[ig+1] );  
        O[ig-nx] = 
            ( Lm.x + Lm.y + Lm.z +
              L0.x + L0.y + L0.z +
              Lp.x + Lp.y + Lp.z ) * 0.11111111111;
        Lm=L0; L0=Lp; 
    }
}





