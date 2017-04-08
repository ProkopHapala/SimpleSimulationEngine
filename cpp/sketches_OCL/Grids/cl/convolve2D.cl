
// http://www.cmsoft.com.br/opencl-tutorial/case-study-high-performance-convolution-using-opencl-__local-memory/

/*
#### performance [CPU ticks/pixel] (lower is better):
 - CPU-naive            9.5
 - GPU-naive            1.64
 - GPU-local            2.56
 - GPU-scanline private 7.35

The test was done on [my laptop](https://www.alza.cz/lenovo-ideapad-y700-15isk?kampan=adw1_notebooky_dsa-produkty&gclid=Cj0KEQjwiI3HBRDv0q_qhqXZ-N4BEiQAOTiCHhVp7gaLjuPEOwdkyfm02VBReROW9J9Kp3sZvFjyEYUaAi8q8P8HAQ) with CPU Intel Core i7 6700HQ Skylake and GPU nVidia 960M by running the kernels 64x/frame on floating point array of 256x256 pixels. 
*/

    __kernel void blur2D_naive(
        __global float* I, 
        __global float* O
    ){
        const int ix = get_global_id (0)+1;
        const int iy = get_global_id (1)+1;
        const int nx = get_global_size(0)+2;
        
        int i = iy * nx + ix;
        
        O[i] =( I[i-nx-1] + I[i-nx] + I[i-nx+1] +
                I[i   -1] + I[i   ] + I[i   +1] +
                I[i+nx-1] + I[i+nx] + I[i+nx+1] ) * 0.11111111111f;
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
        O[i] =( 0.0625*I[i-nx-1] + 0.125*I[i-nx] + 0.0625*I[i-nx+1] +
                0.125 *I[i   -1] + 0.25 *I[i   ] + 0.125 *I[i   +1] +
                0.0625*I[i+nx-1] + 0.125*I[i+nx] + 0.0625*I[i+nx+1] );
        //O[i] = 0.25  *  I[i] 
        //     + 0.0625*( I[i-nx-1] + I[i+nx-1] + I[i+nx+1] + I[i-nx+1] )
        //     + 0.125 *( I[i-nx]   + I[i-1]    + I[i+1]    + I[i+nx]   );
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
        
        switch( get_local_id(1) ){ // some threads copy one more of boundary (halo) pixels
            case 4: 
            switch( get_local_id(0) ){ // copy corner points
                case 0: L[        0      ] = I[ nG.x* get_group_id(1)*get_local_size(1)          + get_group_id(0)*get_local_size(0)         ]; break; // upper-left
                case 1: L[         NBx-1 ] = I[ nG.x* get_group_id(1)*get_local_size(1)          + get_group_id(0)*get_local_size(0)+(NBx-1) ]; break; // upper-right
                case 2: L[ (NBy-1)*NBx   ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+(NBy-1)) + get_group_id(0)*get_local_size(0)         ]; break; // lower-left
                case 3: L[ NBy*    NBx-1 ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+(NBy-1)) + get_group_id(0)*get_local_size(0)+(NBx-1) ]; break; // lower-rigth
            }
            // copy border lines 
            case 0: L[               iL.x    ] = I[ nG.x* get_group_id(1)*get_local_size(1)                   + iG.x                                        ]; break; // top    line
            case 1: L[ NBx*(NBy-1) + iL.x    ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+(NBy-1)         ) + iG.x                                        ]; break; // botton line
            case 2: L[ NBx*iL.x              ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+get_local_id(0) ) +  get_group_id(0)*get_local_size(0)          ]; break; // left   line
            case 3: L[ NBx*iL.x    + (NBx-1) ] = I[ nG.x*(get_group_id(1)*get_local_size(1)+get_local_id(0) ) + (get_group_id(0)*get_local_size(0)+(NBx-1)) ]; break; // right  line
        } // each thread coppied at max. 2 pixels
        
        int ig = iG.y*nG.x + iG.x;
        int il = iL.y*nL.x + iL.x;
        L[il] = I[ig];             // each thread copy his pixel to local memory
        
        barrier(CLK_LOCAL_MEM_FENCE);
        
        O[ig] =( L[il-NBx-1] + L[il-NBx] + L[il-NBx+1] +
                 L[il    -1] + L[il    ] + L[il    +1] +
                 L[il+NBx-1] + L[il+NBx] + L[il+NBx+1] ) / 9.0f;
    }
    
    #define nTiles 16
    #define NBx 18
    #define NBy 18 
    #define copy_tile(event,ig0,I,L) { int ig_=ig0; int il_=0; for(int i=0; i<NBy; i++){   event = async_work_group_copy( L+il_, I+ig_, NBx, event ); ig_+=nx; il_+=NBx; } }
    // https://streamcomputing.eu/blog/2014-06-19/using-async_work_group_copy-on-2d-data/
    __kernel void blur2D_local_async(
        __global float* I, 
        __global float* O
    ){
        const int nx = get_global_size(0)+2;        
        __local float LI[NBx*NBy*2];
        int iL0 = 0;
        int iL1 = NBx*NBy;        
        event_t event = 0;
        int ig0 = get_group_id(0)*get_local_size(0);
        copy_tile(event,ig0,I,LI);
        for( int it=0; it<nTiles; it++ ){
            int ig   = ig0 + (get_local_id(1)+1)*nx  + get_local_id(0)+1;
            int il   =       (get_local_id(1)+1)*NBx + get_local_id(0) + iL0;
            ig0     += get_local_size(1)*nx;
            event_t event_ = 0;
            copy_tile(event_,ig0,I,LI+iL1);
            wait_group_events(1, &event);
            //barrier(CLK_LOCAL_MEM_FENCE);
            O[ig] =( LI[il-NBx] + LI[il-NBx+1] + LI[il-NBx+2] +
                     LI[il    ] + LI[il    +1] + LI[il    +2] +
                     LI[il+NBx] + LI[il+NBx+1] + LI[il+NBx+2] ) * 0.11111111111;
            int iLtmp=iL0; iL0=iL1; iL1=iLtmp;
            event = event_;
        }
    }
    
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
    
    
    
     #define NB 66
    __kernel void blur2D_scanline_async(
        int nx, int ny,
        __global float* I, 
        __global float* O
    ){
        __local float  L[NB*4];
        int i0=0;
        int i1=NB;
        int i2=NB*2;
        int i3=NB*3;
        event_t event = 0;
        int ig0 = get_group_id(0)*get_local_size(0);
        event = async_work_group_copy(  L     , I+ig0, NB, event );    ig0 += nx;
        event = async_work_group_copy(  L+NB  , I+ig0, NB, event );    ig0 += nx;   
        event = async_work_group_copy(  L+NB*2, I+ig0, NB, event );    ig0 += nx;
        const int il = get_local_id(0);
        int ig = get_global_id(0)+1;
        for(int iy=1; iy<(ny-2); iy++ ){
            wait_group_events(1, &event);
            event = async_work_group_copy(  L+i3, I+ig0, NB, event ); ig0 += nx;
            ig += nx;
            O[ig] =  
                ( L[i0+il] + L[i0+il+1] + L[i0+il+2] +
                  L[i1+il] + L[i1+il+1] + L[i1+il+2] +
                  L[i2+il] + L[i2+il+1] + L[i2+il+2] ) * 0.11111111111;
            __local float *Ltmp;
            int itmp=i0; i0=i1; i1=i2; i2=i3; i3=itmp;
        }
    }
    
    
    

    __kernel void blur2D_scanline_priv(
        int nx, int ny,
        __global float* I, 
        __global float* O
    ){ 
        int ig    = get_global_id(0)+1;
        float3 Lm = (float3)( I[ig-1], I[ig], I[ig+1] );  ig += nx;
        float3 L0 = (float3)( I[ig-1], I[ig], I[ig+1] ); 
        for(int iy=1; iy<(ny-1); iy++ ){
            ig += nx;
            float3 Lp= (float3)( I[ig-1], I[ig], I[ig+1] );  
            O[ig-nx] = 
                ( Lm.x + Lm.y + Lm.z +
                  L0.x + L0.y + L0.z +
                  Lp.x + Lp.y + Lp.z ) * 0.11111111111f;              
            Lm=L0; L0=Lp; 
        }
    }


