
__kernel void blur2D_naive(
    __global float* I, 
    __global float* O
){
    const int ix = get_global_id (0)+1;
    const int iy = get_global_id (1)+1;
    const int nx = get_global_size(0)+2;
    
    int i = iy * nx + ix;
    
    const float renorm = 1.0/9.0;
    O[i] =( I[i-nx-1] + I[i-nx] + I[i-nx+1] +
            I[i   -1] + I[i   ] + I[i   +1] +
            I[i+nx-1] + I[i+nx] + I[i+nx+1] ) / 9.0;
}


__kernel void blur2D_local(
    __global float* I, 
    __global float* O
){
    __local float L[18*18];
    const int ix = get_global_id (0)+1;
    const int iy = get_global_id (1)+1;
    const int lx = get_local_id  (0)+1;
    const int ly = get_local_id  (1)+1;
    const int nx = get_global_size(0)+2;
    const int ny = get_global_size(1)+2;
    
    //printf("%i %i\n", get_global_size(0), get_global_size(1) );
    printf("(%i,%i) (%i,%i)\n", lx, ly, get_local_size(0), get_local_size(1) );
    int i = iy*nx+ix;
    //int l = ly*18+lx;
    int l = ly*16+lx;
    
    L[l] = I[i];
    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    const float renorm = 1.0/9.0;
    //O[i] =( L[l-18-1] + L[l-18] + L[l-18+1] +
    //        L[l   -1] + L[l   ] + L[l   +1] +
    //        L[l+18-1] + L[l+18] + L[l+18+1] ) / 9.0;
    O[i] = l;
    
}



