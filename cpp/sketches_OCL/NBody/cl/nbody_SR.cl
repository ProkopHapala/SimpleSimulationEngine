
#define R2SAFE 1.0e-2
#define F2MAX  10.0

#define R_MAX  1.8
#define R2MAX  3.24

float2 pair_force( float2 p1, float2 p2 ){
    float2 d   = p2 - p1;
    float  r2  = dot(d,d);
    if( r2 > R2MAX ) return (float2) (0.0f, 0.0f);
    float ir2 = 1/( r2 + R2SAFE );
    float fr  = (0.7-ir2)*(R2MAX-r2);
    return d * fr;
}

__kernel void NBody_force_naive(
    unsigned int num,
    __global float2* pos, 
    __global float2* force
){
    const int i = get_global_id (0);
    float2 p = pos[i];
    float2 f = (float2) (0.0f, 0.0f);
    for (int j=0; j < num; j++){
        float2 pj = pos[j];
        f += pair_force( p, pj );
    }
    force[i] = f;
}


__kernel void NBody_force(
    unsigned int num,
    __global float2* pos, 
    __global float2* force
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    __local  float2 pos_shared[64];
    float2 p = pos[iG];
    float2 f = (float2) (0.0f, 0.0f);
    for (int i0=0; i0<num; i0+= nL ){ 
        pos_shared[iL] = pos[i0 + iL]; 
        barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
        for (int j=0; j<nL; j++){
            float2 pj  = pos_shared[j];
            f         += pair_force( p, pj );
        }
        barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
    }
    force[iG] = f;
}

// my local memory is 49152 bytes 49kb
// http://stackoverflow.com/questions/31197564/is-cl-device-local-mem-size-for-the-entire-device-or-per-work-group
//
/*



TL;DR: Per single processing unit, hence also the maximum allotable to a work unit.

This value is the amount of local memory available on each compute unit in the device. Since a work-group is assigned to a single compute unit, this is also the maximum amount of local memory that any work-group can have.

For performance reasons on many GPUs, it is usually desirable to have multiple work-groups running on each compute unit concurrently (to hide memory access latency, for example). If one work-group uses all of the available local memory, the device will not be able to schedule any other work-groups onto the same compute unit until it has finished. If possible, it is recommended to limit the amount of local memory each work-group uses (to e.g. a quarter of the total local memory) to allow multiple work-groups to run on the same compute unit concurrently.

*/

__kernel void force_Tiled(
    unsigned int     ncell,
    __global float2* pos, 
    __global float2* force,
    __global int2*   bounds,
    __constant int*  neighCells
){
    __local  float2 PIs[16];
    __local  float2 PJs[256];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    // --- copy content of neighboring cells to local memory
    
    int njs = 0;
    int  icell = iG/nL; 
    for(int ineigh=0; ineigh<9; ineigh++ ){
        int  ineigh = neighCells[ineigh];  
        int2 IB     = bounds[icell];
        for(int i0=0; i0<IB.y; i0+=nL){
            int i   = i0 + iL; 
            PIs[njs+i] = pos[IB.x+i];
        }
        njs += IB.y;
    }

    int nis = bounds[neighCells[0]].y;
    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // --- compute force using only local memory
    
    float2 f = (float2) (0.0f, 0.0f);
    for (int i0=0; i0<nis; i0 += nL ){
        float2 p = PIs[i0+iL];
        for (int j=0; j<njs; j++ ){
            f += pair_force( p, PJs[j] );
        }
    } 
    force[iG] = f;
}







