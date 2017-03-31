
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
    unsigned int      ncell,
    __global  float2* pos, 
    __global  float2* force,
    __global    int2* bounds,
    __global     int* debug,
    __constant   int* neighCells
){
    __local  float2 Ps[256];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    // --- copy content of neighboring cells to local memory
    
    int  icell = iG/nL; 
    int  nis   = bounds[icell].y;
    if( nis==0 ){ debug[icell] = -1; return; }
    int  njs   = 0;
    for(int ineigh=0; ineigh<9; ineigh++ ){
        int  idcell = neighCells[ineigh];  
        int2 bound  = bounds[icell+idcell];
        for(int i0=0; i0<bound.y; i0+=nL){
            int i     = i0 + iL; 
            if(i>=bound.y) break;  // should not matter; question is if it help performance ?
            Ps[njs+i] = pos[bound.x+i];
        }
        njs += bound.y;
    }
    
    if(iL==0) debug[icell] = njs;

    barrier(CLK_LOCAL_MEM_FENCE);
    
    // --- compute force using only local memory
    
    int i0cell = bounds[icell].x;
    for (int i0=0; i0<nis; i0 += nL ){
        int i = i0+iL;
        if (i>=nis) break;
        float2 pi = Ps[i];
        float2 f  = (float2) (0.0f, 0.0f);
        for (int j=0; j<njs; j++ ){
            f += pair_force( pi, Ps[j] );
        }
        force[i0cell+i] = f;
        //force[i0cell+i] = (float2)(151515.0,888888.0);
    } 
    
}







