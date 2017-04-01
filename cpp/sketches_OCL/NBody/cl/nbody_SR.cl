
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

/*
int loadCell( int2 bound, int iL, int nL, __global float2 * pos, __local float2 * Ps ){
    for(int i0=0; i0<bound.y; i0+=nL){
        int i     = i0 + iL; 
        if(i>=bound.y) break;              // should not matter; question is if it help performance ?
        Ps[i] = pos[bound.x+i];
    }
    return bound.y;
    
    njs+=loadCell( bound0            , iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell   -1], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell   +1], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell-nx-1], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell-nx  ], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell-nx+1], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell+nx-1], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell+nx  ], iL, nL, pos, Ps+njs );
    njs+=loadCell( bounds[icell+nx+1], iL, nL, pos, Ps+njs );
}
*/

/*
 Performance test [Mticks]:
     n     nx   ny   ncell  nloc   cell_sz  CPU       Full     Just-Load   Check-cell    Empty    
   4096  32  32   1024       16     4.0              1.07348       0.51388    0.21731    0.46654 
   4096  24  24    576       16     6.0              1.37183       0.43338    0.45145    0.467181 
   4096  24  24    576        8     6.0              1.59139   
   4096  24  24    576       32     6.0    3.27869   1.09454  
   8192  64  64   4096       32     6.0    6.24658   1.60579       
  16384  64  64   4096       32     6.0   12.14813   2.69230  
  32768  64  64   4096       32     6.0   23.73025   4.98328 
*/


__kernel void force_Tiled(
    unsigned int      nx,
    __global  float2* pos, 
    __global  float2* force,
    __global    int2* bounds
){
    __local  float2 Ps[256];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
 
    // --- copy content of neighboring cells to local memory
    
    int icell = iG/nL; 
    int2 bound0 = bounds[icell];
    
    if( bound0.y==0 )return;
    
    #define LOAD_CELL( BOUND ){ \
        int2 bound = BOUND;                  \
        for(int i0=0; i0<bound.y; i0+=nL){   \
            int i     = i0 + iL;             \
            if(i>=bound.y) break;            \
            Ps[njs+i] = pos[bound.x+i];      \
        }                                    \
        njs+=bound.y;          }              
    int njs = 0; 
    LOAD_CELL( bound0             );   
    LOAD_CELL( bounds[icell   -1] );    
    LOAD_CELL( bounds[icell   +1] );    
    LOAD_CELL( bounds[icell-nx-1] );   
    LOAD_CELL( bounds[icell-nx  ] );   
    LOAD_CELL( bounds[icell-nx+1] );    
    LOAD_CELL( bounds[icell+nx-1] );   
    LOAD_CELL( bounds[icell+nx  ] );   
    LOAD_CELL( bounds[icell+nx+1] ); 
    #undef  LOAD_CELL
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // --- compute force using only local memory
    for (int i0=0; i0<bound0.y; i0 += nL ){
        int i = i0+iL;
        if (i>=bound0.y) break;
        float2 pi = Ps[i];
        float2 f  = (float2) (0.0f, 0.0f);
        for (int j=0; j<njs; j++ ){
            f += pair_force( pi, Ps[j] );
        }
        force[bound0.x+i] = f;
    }    
}

__kernel void force_Tiled_Sorted(
    unsigned int      nx,
    __global  float2* pos, 
    __global  float2* force,
    __global    int2* bounds,
    __global    int * order
){
    __local  float2 Ps[256];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
 
    // --- copy content of neighboring cells to local memory
    
    int icell    = iG/nL;
        icell    = order [icell]; // may be slow?   should we do if(iL==0) icell 
    int2 bound0  = bounds[icell];
        
    #define LOAD_CELL( BOUND ){ \
        int2 bound = BOUND;                  \
        for(int i0=0; i0<bound.y; i0+=nL){   \
            int i     = i0 + iL;             \
            if(i>=bound.y) break;            \
            Ps[njs+i] = pos[bound.x+i];      \
        }                                    \
        njs+=bound.y;          }              
    int njs = 0; 
    LOAD_CELL( bound0             );   
    LOAD_CELL( bounds[icell   -1] );       // can be slow ... maybe we should read  bounds[icell-1] to local memory ?
    LOAD_CELL( bounds[icell   +1] );    
    LOAD_CELL( bounds[icell-nx-1] );   
    LOAD_CELL( bounds[icell-nx  ] );   
    LOAD_CELL( bounds[icell-nx+1] );    
    LOAD_CELL( bounds[icell+nx-1] );   
    LOAD_CELL( bounds[icell+nx  ] );   
    LOAD_CELL( bounds[icell+nx+1] ); 
    #undef  LOAD_CELL
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // --- compute force using only local memory
    for (int i0=0; i0<bound0.y; i0 += nL ){
        int i = i0+iL;
        if (i>=bound0.y) break;
        float2 pi = Ps[i];
        float2 f  = (float2) (0.0f, 0.0f);
        for (int j=0; j<njs; j++ ){
            f += pair_force( pi, Ps[j] );
        }
        force[bound0.x+i] = f;
    }    
    
}

/*
__kernel void force_Inters(
    unsigned             int    ninters,
             __global float2*   pos, 
    volatile __global  float2*  force,
             __global    int4*  inters,
){
    __local  float2 pjs[256];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
 
    if(iG<(ninters*nL)) return;
    
    int4 inds = inters[0];
 
    // probably it cannot be done
    // http://stackoverflow.com/questions/12535158/how-to-correclty-sum-results-from-local-to-global-memory-in-opencl
    // https://streamcomputing.eu/blog/2016-02-09/atomic-operations-for-floats-in-opencl-improved/
    for (int i0=0; i0<inds.y; i0 += nL ){
        int i = i0+iL;
        if (i>=inds.y) break;
        float2 pi = pos[inds.x+i];
        float2 f  = (float2) (0.0f, 0.0f);
        for (int j=0; j<inds.w; j++ ){
            f += pair_force( pi, pjs[j] );
        }
        force[inds.x+i] = f;
    }    
}
*/







