
int rand(int seed){
    //int const a = 16807; //ie 7**5
    //int const m = 2147483647; //ie 2**31-1
    return (long(seed* 16807))%2147483647;
    return seed;
}

static float noise3D(float x, float y, float z) {
    float ptr = 0.0f;
    return fract(sin(x*112.9898f + y*179.233f + z*237.212f) * 43758.5453f, &ptr);
}

__kernel void scatterFromPoint( 
    __global float4* triangles, // triangles of scattering material; each float4 one vertex including thickness
    __global float4* materials, // material for each triangle 
    __global float4* targets,   // positions and radius of spherical targets
    float*           count,     // count of particles which hit the target
    int ntris,
    int nbounce,
    float4 source,
){
    __local float4 LTRIS[32];
    
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    float3 rd = randomDirection(iG);
    for (int i=0;i<nbounce;i++){
        
        for (int i0=0; i0<nAtoms; i0+= nL ){
            int i = i0 + iL;
            if(i>=nAtoms) break;

            LATOMS[iL] = atoms[i];
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int j=0; j<nL; j++){
                fe += getCoulomb( LATOMS[j], pos );
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }
}
