// Normalize a float2 vector
float2 normalize2(float2 v) {
    float len = length(v);
    return (len > 1e-10f) ? v / len : (float2)(0.0f, 0.0f);
}

// One iteration of Jacobi solver for position-based dynamics
__kernel void updateJacobi_neighs(
    const int n,              // number of particles
    const int max_neighs,     // maximum number of neighbors per particle
    __global float4* ps,      // particle positions and masses (x,y,z,mass)
    __global float4* ps_out,  // output positions
    __global int* neighs,     // neighbor indices
    __global float* kngs,     // stiffness coefficients
    __global float* l0s,      // rest lengths
    const float inv_dt2,      // 1/dt^2
    const float Rd            // relaxation distance
) {
    int i = get_global_id(0);
    if (i >= n) return;
    
    // Load current position
    float4 pi = ps[i];
    float2 pos_i = (float2)(pi.x, pi.y);
    float mass_i = pi.w;
    
    // Initialize right-hand side with inertial term
    float2 bi = pos_i * (mass_i * inv_dt2);
    
    // Sum over neighbors
    for (int k = 0; k < max_neighs; k++) {
        int j = neighs[i * max_neighs + k];
        if (j < 0) break;  // No more neighbors
        
        float kng = kngs[i * max_neighs + k];
        float l0 = l0s[i * max_neighs + k];
        
        // Get neighbor position
        float4 pj = ps[j];
        float2 pos_j = (float2)(pj.x, pj.y);
        
        // Get direction vector
        float2 dir_ij = pos_j - pos_i;
        dir_ij = normalize2(dir_ij);
        
        // Add contribution from this neighbor
        bi += kng * (pos_j + dir_ij * l0);
    }
    
    // Update position
    float2 pos_new = bi / (mass_i * inv_dt2 + 1.0f);
    ps_out[i] = (float4)(pos_new.x, pos_new.y, 0.0f, mass_i);
}

/// this function updates neighs, kngs, l0s based on collisions found considering current positions of points, it use groups to speed up
/// NOTE: we should compare speed of updateBondsAndCollison_neighs which use pre-calculated collision neighbors vs updateCollisonNeighbors_groups search them on the fly
//    * updateBondsAndCollison_neighs is simpler therefore faster
//    * updateCollisonNeighbors_groups however does not need separate execution of  updateCollisonNeighbors, and storing them in global memory
//        * since updateCollisonNeighbors is less well parallelized, it may be faster to use  updateCollisonNeighbors_groups with on-the-fly search of colision neighbors
__kernel void updateCollisonNeighbors(
    int ngroup,                      // number of groups
    int npoint,                      // number of points
    int nmax_neigh,                  // max number of neighbors
    __global const int4*    granges, // [ngroup] ranges of inds[] for each group {i0,n0,nng}, i0 start of group, n0 number of atoms in the group, nng number of neighbor candidates (r<rc)
    __global const int*     inds,    // [ngroup] ranges of atoms which can corresponding to each group {i0,n}
    __global const float4*  ps,      // [npoint] x,y,z,mass
    __global const float *  Rs,      // [npoint] collision radius
    __global int*     neighs,        // [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global float*   kngs,          // [npoint,nmax_neigh] stiffnesses for bonds to neighbors
    __global float*   l0s,           // [npoint,nmax_neigh] rest lengths for bonds to neighbors
    const float Kcol,                // collision stiffness (TODO: we should define this as per-particle, but then we need to define some mixing rule)
    const float Rp                   // extened potential collision radius, we detect poential neighbors if they are r<(rc+Rd) because during Jacobi iteration they may potential get into collision ( safety margin )
){

    __local float4 ps_loc[32];  // local memory for position of atoms j to find neighbors

    const int ig = get_group_id(0);
    const int iL = get_local_id(0);
    const int nL = get_local_size(0); // we assume nL is 32
    
    const int4 g = granges[ig];
    
    int i = iL;
    while( i < g.y ){
        const int    ia = inds[g.x+i];
        const float3 pi = ps[ia].xyz;
        const float  ri = Rs[ia] + Rp; // extended potential collision radius of particle i

        // count filled neighbors for this point  // TODO: if we have nneigh array (stroing ni) this would be much simpler
        int        ni = 0; 
        const int in0 = ia*nmax_neigh;
        for( int j=0; j<nmax_neigh; j++ ){
            if(neighs[in0+j] < 0) break;
            ni++;
        }

        for( int j0=0; j0<g.z; j0+=nL ){

            { // load j-atom to local memory
                int ja = g.x+j0+iL;
                ps_loc[iL]   = ps[ja];
                ps_loc[iL].w = Rs[ja];
                barrier(CLK_LOCAL_MEM_FENCE);
            }
            for( int jl=0; jl<nL; jl++ ){
                const int j = j0 + jl;
                if( j >= g.z ) break;   // 
                const float3 pj  = ps_loc[jl].xyz;
                const float3 dij = pi.xyz - pj;
                const float  r2  = dot(dij,dij);
                const float  rc  = ri + ps_loc[jl].w; // ri+rj (+Rp) - sum of (potential) collision radius of both particles
                if( r2>(rc*rc) ) continue;

                const int  ja = inds[g.x+j];
                const int ing = in0+ni;
                neighs[ing]   = ja;
                kngs  [ing]   = Kcol;    // TODO: we should define some mixing rule rather than using constant stiffness
                l0s   [ing]   = rc-Rp;   // we store true collision radius ( ri+rj, ommiting Rp)
                ni++;

            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }
}