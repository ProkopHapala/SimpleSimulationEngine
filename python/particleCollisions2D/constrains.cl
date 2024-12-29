

__kernel void updateJacobi_neighs( 
    int npoint,                     // number of points
    int nmax_neigh,                 // max number of neighbors
    __global const float4*  ps,     // [npoint] x,y,z,mass
    __global       float4*  ps_out, // [npoint] x,y,z,mass
    //__global const float* Aii,    // [npoint] diagonal of the Projective-Dynamics matrix of indexes
    __global const int*     neighs, // [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float*   kngs,   // [npoint,nmax_neigh] stiffnesses for bonds to neighbors
    __global const float*   l0s,    // [npoint,nmax_neigh] rest lengths for bonds to neighbors
    float inv_dt2,                  // 1/dt^2, controls scale of inertial term
    const float Rd                  // reduced collision radius, we detect only hard collision if r<(r-Rd) (i.e. if we are deep in penetration) to prevent instability, soft collison should be solved by non-covalent force-field
){

    const int iG     =  get_global_id(0);
    const int j0     = iG * nmax_neigh;
    const float4 pi  = ps[iG];                       
     
    // inertial term
    const float mi = pi.w*inv_dt2;
    float3 pi_new  = pi.xyz * mi;
    float  Aii     = mi;

    for( int j = 0; j < nmax_neigh; j++ ){
        int jG = neighs[j];
        if( jG == -1 ) break;
        const int jng = j0 + j;
        const float  k   = kngs[jng];
        float        l0  = l0s [jng];
        const float3 pj  = ps[jG].xyz;
        const float3 dij = pi.xyz - pj;
        const float  l2  = dot(dij,dij);

        if(l0<0){ // collision - if not in repulsion we skip this neighbor
            const float lc = l0 - Rd;  //  we cut-off collision if we are noit deep enough in penetration 
            if( l2<(lc*lc) ) continue;
        }
        const float l    = sqrt(l2);
        const float3 pij = pj + dij * l0/l;  // p'_{ij} is ideal position of particle i which satisfy the constraint between i and j
        // ToDo: we should make sure pij originating from collision does not reach beyond l0
        pi_new += pij * k;   // 
        Aii    += k;         // diagonal of the Projective-Dynamics matrix, \sum_j k_{ij}       (+ inertial term M_i/dt^2)
    }

    // NOTE 1: pi_new is basically weighted average of predicted positions of particles and due to all the constraints weighted by stiffness (and due to intertial term M_i/dt^2)
    // pi_new   = \sum_j pij * k_{ij} / \sum_j k_{ij}
    // NOTE 2: at the same time, pi_new can be seen as one step of Jacobi iteration for Projective-Dynamics equation Ap=b
    // pi_new   = (bi + sum_j) / A_{ii} 
    // where  bi     = \sum_j K_{ij} d_{ij}    (+ inertial term M_i/dt^2 pi )
    // and    A_{ii} = \sum_j K_{ij}           (+ inertial term M_i/dt^2    )

    pi_new /= Aii;

    ps_out[iG] = (float4)(pi_new, pi.w);

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