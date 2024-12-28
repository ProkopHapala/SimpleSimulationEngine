

__kernel void updateBondsAndCollison_v1( 
    int npoint,                     // number of points
    int nmax_neigh,                 // max number of neighbors
    __global const float4*  ps,     // [npoint] x,y,z,mass
    __global       float4*  ps_out, // [npoint] x,y,z,mass
    //__global const float* Aii,    // [npoint] diagonal of the Projective-Dynamics matrix of indexes
    __global const int*     neighs, // [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float*   kngs,   // [npoint,nmax_neigh] stiffnesses for bonds to neighbors
    __global const float*   l0s,    // [npoint,nmax_neigh] rest lengths for bonds to neighbors
    float inv_dt2                   // 1/dt^2, controls scale of inertial term
){

    const int iG     =  get_global_id(0);
    const int j0     = iG * nmax_neigh;
    const float4 pi  = ps[iG];                       
     
    // inertial term
    const float mi = pi.w*inv_dt2;
    float3 pi_new  = pi.xyz * mi;
    float  Aii     = mi;

    for( int j = 0; j < npoint; j++ ){
        int jG = neighs[j];
        if( jG == -1 ) continue;
        const int jng = j0 + j;
        const float  k   = kngs[jng];
        const float  l0  = l0s [jng];
        const float3 pj  = ps[jG].xyz;
        const float3 dij = pi.xyz - pj;
        const float  l2  = dot(dij,dij);

        if(l0<0){ // collision - if not in repulsion we skip this neighbor
            if( l2<(l0*l0) ) continue;
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
    // where bi = \sum_j K_{ij} d_{ij}    (+ inertial term M_i/dt^2 pi )
    // annd A_{ii} = \sum_j K_{ij}        (+ inertial term M_i/dt^2 )

    pi_new /= Aii;

    ps_out[iG] = (float4)(pi_new, pi.w);

}