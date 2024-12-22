// def jacobi_iteration_sparse(x, b, neighs, kngs, Aii ):
//     n = len(x)
//     x_out = np.zeros_like(x)
//     r     = np.zeros_like(x)
//     for i in range(n):
//         sum_j  = 0  # RHS term
//         ngsi = neighs[i]
//         ksi  = kngs[i] 
//         ni = len(ngsi)
//         for jj in range(ni):
//             j      = ngsi[jj]
//             k      = ksi[jj]
//             sum_j += k * x[j]                 # Off-diagonal contribution
//         x_out[i] =  (b[i] + sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
//         r[i]     = b[i] + sum_j - Aii[i]*x[i] # Residual r = b - Ax ;  Ax = Aii * x[i] + sum_(j!=i){ Aij * x[j] }
//     return x_out, r

__kernel void jacobi_iteration_sparse(
    __global const float4* x,        // [npoint,4] solution candiate for position + mass {x,y,z,mass} 
    __global const float4* b,        // [npoint,4] RHS vector b {x,y,z, ? }
    __global const int*   neighs,    // [npoint,nneigh_max] Neighbor indices, padded with -1
    __global const float* kngs,      // [npoint,nneigh_max] Non-zero matrix elements (stiffness)
    __global const float* Aii,       // [npoint] Diagonal elements of the matrix A
    __global       float4* x_out,    // [npoint,4] Output solution vector {x,y,z,mass}
    __global       float4* r,        // [npoint,4] Output residual vector {x,y,z,?}
    const int nmax_neigh,            // Maximum number of neighbors (padding size)
    const int n                      // Number of equations
){
    int iG = get_global_id(0);
    if(iG >= n) return;
    
    // Calculate sum of off-diagonal terms
    float3 sum_j   = (float3){0.0f, 0.0f, 0.0f};
    float4 xi      = x[iG];
    int row_offset = iG * nmax_neigh;
    
    // Loop through neighbors
    for(int jj = 0; jj < nmax_neigh; jj++){
        int j = neighs[row_offset + jj];
        if(j < 0) break;  // Stop at sentinel value (-1)
        float k = kngs[row_offset + jj];
        sum_j += k * x[j].xyz;  // Off-diagonal contribution
    }
    
    float3 bi = b[iG].xyz;
    x_out[iG] = (float4){  ( bi - sum_j) / Aii[iG],         xi.w };   // Calculate new solution
    r[iG]     = (float4){  bi + sum_j - Aii[iG] * xi.xyz,   0.0f };   // Calculate residual
}


















// =========== Optimized kernel using local memory and group-wise processing ===========
//    Not sure if this is faster than the original kernel

// Optimized Jacobi iteration kernel using local memory and group-wise processing
__kernel void jacobi_iteration_sparse_local(
    __global const float4* x,              // [npoint,4] solution candidate for position + mass {x,y,z,mass}
    __global const float4* b,              // [npoint,4] RHS vector b {x,y,z, ?}
    __global const int2*   group2params,   // [ngroup,2] Range of parameters for each group {iGroup0, nGroupNeighs}
    __global const int2*   point2params,   // [npoint,2] Range of parameters for each point {iPoint0, nPointNeighs}
    __global const int2*   group2points,   // Points in each group
    __global const int*    group_points,   // Points in each group
    __global const int*    group2neighs,   // Folded neighbor indices for each group
    __global const float*  group_ngs,      // Folded stiffness values for each group
    __global const float*  group_kngs,     // Folded stiffness values for each group
    __global const float*  Aii,            // Diagonal elements
    __global float4* x_out,                // Output solution vector
    __global float4* r,                    // Output residual vector
    const int max_nG,             // Maximum number of points in a group
    const int max_group_neighs            // Maximum number of total neighbors for a group
){
    // Get work-group and local IDs
    int iGroup   = get_group_id(0);
    int iL       = get_local_id(0);
    int nL       = get_local_size(0); 
    int iG       = get_global_size(0);
    
    __local  float4* loc_xs;         // packed array of all the points which are neighbors of any of the points in this group
    __local  int*    loc_ngs;        // Local cache for neighbor indice
    __local  float*  loc_kngs;       // Local cache for stiffness values

    // load points to local memory
    const int2 gpoint = group2points[iGroup];
    for(int i = iL; i < gpoint.y; i += nL){
        int ii     = group_points[gpoint.x + i];
        loc_xs[i]  = x[ii];
    }

    // load neighbor params into local memory
    const int2 gpar = group2params[iGroup];
    for(int i = iL; i < gpar.y; i += nL){
        if(i < gpar.y){
            int ii = i + gpar.x; 
            loc_ngs [i] = group_ngs [ ii ];
            loc_kngs[i] = group_kngs[ ii ];
        }
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);

    int2 lrange = point2params[iG];

    float4 xi    = x[iG];            // Current point position + mass
    float4 bi    = b[iG];            // RHS vector for current point
    float3 sum_j = (float3)(0.0f);   // Accumulator for neighbor contributions
    
    // Process points assigned to this work-item
    for(int i = iL; i < lrange.y; i += nL ){
        int il  = i + lrange.x;
        int j   = loc_ngs [il];
        float k = loc_kngs[il];
        sum_j += k * loc_xs[j].xyz;  // Only use xyz components
    }

    const float aii = Aii[iG];
    x_out[iG] = (float4)( (bi.xyz + sum_j) / aii,          xi.w );
    r    [iG] = (float4)(  bi.xyz + sum_j  - aii * xi.xyz, 0.0f );
    
}