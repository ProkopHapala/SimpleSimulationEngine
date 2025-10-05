// def jacobi_iteration_sparse(x, b, neighs, kngs, Aii ):
//     """One iteration of Jacobi method using sparse operations"""
//     n = len(x)
//     x_out = np.zeros_like(x)
//     r     = np.zeros_like(x)
//     for i in range(n):
//         sum_j  = 0  # RHS term
//         ngsi = neighs[i]
//         ksi  = kngs[i] 
//         ni = len(ngsi)
//         #print(f"\nCPU Point {i}:")
//         #print(f"  Initial sum_j: {sum_j}")
//         #print(f"  Neighbors: {ngsi[:ni]}")
//         #print(f"  Stiffness: {ksi[:ni]}")
//         for jj in range(ni):
//             j      = ngsi[jj]
//             k      = ksi[jj]
//             sum_j += k * x[j]   # Off-diagonal contribution
//             #print(f"    j={j}, k={k}, x[j]={x[j]}, sum_j={sum_j}")
//         x_out[i] =  (b[i] + sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
//         r[i]     = b[i] - (-sum_j + Aii[i]*x[i]) # Residual r = b - Ax ;  Ax = Aii * x[i] + sum_(j!=i){ Aij * x[j] }

// TODO(VBD): Vertex Block Descent kernels must be launched per color to preserve Gauss-Seidel ordering.
// Each launch will iterate over pre-packed workgroups of size 32, where the host provides
// (a) up to 32 vertex indices for the block and (b) the full set of their neighbor indices
// (<=128) so that no entries from the `neighs` table are missing within shared local memory.
// Python driver is responsible for building these chunk buffers before dispatch.
//         #print(f"  Final: b[i]={b[i]}, Aii[i]={Aii[i]}, x_out[i]={x_out[i]}, r[i]={r[i]}")
//     return x_out, r


#define VBD_WORKGROUP_SIZE 32
#define VBD_NEIGHBOR_LIMIT 128

inline void init_hessian(float* H, float diag){
    H[0]=H[4]=H[8]=diag;
    H[1]=H[2]=H[3]=H[5]=H[6]=H[7]=0.0f;
}

inline void accum_vertex_hessian(float* H, float3 dir, float k, float rest, float len){
    float coeff_iso, coeff_dir;
    if (len > 1e-7f){
        float inv_len = 1.0f / len;
        coeff_iso = k * (1.0f - rest * inv_len);
        coeff_dir = k * (rest * inv_len);
    }else{
        coeff_iso = k;
        coeff_dir = k;
        dir = (float3)(0.0f);
    }

    float3 row0 = (float3)(dir.x * dir.x, dir.x * dir.y, dir.x * dir.z);
    float3 row1 = (float3)(dir.y * dir.x, dir.y * dir.y, dir.y * dir.z);
    float3 row2 = (float3)(dir.z * dir.x, dir.z * dir.y, dir.z * dir.z);

    H[0] += coeff_iso + coeff_dir * row0.x;
    H[1] += coeff_dir * row0.y;
    H[2] += coeff_dir * row0.z;
    H[3] += coeff_dir * row1.x;
    H[4] += coeff_iso + coeff_dir * row1.y;
    H[5] += coeff_dir * row1.z;
    H[6] += coeff_dir * row2.x;
    H[7] += coeff_dir * row2.y;
    H[8] += coeff_iso + coeff_dir * row2.z;
}

inline float3 solve3x3(const float* H, float3 grad, const float det_eps){
    const float c00 = H[4]*H[8] - H[5]*H[7];
    const float c01 = H[2]*H[7] - H[1]*H[8];
    const float c02 = H[1]*H[5] - H[2]*H[4];
    const float c10 = H[5]*H[6] - H[3]*H[8];
    const float c11 = H[0]*H[8] - H[2]*H[6];
    const float c12 = H[2]*H[3] - H[0]*H[5];
    const float c20 = H[3]*H[7] - H[4]*H[6];
    const float c21 = H[1]*H[6] - H[0]*H[7];
    const float c22 = H[0]*H[4] - H[1]*H[3];
    const float det = H[0]*c00 + H[1]*c10 + H[2]*c20;
    if (fabs(det) < det_eps) return (float3)(0.0f, 0.0f, 0.0f);
    const float inv_det = 1.0f / det;
    float3 dx;
    dx.x = (c00*grad.x + c01*grad.y + c02*grad.z) * inv_det;
    dx.y = (c10*grad.x + c11*grad.y + c12*grad.z) * inv_det;
    dx.z = (c20*grad.x + c21*grad.y + c22*grad.z) * inv_det;
    return dx;
}

__kernel void jacobi_iteration_sparse(
    __global const float4* x,        // 1 : [npoint,4] solution candiate for position + mass {x,y,z,mass} 
    __global const float4* b,        // 2 : [npoint,4] RHS vector b {x,y,z, ? }
    __global const int*   neighs,    // 3 : [npoint,nneigh_max] Neighbor indices, padded with -1
    __global const float* kngs,      // 4 : [npoint,nneigh_max] Non-zero matrix elements (stiffness)
    __global const float* Aii,       // 5 : [npoint] Diagonal elements of the matrix A
    __global       float4* x_out,    // 6 : [npoint,4] Output solution vector {x,y,z,mass}
    __global       float4* r,        // 7 : [npoint,4] Output residual vector {x,y,z,?}
    const int nmax_neigh,            // 8 : Maximum number of neighbors (padding size)
    const int n                      // 9 : Number of equations
){
    const int iG = get_global_id(0);
    if(iG >= n) return;

    // if(iG==0){  
    //     //for(int i=0; i<n; i++){  printf("GPU i: %i Aii[i]: %f b[i]: %f \n", i, Aii[i], b[i].x ); }
    //     for(int iG=0; iG<n; iG++){ 
    //         // Calculate sum of off-diagonal terms
    //         float3 sum_j   = (float3){0.0f, 0.0f, 0.0f};
    //         float4 xi      = x[iG];
    //         int j0 = iG * nmax_neigh;
    //         // printf("\nGPU Point %d:\n", iG);
    //         // printf("  Initial sum_j: (%f, %f, %f)\n", sum_j.x, sum_j.y, sum_j.z);
    //         // Loop through neighbors
    //         for(int jj = 0; jj < nmax_neigh; jj++){
    //             int j = neighs[j0 + jj];
    //             if(j < 0) break;  // Stop at sentinel value (-1)
    //             float k = kngs[j0 + jj];
    //             float4 xj = x[j];
    //             sum_j += k * xj.xyz;  // Off-diagonal contribution
    //             //printf("    j=%d, k=%f, x[j]=(%f, %f, %f), sum_j=(%f, %f, %f)\n",    j, k, xj.x, xj.y, xj.z, sum_j.x, sum_j.y, sum_j.z);
    //         }
    //         float3 bi = b[iG].xyz;
    //         float3 x_out_i = (bi + sum_j) / Aii[iG];  // Calculate new solution
    //         float3 r_i     = bi +  sum_j - Aii[iG] * xi.xyz;  // Calculate residual
    //         //printf("GPU i: %i Aii[i]: %f b[i]: %f sum_j: %f  x_out[i]: %f r[i]: %f\n", iG, Aii[iG], b[iG].x, sum_j.x,  x_out_i.x, r_i.x);
    //         // float3 bi = b[iG].xyz;
    //         // x_out[iG] = (float4){ (bi + sum_j) / Aii[iG], xi.w };  // Calculate new solution
    //         // r[iG]     = (float4){ bi +  sum_j - Aii[iG] * xi.xyz, 0.0f };  // Calculate residual
    //     }
    // }
    
    
    const float4 xi     = x[iG];
    const int    j0     = iG * nmax_neigh;  // offset of the first neighbor for this point
    float3       sum_j  = (float3){0.0f, 0.0f, 0.0f};
    
    // printf("\nGPU Point %d:\n", iG);
    // printf("  Initial sum_j: (%f, %f, %f)\n", sum_j.x, sum_j.y, sum_j.z);
    
    for(int jj = 0; jj < nmax_neigh; jj++){
        const int j      = neighs[j0 + jj];
        if(j < 0) break;  // Stop at sentinel value (-1)
        const float k    = kngs[j0 + jj];
        sum_j           += k * x[j].xyz;
        //printf("    j=%d, k=%f, x[j]=(%f, %f, %f), sum_j=(%f, %f, %f)\n",    j, k, xj.x, xj.y, xj.z, sum_j.x, sum_j.y, sum_j.z);
    }
    
    const float3 bi = b[iG].xyz;
    x_out[iG] = (float4){ (bi + sum_j) / Aii[iG], xi.w };           // Calculate new solution
    r[iG]     = (float4){  bi + sum_j -  Aii[iG] * xi.xyz, 0.0f };  // Calculate residual

    //x_out[i] =  (b[i] + sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
    //r[i]     = b[i] - (-sum_j + Aii[i]*x[i]) # Residual r = b - Ax ;  Ax = Aii * x[i] + sum_(j!=i){ Aij * x[j] }
    //printf("  Final: b[i]=(%f, %f, %f), Aii[i]=%f\n", bi.x, bi.y, bi.z, Aii[iG]);
    //printf("  x_out[i]=(%f, %f, %f), r[i]=(%f, %f, %f)\n",    x_out[iG].x, x_out[iG].y, x_out[iG].z,  r[iG].x, r[iG].y, r[iG].z);
    //printf( "  Final: b[i]=%f, Aii[i]=%f, x_out[i]=%f, r[i]=%f \n", b[iG].x, Aii[iG], x_out[iG].x, r[iG].x );
}

__kernel void gauss_seidel_iteration_colored(
    __global       float4* x,           // 1 : [npoint,4] solution candidate for position + mass {x,y,z,mass}, both input and output as Gauss-Seidel operates in-place
    __global const float4* b,           // 2 : [npoint,4] RHS vector b {x,y,z, ?}
    __global const int*   neighs,       // 3 : [npoint,nneigh_max] Neighbor indices, padded with -1
    __global const float* kngs,         // 4 : [npoint,nneigh_max] Non-zero matrix elements (stiffness)
    __global const float* Aii,          // 5 : [npoint] Diagonal elements of the matrix A
    __global const int*   color_group,  // 6 : [npoint] Node indices for the current color group
    __global       float4* r,           // 7 : [npoint,4] Output residual vector {x,y,z,?}
    const int nmax_neigh,               // 8 : Maximum number of neighbors (padding size)
    const int n,                        // 9 : Total number of nodes
    const int group_start,              // 10 : Starting index in the color group buffer
    const int group_count               // 11 : Number of nodes in this color group
){
    const int iG = get_global_id(0);
    if(iG >= group_count) return;
    const int i = color_group[group_start + iG];   // Get the actual node index from the color group

    const float4 xi    = x[i];
    const int    j0    = i * nmax_neigh;   // offset of the first neighbor for this point
    float3       sum_j = (float3){0.0f, 0.0f, 0.0f};
    
    for(int jj = 0; jj < nmax_neigh; jj++){
        const int j = neighs[j0 + jj];
        if(j < 0) break;  // Stop at sentinel value (-1)
        const float k = kngs[j0 + jj];
        sum_j += k * x[j].xyz;  // Off-diagonal contribution
    }
    const float3 bi = b[i].xyz;  // Fixed: Use i instead of iG
    x[i] = (float4){ (bi + sum_j) / Aii[i], xi.w };  // Update solution in-place
    r[i] = (float4){ bi + sum_j - Aii[i] * xi.xyz, 0.0f };  // Calculate residual
}



__kernel void jacobi_iteration_sparse_local(
    __global const float4* x,              // 1 : [npoint,4] solution candidate for position + mass {x,y,z,mass}
    __global const float4* b,              // 2 : [npoint,4] RHS vector b {x,y,z, ?}
    __global const int2*   group2params,   // 3 : [ngroup,2] Range of parameters for each group {iGroup0, nGroupNeighs}
    __global const int2*   point2params,   // 4 : [npoint,2] Range of parameters for each point {iPoint0, nPointNeighs}
    __global const int2*   group2points,   // 5 : [ngroup,2] Range of points for each group {iPoint0, nPointNeighs}
    __global const int*    group_points,   // 6 : [ngroup,2] Range of points for each group {iPoint0, nPointNeighs}
    __global const int*    group2neighs,   // 7 : Folded neighbor indices for each group
    __global const float*  group_ngs,      // 8 : Folded stiffness values for each group
    __global const float*  group_kngs,     // 9 : Folded stiffness values for each group
    __global const float*  Aii,            // 10 : Diagonal elements
    __global float4* x_out,                // 11 : Output solution vector
    __global float4* r,                    // 12 : Output residual vector
    const int max_nG,                      // 13 : Maximum number of points in a group
    const int max_group_neighs             // 14 : Maximum number of total neighbors for a group
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
    
    // printf("\nGPU Point %d:\n", iG);
    // printf("  Initial sum_j: (%f, %f, %f)\n", sum_j.x, sum_j.y, sum_j.z);
    
    // Process points assigned to this work-item
    for(int i = iL; i < lrange.y; i += nL ){
        int il  = i + lrange.x;
        int j   = loc_ngs [il];
        float k = loc_kngs[il];
        sum_j += k * loc_xs[j].xyz;  // Only use xyz components
       // printf("    j=%d, k=%f, x[j]=(%f, %f, %f), sum_j=(%f, %f, %f)\n",    j, k, loc_xs[j].x, loc_xs[j].y, loc_xs[j].z, sum_j.x, sum_j.y, sum_j.z);
    }

    const float aii = Aii[iG];
    x_out[iG] = (float4)( (bi.xyz - sum_j) / aii,          xi.w );
    r    [iG] = (float4)(  bi.xyz - (sum_j + aii * xi.xyz), 0.0f );
    
    //printf("  Final: b[i]=(%f, %f, %f), Aii[i]=%f\n", bi.x, bi.y, bi.z, Aii[iG]);
    //printf("  x_out[i]=(%f, %f, %f), r[i]=(%f, %f, %f)\n",     x_out[iG].x, x_out[iG].y, x_out[iG].z,  r[iG].x, r[iG].y, r[iG].z);
}


// =========== Optimized kernel using local memory and group-wise processing ===========
//    Not sure if this is faster than the original kernel
//    TODO(VBD): Vertex Block Descent kernels must be launched per color to preserve Gauss-Seidel ordering.
//    Each launch will iterate over pre-packed workgroups of size 32, where the host provides
//    (a) up to 32 vertex indices for the block and (b) the full set of their neighbor indices
//    (<=128) so that no entries from the `neighs` table are missing within shared local memory.
//    Python driver is responsible for building these chunk buffers before dispatch.

__kernel void vbd_vertex_chunk(
    __global float4* x,                  // 1 : [nverts,4] positions
    __global const float4* y,            // 2 : [nverts,4] rest positions
    __global const int* chunk_vertices,  // 3 : [VBD_WORKGROUP_SIZE] vertex indices
    __global const int* chunk_neighbors, // 4 : [VBD_WORKGROUP_SIZE*MAX_NEIGHBORS] neighbor indices
    __global const int* slot_table,      // 5 : [VBD_WORKGROUP_SIZE*MAX_NEIGHBORS] slot table
    __global const int* neighs,          // 6 : [nverts*MAX_NEIGHBORS] neighbors
    __global const float* kngs,          // 7 : [nverts*MAX_NEIGHBORS] stiffness
    __global const float* rest_lengths,  // 8 : [nverts*MAX_NEIGHBORS] rest length
    const int vcount,                    // 9 : number of vertices
    const int ncount,                    // 10 : number of neighbors
    const int nmax,                      // 11 : maximum number of neighbors
    const float inv_h2,                  // 12 : inverse of h^2
    const float det_eps                  // 13 : determinant epsilon
){
    const int lid = get_local_id(0);
    const int gid = get_global_id(0);
    if (lid >= VBD_WORKGROUP_SIZE || gid >= VBD_WORKGROUP_SIZE) return;

    __local float4 loc_x[VBD_WORKGROUP_SIZE];
    __local float4 loc_y[VBD_WORKGROUP_SIZE];
    __local float4 loc_neighbors[VBD_NEIGHBOR_LIMIT];
    __local int loc_neighbor_index[VBD_WORKGROUP_SIZE * 128];

    if (lid < vcount){
        const int vid = chunk_vertices[lid];
        loc_x[lid] = x[vid];
        loc_y[lid] = y[vid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int idx = lid; idx < ncount; idx += VBD_WORKGROUP_SIZE){
        const int nid = chunk_neighbors[idx];
        loc_neighbors[idx] = x[nid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int idx = lid; idx < vcount * nmax; idx += VBD_WORKGROUP_SIZE){
        loc_neighbor_index[idx] = slot_table[idx];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (lid >= vcount) return;

    const int vid = chunk_vertices[lid];
    float4 xi = loc_x[lid];
    const float4 yi = loc_y[lid];

    const float mass = xi.w;
    const float weight = mass * inv_h2;

    float3 grad = weight * (xi.xyz - yi.xyz);
    float H[9];
    init_hessian(H, weight);

    const int base = vid * nmax;
    const int slot_base = lid * nmax;

    for (int jj = 0; jj < nmax; ++jj){
        const int packed = loc_neighbor_index[slot_base + jj];
        if (packed < 0) break;

        float4 xj = loc_neighbors[packed];
        const float k = kngs[base + jj];
        const float rest = rest_lengths[base + jj];

        float3 d   = xi.xyz - xj.xyz;
        float len  = length(d);
        float3 dir = (float3)(0.0f);
        if (len > 1e-7f) dir = d / len;

        grad += k * (len - rest) * dir;
        accum_vertex_hessian(H, dir, k, rest, len);
    }

    grad += inv_h2 * (xi.xyz - yi.xyz);
    const float3 dx = solve3x3(H, grad, det_eps);
    xi.xyz -= dx;
    x[vid] = xi;
}

__kernel void vbd_vertex_serial(
    __global float4* x,                 // 1 : [nverts,4] positions
    __global const float4* y,           // 2 : [nverts,4] rest positions
    __global const int*   neighs,       // 3 : [nverts*MAX_NEIGHBORS] neighbors
    __global const float* kngs,         // 4 : [nverts*MAX_NEIGHBORS] stiffness
    __global const float* rest_lengths, // 5 : [nverts*MAX_NEIGHBORS] rest length
    const int n_points,                 // 6 : number of vertices
    const int nmax,                     // 7 : maximum number of neighbors
    const float inv_h2,                 // 8 : inverse of h^2
    const float det_eps                 // 9 : determinant epsilon
){
    if (get_global_id(0) != 0) return;

    for (int vid = 0; vid < n_points; ++vid){
        float4 xi = x[vid];
        float4 yi = y[vid];

        const float mass = xi.w;
        const float weight = mass * inv_h2;

        float3 grad = weight * (xi.xyz - yi.xyz);
        float H[9];
        init_hessian(H, weight);

        const int base = vid * nmax;
        for (int jj = 0; jj < nmax; ++jj){
            const int nb = neighs[base + jj];
            if (nb < 0) break;
            const float k = kngs[base + jj];
            const float rest = rest_lengths[base + jj];

            const float3 xj = x[nb].xyz;
            float3 d = xi.xyz - xj;
            float len = length(d);
            const float safe_len = fmax(len, 1e-7f);
            float3 dir = (float3)(0.0f);
            if (len > 1e-7f){
                dir = d / safe_len;
            }

            const float stretch = len - rest;
            grad += k * stretch * dir;

            accum_vertex_hessian(H, dir, k, rest, len);
        }

        // if (vid == 1){
        //     printf("[GPU] grad[1] = (%10.4g, %10.4g, %10.4g)\n", grad.x, grad.y, grad.z);
        //     printf("[GPU] H[1]    = [%10.4g %10.4g %10.4g; %10.4g %10.4g %10.4g; %10.4g %10.4g %10.4g]\n", H[0],H[1],H[2], H[3],H[4],H[5], H[6],H[7],H[8]);
        // }

        const float3 dx = solve3x3(H, grad, det_eps);
        xi.xyz -= dx;
        x[vid] = xi;
    }
}


// =============================================================================
// ===  Displacement ("Diff") Projective Dynamics helpers ======================
// =============================================================================

// Serial Jacobi iteration operating on displacement variables (dp).
__kernel void jacobi_iteration_diff_serial(
    __global const float4* x_in,    // 1 : [npoint,4] current displacement iterate
    __global const float4* bvec,    // 2 : [npoint,4] RHS xyz + diagonal w
    __global const int*    neighs,  // 3 : [npoint,nmax]
    __global const float*  kngs,    // 4 : [npoint,nmax]
    __global       float4* x_out,   // 5 : [npoint,4] output displacement iterate
    const int nmax,                 // 6 : Maximum number of neighbors (padding size)
    const int npoint                // 9 : Number of equations
){
    if(get_global_id(0) != 0) return;
    for(int i=0; i<npoint; i++){
        float3 sum_j = (float3)(0.0f);
        const int j0 = i * nmax;
        for(int jj=0; jj<nmax; jj++){
            const int j = neighs[j0 + jj];
            if(j < 0) break;
            const float k = kngs[j0 + jj];
            sum_j += k * x_in[j].xyz;
        }
        const float4 bi  = bvec[i];
        const float  Aii = bi.w;
        float3 xi = (bi.xyz + sum_j) / Aii;
        x_out[i] = (float4)(xi, x_in[i].w);
    }
}

// Serial momentum mixing stage mirroring updateIterativeMomentumDiff().
__kernel void pd_momentum_mix_serial(
    __global       float4* x_curr,     // 1 : [npoint,4] current displacement iterate (psa)
    __global const float4* x_new,      // 2 : [npoint,4] freshly computed iterate (psb)
    __global       float4* dps_store,  // 3 : [npoint,4] momentum buffer (linsolve_yy)
    const float bmix,                  // 4 : Mixing parameter
    const int   npoint                 // 5 : Number of equations
){
    if(get_global_id(0) != 0) return;
    for(int i=0; i<npoint; i++){
        float3 p_new = x_new[i].xyz;
        float3 mom   = dps_store[i].xyz;
        float3 p_mixed = p_new + mom * bmix;
        float3 delta   = p_mixed - x_curr[i].xyz;
        x_curr   [i] = (float4)(p_mixed, x_curr[i].w);
        dps_store[i] = (float4)(delta,  dps_store[i].w);
    }
}

#define MAX_NEIGHBORS 16

/**
 * @brief kernel jacobi_fly : Iterative Jacobi solver for truss dynamics with non-linear "on-the-fly" updates.
 *
 * Implements the Jacobi relaxation method where the system matrix (A) and right-hand-side (b)
 * are re-evaluated in every iteration based on the current vertex positions. This makes it
 * suitable for large, non-linear deformations.
 *
 * Algorithm:
 *   - Jacobi: All vertex updates in an iteration are based on the state from the *previous* complete iteration.
 *   - On-the-Fly: Forces and matrix diagonals (A_ii) are recomputed from scratch each iteration.
 *
 * Design & Performance:
 *   - The entire iterative process is contained within a single kernel launch to minimize overhead.
 *   - The Jacobi method's requirement for a separate read/write state is handled by reading a snapshot
 *     of global memory into `__local` memory each iteration. This is a trade-off to strictly conserve
 *     `__local` memory at the cost of higher memory bus traffic compared to a double-buffered approach.
 *   - Assumes a single workgroup is used for all vertices.
 */
__kernel void jacobi_fly(
    __global float4*      pos,       // 1 : [nverts,4] positions
    __global const int*   neighs,    // 2 : [nverts*MAX_NEIGHBORS] neighbors
    __global const float* kngs,      // 3 : [nverts*MAX_NEIGHBORS] stiffness
    __global const float* l0ngs,     // 4 : [nverts*MAX_NEIGHBORS] rest length
    const int nverts,                // 5 : number of vertices
    const float dt,                  // 6 : time step
    const int nIters                 // 7 : number of iterations
){
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);
    __local float4 lpos[1024];
    float inv_dt2 = 1.0f / (dt*dt);
    for (int iter=0; iter<nIters; ++iter) {
        for (int i = lid; i < nverts; i += lsize) {lpos[i] = pos[i];} 
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int i=lid; i<nverts; i+=lsize ){
            float3 pi = lpos[i].xyz;
            float  mi = lpos[i].w;
            // inertial contribution: M_i/dt^2 * p_i
            float  Ii = mi*inv_dt2;
            float3 bi = pi*Ii;
            float Aii = Ii;
            for (int jj=0; jj<MAX_NEIGHBORS; ++jj) {
                int idx = i * MAX_NEIGHBORS + jj;
                int j      = neighs[idx];
                if (j < 0) break;
                // load parameters (float)
                float  k   = kngs [idx];
                float  l0  = l0ngs[idx];
                float3 pj  = lpos[j].xyz;
                float3 dij = pi - pj;   
                float len  = length(dij);
                float c    = k * ( l0/len - 1.0f );
                bi  += dij * c;      //  b_i    +=  \sum_j ( K_{ij} d_{ij} )   
                bi  += pj  * k;      //  s_j    +=  \sum_j ( K_{ij} p_j    )
                Aii += k;            //  A_{ii} +=  \sum_j   K_{ij} 
            }
            // Jacobi update
            float invA    = 1.0f / Aii;
            float3 pi_new = bi * invA;
            pos[i] = (float4)(pi_new, mi);  // store to global memory, so local memory is still in original state (cosistent with jacobi update)
        }
        barrier(CLK_LOCAL_MEM_FENCE); // This is now SUFFICIENT
    }
}

/**
 * @brief kernel GS_fly : Iterative Gauss-Seidel solver for truss dynamics with non-linear "on-the-fly" updates.
 *
 * Implements the Gauss-Seidel relaxation method using graph coloring to enable parallelism.
 * The system is re-evaluated each iteration, making it robust for non-linear dynamics.
 * It typically converges faster than the equivalent Jacobi method.
 *
 * Algorithm:
 *   - Gauss-Seidel: Updates use the most recent values from neighboring vertices within the same iteration.
 *   - Parallelism is achieved by partitioning vertices into independent sets (colors).
 *   - On-the-Fly: Forces and matrix diagonals (A_ii) are recomputed each iteration.
 *
 * Design & Performance:
 *   - A single kernel launch performs all iterations, loading vertex data into `__local` memory once.
 *   - A single `__local` memory buffer is used for fast, in-place updates.
 *   - A `barrier` after each color sweep ensures updated values are visible for the next color set.
 *   - This design minimizes global memory traffic for maximum performance.
 */
__kernel void GS_fly(
    __global float4*      pos,           // 1 : [nverts] positions
    __global const int*   neighs,        // 2 : [nverts*MAX_NEIGHBORS] neighbors
    __global const float* kngs,          // 3 : [nverts*MAX_NEIGHBORS] stiffness
    __global const float* l0ngs,         // 4 : [nverts*MAX_NEIGHBORS] rest length
    const int nverts,                    // 5 : number of vertices
    const float dt,                      // 6 : time step
    const int nIters,                    // 7 : number of iterations
    __global const int*   color_list,    // 8 : [num_colors] color list
    __global const int*   color_offsets, // 9 : [num_colors+1] color offsets
    const int             num_colors     // 10 : number of colors
){
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);

    __local float4 lpos[1024]; // adjust to your local mem budget
    for (int i=lid; i<nverts; i+=lsize) { lpos[i] = pos[i]; }
    barrier(CLK_LOCAL_MEM_FENCE);
    const float inv_dt2 = 1.0f / (dt * dt);  // precompute inv_dt2
    for (int iter = 0; iter < nIters; ++iter) {
        // sweep colors
        for (int c = 0; c < num_colors; ++c) {
            int off         = color_offsets[c];
            int off_next    = color_offsets[c+1];
            int color_count = off_next - off;

            // each thread processes multiple entries of this color
            for (int local_idx = lid; local_idx < color_count; local_idx += lsize) {
                int vid   = color_list[off + local_idx];
                float3 pi = lpos[vid].xyz;
                float  mi = lpos[vid].w;

                // inertial contribution 
                float  Aii  = mi * inv_dt2; // A_{ii} =  M_i/dt^2        
                float3 bi   = pi * Aii;     // b_i    =  M_i/dt^2 p'_i

                for (int e=0; e<MAX_NEIGHBORS; ++e) {
                    int idx = vid * MAX_NEIGHBORS + e;
                    int j      = neighs[idx];
                    if (j < 0) break;
                    // load parameters (float)
                    float k    = kngs [idx];
                    float l0   = l0ngs[idx];
                    float3 pj  = lpos[j].xyz;
                    float3 dij = pi - pj;  
                    float len  = length(dij);
                    float c    = k * ( l0/len - 1.0f );
                    bi  += dij * c;      //  b_i    +=  \sum_j ( K_{ij} d_{ij} )   
                    bi  += pj  * k;      //  s_j    +=  \sum_j ( K_{ij} p_j    )
                    Aii += k;            //  A_{ii} +=  \sum_j   K_{ij} 
                } // neighbor loop
                // Gauss-Seidel update
                float invA = 1.0f / Aii;
                float3 pi_new = bi * invA;
                lpos[vid]  = (float4)( pi_new, mi ); // store to local memory so we have fresh values for next iteration (consistent with Gauss-Seidel update)
            } // per-thread local_idx
            barrier(CLK_LOCAL_MEM_FENCE);  // ensure updated color values visible to other threads
        } // colors
    } // iterations
    for (int i=lid; i<nverts; i+=lsize) { pos[i]=lpos[i];} // write back local dp to global memory
}

/**
 * @brief kernel jacobi_truss_diff : Iterative linear Jacobi solver for pre-computed truss dynamics.
 *
 * Solves the linear system `A * dp = b'` for vertex displacements `dp`, where the system matrix
 * diagonal (A_ii) and the right-hand-side (b') are pre-computed and passed in via `bvec`.
 * This version is faster per-iteration but assumes a linearized system around a reference state.
 *
 * Algorithm:
 *   - Jacobi: All vertex updates in an iteration depend only on the results of the previous iteration.
 *   - Linear ("Diff"): Operates on a fixed system; only performs a sum of `k_ij * dp_j`.
 *
 * Design & Performance:
 *   - All iterations are performed within one kernel launch.
 *   - Uses a `__local` memory ping-pong buffer scheme (`ldp_read`, `ldp_write`) to keep all
 *     iterative updates on-chip, maximizing performance by avoiding global memory access within the loop.
 *   - Final results are written back once after all iterations are complete.
 */
__kernel void jacobi_truss_diff(
    __global float4*      dp,      // 1 : [nverts,4] in/out displacements dp = p - p0
    __global const int*   neighs,  // 2 : [nverts,MAX_NEIGHBORS] neighbor indices
    __global const float* kngs,    // 3 : [nverts,MAX_NEIGHBORS] stiffness k_ij for each neighbor
    __global const float4* bvec,   // 4 : [nverts,4] Precomputed b'.xyz and Aii'.w
    const int nverts,              // 5 : Number of vertices
    const int nIters               // 6 : Number of iterations
){
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);

    __local float4 ldp_A[1024]; // Ping-pong buffer A
    __local float4 ldp_B[1024]; // Ping-pong buffer B

    for (int i = lid; i < nverts; i += lsize) { ldp_A[i] = dp[i]; } // Load initial dp once
    barrier(CLK_LOCAL_MEM_FENCE);

    __local float4* ldp_read = ldp_A;
    __local float4* ldp_write = ldp_B;

    for (int iter = 0; iter < nIters; ++iter) {
        for (int i=lid; i<nverts; i+=lsize ){
            float3 sum_j = (float3)(0.0f);
            for (int jj=0; jj<MAX_NEIGHBORS; ++jj) {
                int idx = i * MAX_NEIGHBORS + jj;
                int j   = neighs[idx];
                if (j < 0) break;
                float  k   = kngs[idx];
                float3 dpj = ldp_read[j].xyz;
                sum_j += dpj * k; // sum_j += k_ij * dp_j
            }
            const float4 bi = bvec[i]; // bi.xyz is RHS, bi.w is A_ii
            float3 dpi_new = (bi.xyz + sum_j) * (1.0f / bi.w);
            ldp_write[i] = (float4)(dpi_new, 0.0f); // Store new displacement
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        __local float4* temp = ldp_read; // Swap buffers for next iteration
        ldp_read = ldp_write;
        ldp_write = temp;
    }

    for (int i=lid; i<nverts; i+=lsize) { dp[i] = ldp_read[i]; } // Write final result back to global
}

/**
 * @brief kernel GS_diff : Iterative linear Gauss-Seidel solver for pre-computed truss dynamics.
 *
 * Solves the linear system `A * dp = b'` using the Gauss-Seidel method with graph coloring.
 * The system diagonal (A_ii) and RHS (b') are passed in via `bvec`. It offers faster
 * convergence than the linear Jacobi method for the same problem.
 *
 * Algorithm:
 *   - Gauss-Seidel: Vertex updates use the newest available neighbor displacements from the same iteration.
 *   - Graph coloring enables parallel updates of independent vertex sets.
 *   - Linear ("Diff"): Performs a simple weighted sum of neighbor displacements.
 *
 * Design & Performance:
 *   - A single kernel launch contains the entire iterative loop, loading data to `__local` memory once.
 *   - A single `__local` memory buffer is used for fast, in-place updates.
 *   - Barriers between color sweeps synchronize the workgroup to ensure correctness.
 *   - This design minimizes global memory traffic for maximum performance.
 */

__kernel void GS_diff(
    __global float4*      dp,            // 1 : [nverts,4] in/out displacements dp = p - p0
    __global const int*   neighs,        // 2 : [nverts,MAX_NEIGHBORS] neighbor indices
    __global const float* kngs,          // 3 : [nverts,MAX_NEIGHBORS] stiffness k_ij for each neighbor
    __global const float4* bvec,         // 4 : [nverts,4] Precomputed b'.xyz and Aii'.w
    const int nverts,                    // 5 : Number of vertices
    const int nIters,                    // 6 : Number of iterations
    __global const int*   color_list,
    __global const int*   color_offsets,
    const int             num_colors
){
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);

    __local float4 ldp[1024]; // Single buffer for in-place updates
    for (int i=lid; i<nverts; i+=lsize) { ldp[i] = dp[i]; }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int iter = 0; iter < nIters; ++iter) {
        for (int c = 0; c < num_colors; ++c) { // Sweep through colors
            int off         = color_offsets[c];
            int off_next    = color_offsets[c+1];
            int color_count = off_next - off;

            for (int local_idx = lid; local_idx < color_count; local_idx += lsize) {
                int vid = color_list[off + local_idx];
                float3 sum_j = (float3)(0.0f);
                for (int e=0; e<MAX_NEIGHBORS; ++e) {
                    int idx = vid * MAX_NEIGHBORS + e;
                    int j   = neighs[idx];
                    if (j < 0) break;
                    float  k   = kngs[idx];
                    float3 dpj = ldp[j].xyz; // Read most recent value from this iteration
                    sum_j += dpj * k; // sum_j += k_ij * dp_j
                }
                const float4 bi = bvec[vid]; // bi.xyz is RHS, bi.w is A_ii
                float3 dpi_new = (bi.xyz + sum_j) * (1.0f / bi.w);
                ldp[vid] = (float4)(dpi_new, 0.0f); // In-place update
            }
            barrier(CLK_LOCAL_MEM_FENCE); // Sync before processing next color
        }
    }
    for (int i=lid; i<nverts; i+=lsize) { dp[i]=ldp[i]; } // Write final result back to global
}




/**
 * @brief kernel precompute_dRHS : Pre-computes RHS and diagonal for the linear solver (float, OPTIMIZED).
 *
 * This high-performance version loads all vertex positions into `__local` memory once at the start.
 * All subsequent neighbor position lookups during force calculation are served from this fast on-chip
 * cache, drastically reducing global memory traffic and latency.
 *
 * @param bvec      Output buffer for computed results (float4 per vertex).
 * @param pos       Input buffer of current vertex positions and mass (pos.xyz, pos.w = mass).
 * @param kfix      Additional stiffness for fixed/constrained vertices (pass a zeroed buffer if not used).
 */
__kernel void precompute_dRHS(
    __global float4*      bvec,          // 1 : [nverts,4] Precomputed b'.xyz and Aii'.w
    __global const float4* pos,          // 2 : [nverts,4] Current vertex positions and mass (pos.xyz, pos.w = mass)
    __global const int*   neighs,        // 3 : [nverts,MAX_NEIGHBORS] neighbor indices
    __global const float* kngs,          // 4 : [nverts,MAX_NEIGHBORS] stiffness k_ij for each neighbor
    __global const float* l0ngs,         // 5 : [nverts,MAX_NEIGHBORS] rest length l0_ij for each neighbor
    __global const float* kfix,          // 6 : [nverts] Additional stiffness for fixed/constrained vertices
    const int nverts,                    // 7 : Number of vertices
    const float dt                       // 8 : Time step
){
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);
    __local float4 lpos[1024];

    for (int i=lid; i<nverts; i+=lsize) { lpos[i] = pos[i]; } // Load all positions into fast local memory
    barrier(CLK_LOCAL_MEM_FENCE); // Ensure all loads are complete before proceeding

    for (int i=lid; i<nverts; i+=lsize) { // Each thread computes its assigned vertices
        const float4 pi_mass = lpos[i]; // Read own position from local memory
        const float3 pi = pi_mass.xyz;
        const float  mi = pi_mass.w;

        const float inv_dt2 = 1.0f / (dt * dt);
        float Ii = mi * inv_dt2 + kfix[i];
        float3 bi = (float3)(0.0f);
        float  Aii = Ii;

        for (int jj=0; jj<MAX_NEIGHBORS; ++jj) {
            int idx = i * MAX_NEIGHBORS + jj;
            int j   = neighs[idx];
            if (j < 0) break;
            float  k   = kngs[idx];
            float  l0  = l0ngs[idx];
            float3 pj  = lpos[j].xyz; // Read neighbor position from fast local memory
            float3 dij = pi - pj;
            float  len = length(dij);
            bi  += dij * (k * (l0 / len - 1.0f));
            Aii += k;
        }
        bvec[i] = (float4)(bi, Aii);
    }
}



#ifdef cl_khr_fp64

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

/**
 * @brief kernel precompute_dRHS_d : Pre-computes RHS and diagonal for the linear solver (double, OPTIMIZED).
 *
 * High-performance version using `__local` memory to cache vertex positions. All critical
 * geometry calculations are done in double precision to avoid numerical issues, reading from
 * a `double4` position buffer.
 *
 * @param bvec      Output buffer for computed results (float4 per vertex).
 * @param pos       Input buffer of current vertex positions and mass (double4: pos.xyz, pos.w = mass).
 * @param kfix      Additional stiffness for fixed/constrained vertices (pass a zeroed buffer if not used).
 */
__kernel void precompute_dRHS_d(
    __global float4*        bvec,          // 1 : [nverts,4] Precomputed b'.xyz and Aii'.w
    __global const double4*  pos,          // 2 : [nverts,4] Current vertex positions and mass (pos.xyz, pos.w = mass)
    __global const int*     neighs,        // 3 : [nverts,MAX_NEIGHBORS] neighbor indices
    __global const float*   kngs,          // 4 : [nverts,MAX_NEIGHBORS] stiffness k_ij for each neighbor
    __global const float*   l0ngs,         // 5 : [nverts,MAX_NEIGHBORS] rest length l0_ij for each neighbor
    __global const float*   kfix,          // 6 : [nverts] Additional stiffness for fixed/constrained vertices
    const int nverts,                      // 7 : Number of vertices
    const float dt                         // 8 : Time step
){
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);
    __local double4 lpos[1024];

    for (int i=lid; i<nverts; i+=lsize) { lpos[i] = pos[i]; } // Load all double-precision positions into local memory
    barrier(CLK_LOCAL_MEM_FENCE); // Ensure all loads complete

    for (int i=lid; i<nverts; i+=lsize) { // Each thread computes its assigned vertices
        const double4 pi_mass = lpos[i]; // Read own position from local memory
        const double3 pi = pi_mass.xyz;
        const double  mi = pi_mass.w;

        const double inv_dt2 = 1.0 / ((double)dt * (double)dt);
        double Ii = mi * inv_dt2 + (double)kfix[i];
        double3 bi = (double3)(0.0);
        double  Aii = Ii;

        for (int jj=0; jj<MAX_NEIGHBORS; ++jj) {
            int idx = i * MAX_NEIGHBORS + jj;
            int j   = neighs[idx];
            if (j < 0) break;
            double k   = (double)kngs[idx];
            double l0  = (double)l0ngs[idx];
            double3 pj = lpos[j].xyz; // Read neighbor position from fast local memory
            double3 dij = pi - pj;
            double  len = length(dij);
            bi  += dij * (k * (l0 / len - 1.0));
            Aii += k;
        }
        bvec[i] = (float4)((float)bi.x, (float)bi.y, (float)bi.z, (float)Aii); // Cast back to float4 for storage
    }
}

#endif // cl_khr_fp64
