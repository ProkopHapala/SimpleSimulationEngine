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
    __global       float4* x,           // [npoint,4] solution candidate for position + mass {x,y,z,mass}, both input and output as Gauss-Seidel operates in-place
    __global const float4* b,           // [npoint,4] RHS vector b {x,y,z, ?}
    __global const int*   neighs,       // [npoint,nneigh_max] Neighbor indices, padded with -1
    __global const float* kngs,         // [npoint,nneigh_max] Non-zero matrix elements (stiffness)
    __global const float* Aii,          // [npoint] Diagonal elements of the matrix A
    __global const int*   color_group,  // [npoint] Node indices for the current color group
    __global       float4* r,           // [npoint,4] Output residual vector {x,y,z,?}
    const int nmax_neigh,               // Maximum number of neighbors (padding size)
    const int n,                        // Total number of nodes
    const int group_start,              // Starting index in the color group buffer
    const int group_count               // Number of nodes in this color group
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
    __global float4* x,
    __global const float4* y,
    __global const int* chunk_vertices,
    __global const int* chunk_neighbors,
    __global const int* slot_table,
    __global const int* neighs,
    __global const float* kngs,
    __global const float* rest_lengths,
    const int vcount,
    const int ncount,
    const int nmax,
    const float inv_h2,
    const float det_eps
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
    const float4 xi = loc_x[lid];
    const float4 yi = loc_y[lid];

    float3 grad = (float3)(0.0f);
    float H[9];
    H[0]=H[4]=H[8]=inv_h2 * xi.w;
    H[1]=H[2]=H[3]=H[5]=H[6]=H[7]=0.0f;

    const int base = vid * nmax;
    const int slot_base = lid * nmax;

    for (int jj = 0; jj < nmax; ++jj){
        const int packed = loc_neighbor_index[slot_base + jj];
        if (packed < 0) break;

        const int nslot = packed;
        float4 xj = loc_neighbors[nslot];
        const int nj = chunk_neighbors[nslot];

        const float k = kngs[base + jj];
        const float rest = rest_lengths[base + jj];

        float3 d = xi.xyz - xj.xyz;
        float len = length(d);
        if (len > 1e-7f){
            float3 dir = d / len;
            float coeff = k * (len - rest);
            grad += coeff * dir;
            float3 outer0 = dir * dir.x;
            float3 outer1 = dir * dir.y;
            float3 outer2 = dir * dir.z;
            float sdiag = k * rest / len;
            H[0] += k * (1.0f - dir.x * dir.x) + sdiag * dir.x * dir.x;
            H[4] += k * (1.0f - dir.y * dir.y) + sdiag * dir.y * dir.y;
            H[8] += k * (1.0f - dir.z * dir.z) + sdiag * dir.z * dir.z;
            H[1] -= k * dir.x * dir.y - sdiag * dir.x * dir.y;
            H[2] -= k * dir.x * dir.z - sdiag * dir.x * dir.z;
            H[5] -= k * dir.y * dir.z - sdiag * dir.y * dir.z;
            H[3] = H[1];
            H[6] = H[2];
            H[7] = H[5];
        }
    }

    grad += inv_h2 * (xi.xyz - yi.xyz);
    float c00 = H[4]*H[8] - H[5]*H[7];
    float c01 = H[2]*H[7] - H[1]*H[8];
    float c02 = H[1]*H[5] - H[2]*H[4];
    float c10 = H[5]*H[6] - H[3]*H[8];
    float c11 = H[0]*H[8] - H[2]*H[6];
    float c12 = H[2]*H[3] - H[0]*H[5];
    float c20 = H[3]*H[7] - H[4]*H[6];
    float c21 = H[1]*H[6] - H[0]*H[7];
    float c22 = H[0]*H[4] - H[1]*H[3];

    float det = H[0]*c00 + H[1]*c10 + H[2]*c20;
    if (fabs(det) < det_eps) return;

    float inv_det = 1.0f/det;
    float dx0 = (c00*grad.x + c01*grad.y + c02*grad.z) * inv_det;
    float dx1 = (c10*grad.x + c11*grad.y + c12*grad.z) * inv_det;
    float dx2 = (c20*grad.x + c21*grad.y + c22*grad.z) * inv_det;

    x[vid].xyz = xi.xyz - (float3)(dx0, dx1, dx2);
}

__kernel void vbd_vertex_serial(
    __global float4* x,
    __global const float4* y,
    __global const int* neighs,
    __global const float* kngs,
    __global const float* rest_lengths,
    const int n_points,
    const int nmax,
    const float inv_h2,
    const float det_eps
){
    if (get_global_id(0) != 0) return;

    for (int vid = 0; vid < n_points; ++vid){
        float4 xi = x[vid];
        float4 yi = y[vid];

        float3 grad = (float3)(0.0f);
        float H[9];
        H[0]=H[4]=H[8]=inv_h2 * xi.w;
        H[1]=H[2]=H[3]=H[5]=H[6]=H[7]=0.0f;

        const int base = vid * nmax;
        for (int jj = 0; jj < nmax; ++jj){
            const int nb = neighs[base + jj];
            if (nb < 0) break;
            const float k = kngs[base + jj];
            const float rest = rest_lengths[base + jj];

            float3 d = xi.xyz - x[nb].xyz;
            float len = length(d);
            if (len > 1e-7f){
                float3 dir = d / len;
                grad += k * (len - rest) * dir;
                float sdiag = k * rest / len;
                H[0] += k * (1.0f - dir.x * dir.x) + sdiag * dir.x * dir.x;
                H[4] += k * (1.0f - dir.y * dir.y) + sdiag * dir.y * dir.y;
                H[8] += k * (1.0f - dir.z * dir.z) + sdiag * dir.z * dir.z;
                H[1] -= k * dir.x * dir.y - sdiag * dir.x * dir.y;
                H[2] -= k * dir.x * dir.z - sdiag * dir.x * dir.z;
                H[5] -= k * dir.y * dir.z - sdiag * dir.y * dir.z;
                H[3] = H[1];
                H[6] = H[2];
                H[7] = H[5];
            }
        }

        grad += inv_h2 * (xi.xyz - yi.xyz);

        float c00 = H[4]*H[8] - H[5]*H[7];
        float c01 = H[2]*H[7] - H[1]*H[8];
        float c02 = H[1]*H[5] - H[2]*H[4];
        float c10 = H[5]*H[6] - H[3]*H[8];
        float c11 = H[0]*H[8] - H[2]*H[6];
        float c12 = H[2]*H[3] - H[0]*H[5];
        float c20 = H[3]*H[7] - H[4]*H[6];
        float c21 = H[1]*H[6] - H[0]*H[7];
        float c22 = H[0]*H[4] - H[1]*H[3];

        float det = H[0]*c00 + H[1]*c10 + H[2]*c20;
        if (fabs(det) < det_eps) continue;

        float inv_det = 1.0f / det;
        float dx0 = (c00*grad.x + c01*grad.y + c02*grad.z) * inv_det;
        float dx1 = (c10*grad.x + c11*grad.y + c12*grad.z) * inv_det;
        float dx2 = (c20*grad.x + c21*grad.y + c22*grad.z) * inv_det;

        x[vid].xyz = xi.xyz - (float3)(dx0, dx1, dx2);
    }
}
