// OpenCL Kernel for Spacecraft Simulation (Occlusion & Scattering)

#ifdef DEBUG
#define PRINTF(a) printf a
#else
#define PRINTF(a)
#endif

// -- Helpers --

// Intersection with Thin Triangle
// returns 1.0 if hit, 0.0 if missed.
float ray_triangle_intersect(float3 ray_org, float3 ray_dir, float3 v0, float3 v1, float3 v2, float* t) {
    const float EPSILON = 1e-6f;
    float3 edge1 = v1 - v0;
    float3 edge2 = v2 - v0;
    float3 h = cross(ray_dir, edge2);
    float a = dot(edge1, h);
    
    if (a > -EPSILON && a < EPSILON) return 0.0f; // Parallel
    
    float f = 1.0f / a;
    float3 s = ray_org - v0;
    float u = f * dot(s, h);
    
    if (u < 0.0f || u > 1.0f) return 0.0f;
    
    float3 q = cross(s, edge1);
    float v = f * dot(ray_dir, q);
    
    if (v < 0.0f || u + v > 1.0f) return 0.0f;
    
    *t = f * dot(edge2, q);
    
    if (*t > EPSILON) {
        return 1.0f;
    }
    return 0.0f;
}

// 1. Occlusion Kernel (Brute force for Radiosity)
__kernel void compute_occlusion(
    __global const float4* points, // x,y,z, face_id
    __global const float4* tris,   // A, B, C sequential, w=face_id
    __global float* occ_matrix,
    int n_points,
    int n_tris
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    
    if (i >= n_points || j >= n_points) return;
    
    // Diagonal is no occlusion (or self, don't care for radiosity since F_ii = 0)
    if (i == j) {
        occ_matrix[i * n_points + j] = 0.0f;
        return;
    }
    
    float4 p_i = points[i];
    float4 p_j = points[j];
    
    float3 ray_org = p_i.xyz;
    float3 diff = p_j.xyz - p_i.xyz;
    float dist = length(diff);
    float3 ray_dir = diff / dist;
    
    int face_i = (int)(p_i.w);
    int face_j = (int)(p_j.w);
    
    float is_occluded = 0.0f;
    
    for (int t_idx = 0; t_idx < n_tris; t_idx++) {
        float4 v0 = tris[t_idx * 3 + 0];
        float4 v1 = tris[t_idx * 3 + 1];
        float4 v2 = tris[t_idx * 3 + 2];
        
        int face_t = (int)(v0.w);
        
        // Skip self-occlusion with the source or destination faces
        if (face_t == face_i || face_t == face_j) continue;
        
        float t_hit = 0.0f;
        if (ray_triangle_intersect(ray_org, ray_dir, v0.xyz, v1.xyz, v2.xyz, &t_hit) > 0.5f) {
            // Check if hit is strictly between i and j
            if (t_hit > 1e-4f && t_hit < dist - 1e-4f) {
                is_occluded = 1.0f;
                break; // Early exit
            }
        }
    }
    
    occ_matrix[i * n_points + j] = is_occluded;
}

// Ray vs sphere intersection for a segment
// Returns true if segment [org, org+dir*dist] intersects sphere
bool ray_sphere_intersect(float3 org, float3 dir, float dist, float3 center, float radius) {
    float3 oc = org - center;
    float b = dot(oc, dir);
    float c = dot(oc, oc) - radius * radius;
    float h = b * b - c;
    if (h < 0.0f) return false;
    h = sqrt(h);
    return (-b + h > 0.0f) && (-b - h < dist);
}

// 1b. Tiled Occlusion Kernel with cluster culling + local memory
// Single-level bounding sphere clusters — GPU-friendly alternative to multi-level BVH.
// Workgroup TILE×TILE processes a block of (i,j) pairs.
// Per cluster: vote-based broad-phase (ray vs sphere), cooperative local-memory load, narrow-phase.
#define TILE 8
#define CLUSTER_CAPACITY 32

__kernel void compute_occlusion_tiled(
    __global const float4* points,          // x,y,z, face_id
    __global const float4* tris,            // A,B,C sequential, w=face_id (spatially sorted)
    __global const float4* cluster_spheres, // x,y,z, radius
    __global const int2*   cluster_ranges,  // {start_tri_idx, count}
    __global float* occ_matrix,
    int n_points,
    int n_clusters
) {
    int li = get_local_id(0);
    int lj = get_local_id(1);
    int l_idx = lj * TILE + li;

    int gi = get_group_id(0) * TILE + li;
    int gj = get_group_id(1) * TILE + lj;

    // Load points for this tile into local memory
    __local float4 l_pi[TILE];
    __local float4 l_pj[TILE];
    if (lj == 0 && gi < n_points) l_pi[li] = points[gi];
    if (li == 0 && gj < n_points) l_pj[lj] = points[gj];
    barrier(CLK_LOCAL_MEM_FENCE);

    // Valid flag — no early returns (all threads must hit all barriers)
    bool valid = (gi < n_points) && (gj < n_points) && (gi != gj);

    float3 ray_org = (float3)(0.0f);
    float3 ray_dir = (float3)(0.0f);
    float dist = 0.0f;
    int face_i = -1, face_j = -1;

    if (valid) {
        float4 p_i = l_pi[li];
        float4 p_j = l_pj[lj];
        ray_org = p_i.xyz;
        float3 diff = p_j.xyz - p_i.xyz;
        dist = length(diff);
        ray_dir = diff / dist;
        face_i = (int)(p_i.w);
        face_j = (int)(p_j.w);
    }

    // Local memory: triangle cache + vote array
    __local float4 l_tris[CLUSTER_CAPACITY * 3];
    __local int l_vote[TILE * TILE];

    float is_occluded = 0.0f;

    for (int c = 0; c < n_clusters; c++) {
        // Broad phase: does my ray hit this cluster's bounding sphere?
        bool my_active = false;
        if (valid && is_occluded < 0.5f) {
            float4 sphere = cluster_spheres[c];
            my_active = ray_sphere_intersect(ray_org, ray_dir, dist, sphere.xyz, sphere.w);
        }

        // Vote: does ANY thread in workgroup need this cluster?
        l_vote[l_idx] = my_active ? 1 : 0;
        barrier(CLK_LOCAL_MEM_FENCE);

        bool any_active = false;
        for (int v = 0; v < TILE * TILE; v++) {
            if (l_vote[v]) { any_active = true; break; }
        }

        if (any_active) {
            // Cooperative load: cluster triangles → local memory
            int2 range = cluster_ranges[c];
            int total_vecs = range.y * 3;
            for (int f = l_idx; f < total_vecs; f += TILE * TILE) {
                l_tris[f] = tris[range.x * 3 + f];
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            // Narrow phase: test my ray against cached triangles
            if (my_active) {
                for (int t = 0; t < range.y; t++) {
                    float4 v0 = l_tris[t * 3 + 0];
                    int face_t = (int)(v0.w);
                    if (face_t == face_i || face_t == face_j) continue;

                    float4 v1 = l_tris[t * 3 + 1];
                    float4 v2 = l_tris[t * 3 + 2];

                    float t_hit = 0.0f;
                    if (ray_triangle_intersect(ray_org, ray_dir, v0.xyz, v1.xyz, v2.xyz, &t_hit) > 0.5f) {
                        if (t_hit > 1e-4f && t_hit < dist - 1e-4f) {
                            is_occluded = 1.0f;
                            break;
                        }
                    }
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

    if (gi < n_points && gj < n_points) {
        occ_matrix[gi * n_points + gj] = (gi == gj) ? 0.0f : is_occluded;
    }
}

#define MAX_CHANNELS 32
#define MAX_CANDIDATES 64
#define SPARSE_WG 64

__kernel void compute_sparse_channels(
    __global const float4* points,
    __global const float4* norm_area,
    __global const float4* tris,
    __global const float4* cluster_spheres,
    __global const int2*   cluster_ranges,
    __global int*   out_idx,
    __global float* out_w,
    int n_points,
    int n_clusters,
    int row_offset,
    int n_rows,
    int kmax,
    float min_weight,
    float max_dist
) {
    int lid = get_local_id(0);
    int lsz = get_local_size(0);
    int row = get_global_id(0);
    int i = row_offset + row;
    bool valid_i = (row < n_rows) && (i < n_points);

    float4 pi = (float4)(0.0f);
    float4 ni_a = (float4)(0.0f);
    int face_i = -1;
    if (valid_i) {
        pi = points[i];
        ni_a = norm_area[i];
        face_i = (int)(pi.w);
    }

    int kcand = kmax * 4;
    if (kcand < kmax) kcand = kmax;
    if (kcand > MAX_CANDIDATES) kcand = MAX_CANDIDATES;
    if (kmax > MAX_CHANNELS) kmax = MAX_CHANNELS;

    int cand_j[MAX_CANDIDATES];
    float cand_w[MAX_CANDIDATES];
    float cand_key[MAX_CANDIDATES];
    int best_j[MAX_CHANNELS];
    float best_w[MAX_CHANNELS];
    float best_key[MAX_CHANNELS];

    for (int s = 0; s < MAX_CANDIDATES; s++) { cand_j[s] = -1; cand_w[s] = 0.0f; cand_key[s] = -1.0f; }
    for (int s = 0; s < MAX_CHANNELS; s++) { best_j[s] = -1; best_w[s] = 0.0f; best_key[s] = -1.0f; }

    float max_d2 = max_dist * max_dist;
    for (int j = 0; j < n_points; j++) {
        if (!valid_i || j == i) continue;
        float4 pj = points[j];
        float3 d = pj.xyz - pi.xyz;
        float d2 = dot(d, d) + 1e-12f;
        if (max_dist > 0.0f && d2 > max_d2) continue;
        float invd = rsqrt(d2);
        float3 rhat = d * invd;
        float ci = dot(ni_a.xyz, rhat);
        float4 nj_a = norm_area[j];
        float cj = dot(nj_a.xyz, -rhat);
        float key = fabs(ci) * fabs(cj) * nj_a.w / (3.14159265358979323846f * d2);
        if (key <= min_weight) continue;
        int smin = 0;
        float vmin = cand_key[0];
        for (int s = 1; s < MAX_CANDIDATES; s++) {
            if (s >= kcand) break;
            if (cand_key[s] < vmin) { vmin = cand_key[s]; smin = s; }
        }
        if (key > vmin) {
            cand_j[smin] = j;
            cand_w[smin] = (ci >= 0.0f) ? key : -key;
            cand_key[smin] = key;
        }
    }

    __local float4 l_tris[CLUSTER_CAPACITY * 3];
    __local int l_vote[SPARSE_WG];

    for (int cs = 0; cs < MAX_CANDIDATES; cs++) {
        bool valid_pair = valid_i && (cs < kcand) && (cand_j[cs] >= 0);
        int j = valid_pair ? cand_j[cs] : 0;
        float4 pj = valid_pair ? points[j] : (float4)(0.0f);
        int face_j = valid_pair ? (int)(pj.w) : -2;
        float3 ray_org = pi.xyz;
        float3 diff = pj.xyz - pi.xyz;
        float dist = length(diff);
        float3 ray_dir = (dist > 1e-12f) ? (diff / dist) : (float3)(0.0f);
        bool visible = valid_pair;

        for (int c = 0; c < n_clusters; c++) {
            bool my_active = false;
            if (visible) {
                float4 sphere = cluster_spheres[c];
                my_active = ray_sphere_intersect(ray_org, ray_dir, dist, sphere.xyz, sphere.w);
            }

            if (lid < SPARSE_WG) l_vote[lid] = my_active ? 1 : 0;
            barrier(CLK_LOCAL_MEM_FENCE);

            bool any_active = false;
            for (int v = 0; v < SPARSE_WG; v++) {
                if (v >= lsz) break;
                if (l_vote[v]) { any_active = true; break; }
            }

            if (any_active) {
                int2 range = cluster_ranges[c];
                int total_vecs = range.y * 3;
                for (int f = lid; f < total_vecs; f += lsz) {
                    l_tris[f] = tris[range.x * 3 + f];
                }
                barrier(CLK_LOCAL_MEM_FENCE);

                if (my_active) {
                    for (int t = 0; t < range.y; t++) {
                        float4 v0 = l_tris[t * 3 + 0];
                        int face_t = (int)(v0.w);
                        if (face_t == face_i || face_t == face_j) continue;
                        float4 v1 = l_tris[t * 3 + 1];
                        float4 v2 = l_tris[t * 3 + 2];
                        float t_hit = 0.0f;
                        if (ray_triangle_intersect(ray_org, ray_dir, v0.xyz, v1.xyz, v2.xyz, &t_hit) > 0.5f) {
                            if (t_hit > 1e-4f && t_hit < dist - 1e-4f) { visible = false; break; }
                        }
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
        }

        if (visible) {
            float key = cand_key[cs];
            int smin = 0;
            float vmin = best_key[0];
            for (int s = 1; s < MAX_CHANNELS; s++) {
                if (s >= kmax) break;
                if (best_key[s] < vmin) { vmin = best_key[s]; smin = s; }
            }
            if (key > vmin) {
                best_j[smin] = cand_j[cs];
                best_w[smin] = cand_w[cs];
                best_key[smin] = key;
            }
        }
    }

    if (valid_i) {
        for (int s = 0; s < MAX_CHANNELS; s++) {
            if (s >= kmax) break;
            out_idx[row * kmax + s] = best_j[s];
            out_w[row * kmax + s] = best_w[s];
        }
    }
}

// 2. Ionizing Radiation Attenuation Kernel
// Simulates rays from a source to a detector grid, computing attenuation through triangles.
__kernel void compute_attenuation(
    __global const float4* tris,    // A,B,C, w=face_id
    __global const float*  mu_vals, // attenuation coefficient per triangle (or material)
    __global float* image,          // detector image
    int n_tris,
    float3 source_pos,
    float3 det_origin,
    float3 det_u,
    float3 det_v,
    int width,
    int height
) {
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    
    if (gx >= width || gy >= height) return;
    
    float3 pixel_pos = det_origin + (float)gx * det_u + (float)gy * det_v;
    float3 ray_dir = normalize(pixel_pos - source_pos);
    
    float total_optical_depth = 0.0f;
    
    for (int t_idx = 0; t_idx < n_tris; t_idx++) {
        float4 v0 = tris[t_idx * 3 + 0];
        float4 v1 = tris[t_idx * 3 + 1];
        float4 v2 = tris[t_idx * 3 + 2];
        
        float t_hit = 0.0f;
        if (ray_triangle_intersect(source_pos, ray_dir, v0.xyz, v1.xyz, v2.xyz, &t_hit) > 0.5f) {
            // Penetration depth through a thin sheet
            // L = thickness / |cos(theta)| = thickness / |dot(ray, normal)|
            float3 edge1 = v1.xyz - v0.xyz;
            float3 edge2 = v2.xyz - v0.xyz;
            float3 normal = normalize(cross(edge1, edge2));
            float cos_theta = fabs(dot(ray_dir, normal));
            
            float thickness = 0.01f; // Default thin sheet
            float path_len = thickness / (cos_theta + 1e-6f);
            
            float mu = mu_vals[t_idx];
            total_optical_depth += mu * path_len;
        }
    }
    
    image[gy * width + gx] = exp(-total_optical_depth);
}
