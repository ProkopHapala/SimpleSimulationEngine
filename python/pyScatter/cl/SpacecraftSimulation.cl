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
