
// -- GEOMETRY HELPERS --

// Intersection with Finite Cylinder (Tube)
// Returns path length inside the tube
float intersect_tube(float3 ray_org, float3 ray_dir, float3 p0, float3 p1, float radius) {
    float3 ba = p1 - p0;
    float3 oc = ray_org - p0;
    
    float baba = dot(ba, ba);
    float bard = dot(ba, ray_dir);
    float baoc = dot(ba, oc);
    
    float k2 = baba - bard*bard;
    float k1 = baba * dot(oc, ray_dir) - baoc*bard;
    float k0 = baba * dot(oc, oc) - baoc*baoc - radius*radius*baba;
    
    float h = k1*k1 - k2*k0;
    
    if (h < 0.0f) return 0.0f;
    
    h = sqrt(h);
    float t_entry = (-k1 - h) / k2;
    float t_exit  = (-k1 + h) / k2;
    
    // Infinite cylinder intersections found. Now clip to segment.
    float y_entry = baoc + t_entry * bard;
    float y_exit  = baoc + t_exit  * bard;
    
    // We only care about the part of the ray strictly inside the segment 0..baba
    // This is a simplified "caps-less" check (open ended tube) for performance
    if (y_entry < 0.0f && y_exit < 0.0f) return 0.0f;
    if (y_entry > baba && y_exit > baba) return 0.0f;
    
    // Clamp t to the segment limits if needed (simple approx)
    // A robust engine would calculate cap intersections, 
    // here we assume long slender tubes where caps don't matter much.
    return (t_exit - t_entry); 
}

// Intersection with Thin Triangle (Sheet)
// Returns path length: Thickness / dot(normal, ray)
float intersect_tri(float3 ray_org, float3 ray_dir, float3 v0, float3 v1, float3 v2, float thickness) {
    const float EPSILON = 1e-7f;
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
    
    float t = f * dot(edge2, q);
    
    if (t > EPSILON) {
        // Geometric hit. Calculate effective path length through the sheet.
        // L = thickness / |cos(theta)| = thickness / |dot(ray, normal)|
        // normal is normalized cross(edge1, edge2)
        float3 normal = normalize(cross(edge1, edge2));
        float cos_theta = fabs(dot(ray_dir, normal));
        return thickness / (cos_theta + 1e-6f); // Avoid div by zero
    }
    return 0.0f;
}

// -- MAIN KERNEL --

__kernel void xray_simulation(
    __global float4* prim_data,      // 3 float4s per object
    __global float4* cluster_spheres,// {x,y,z,R}
    __global int2*   cluster_ranges, // {start_index, count}
    __global float*  output_image,
    int n_clusters,
    float3 source_pos,
    float  source_radius,
    float3 det_origin,
    float3 det_u,
    float3 det_v,
    int width
) {
    // 1. Thread ID and Pixel Ray
    int lx = get_local_id(0);
    int ly = get_local_id(1);
    int l_idx = ly * get_local_size(0) + lx; // 0..255
    
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    
    float3 pixel_pos = det_origin + (float)gx * det_u + (float)gy * det_v;
    float3 ray_dir = normalize(pixel_pos - source_pos);
    float  ray_len = length(pixel_pos - source_pos); // Max distance

    // 2. Define Workgroup Cone (Frustum approx)
    // Center of the pixel tile on detector
    float3 wg_center_det = det_origin 
         + (get_group_id(0) * get_local_size(0) + get_local_size(0)/2.0f) * det_u 
         + (get_group_id(1) * get_local_size(1) + get_local_size(1)/2.0f) * det_v;
    
    // Approximate radius of the tile
    float wg_radius_det = length(det_u) * get_local_size(0) * 0.75f; 

    // Axis of the cone
    float3 cone_axis = normalize(wg_center_det - source_pos);

    // 3. Local Memory Setup
    // Max primitives in cache (must match python cluster_size limit)
    #define CACHE_CAPACITY 64 
    
    // We store 3 float4s per primitive. 
    // Layout: [prim0_row0, prim0_row1, prim0_row2, prim1_row0...]
    __local float4 l_prim_cache[CACHE_CAPACITY * 3]; 
    __local int    l_active_vote[256];

    float total_mu_dist = 0.0f; // Accumulated attenuation (Sigma * Dist)

    // 4. Batched Cluster Processing
    // Loop over all clusters in steps of 256 (workgroup size)
    for (int batch_base = 0; batch_base < n_clusters; batch_base += 256) {
        
        int c_idx = batch_base + l_idx;
        bool active = false;

        // A. Broad Phase: Cone-Sphere Intersection
        if (c_idx < n_clusters) {
            float4 sphere = cluster_spheres[c_idx];
            float3 c_cen = sphere.xyz;
            float  c_rad = sphere.w;

            // Project sphere center onto cone axis
            float3 vec = c_cen - source_pos;
            float t_proj = dot(vec, cone_axis);
            
            // Distance from axis
            float dist_sq = dot(vec, vec) - t_proj*t_proj;
            float dist = sqrt(max(0.0f, dist_sq));

            // Cone radius at this distance (linear interp)
            // t_norm 0.0 at source, 1.0 at detector
            float t_norm = t_proj / length(wg_center_det - source_pos);
            float cone_r = (1.0f - t_norm)*source_radius + t_norm*wg_radius_det;

            if (dist < (cone_r + c_rad)) {
                active = true;
            }
        }
        
        l_active_vote[l_idx] = active ? 1 : 0;
        barrier(CLK_LOCAL_MEM_FENCE);

        // B. Narrow Phase: Process Active Clusters
        for (int k = 0; k < 256; k++) {
            if (l_active_vote[k]) {
                int active_c_idx = batch_base + k;
                int2 range = cluster_ranges[active_c_idx];
                int start_prim = range.x;
                int count_prim = range.y;

                // Cooperative Load: Global -> Local
                // We need to copy count_prim * 3 float4 vectors.
                int total_vectors = count_prim * 3;
                int num_fetches = (total_vectors + 255) / 256;

                for (int f = 0; f < num_fetches; f++) {
                    int vec_id = f * 256 + l_idx;
                    if (vec_id < total_vectors) {
                         // Read simply as a flat array of float4
                         l_prim_cache[vec_id] = prim_data[start_prim * 3 + vec_id];
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE);

                // Ray Trace against Cache
                for (int p = 0; p < count_prim; p++) {
                    // Unpack from local memory (register pressure is low here)
                    float4 r0 = l_prim_cache[p*3 + 0]; // {x,y,z, type}
                    float4 r1 = l_prim_cache[p*3 + 1]; // {x,y,z, R/Thick}
                    float4 r2 = l_prim_cache[p*3 + 2]; // {x,y,z, mu} (xyz used for tri v2)

                    float type = r0.w;
                    float mu   = r2.w;
                    float path = 0.0f;

                    if (type > 0.5f) { 
                        // Triangle
                        path = intersect_tri(source_pos, ray_dir, r0.xyz, r1.xyz, r2.xyz, r1.w);
                    } else {
                        // Tube
                        path = intersect_tube(source_pos, ray_dir, r0.xyz, r1.xyz, r1.w);
                    }

                    if (path > 0.0f) {
                         total_mu_dist += path * mu;
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write Result (Beer-Lambert Law)
    output_image[gy * width + gx] = exp(-total_mu_dist);
}