// Pack (x, y) into uint32
uint pack_coord(int x, int y) {
    return (uint)x | ((uint)y << 16);
}

// Unpack packed coord
int2 unpack_coord(uint p) {
    return (int2)(p & 0xFFFF, (p >> 16) & 0xFFFF);
}

// Local Relaxation Kernel (Single Buffer, Sync via Registers)
__kernel void solve_tiles(
    __global const float* heightmap,
    __global uint* parent_map,
    __global float* cost_map,
    __global uint* tile_sinks,
    const int width,
    const int height,
    const float K)
{
    int lx = get_local_id(0);
    int ly = get_local_id(1);
    int gx = get_group_id(0) * 16 + lx;
    int gy = get_group_id(1) * 16 + ly;
    int tile_id = get_group_id(1) * (width/16) + get_group_id(0);

    __local float local_h[16][16];
    
    // Single Buffer (Memory Optimized)
    __local float lc[16][16]; 
    __local uint lp[16][16]; 

    // 1. Load Data
    if (gx < width && gy < height) {
        local_h[ly][lx] = heightmap[gy * width + gx];
    } else {
        local_h[ly][lx] = 1e10f;
    }

    // 2. Init Sink
    float centerX = 7.5f; float centerY = 7.5f;
    float distSq = (lx - centerX)*(lx - centerX) + (ly - centerY)*(ly - centerY);
    float potential_h = local_h[ly][lx] + K * distSq;
    
    // Write initial state
    lc[ly][lx] = potential_h;
    lp[ly][lx] = pack_coord(gx, gy);
    barrier(CLK_LOCAL_MEM_FENCE);

    // Reduction for sink
    __local int2 min_pos;
    if (lx == 0 && ly == 0) {
        float min_val = 1e10f;
        for(int j=0; j<16; j++) {
            for(int i=0; i<16; i++) {
                if(lc[j][i] < min_val) {
                    min_val = lc[j][i];
                    min_pos = (int2)(i, j);
                }
            }
        }
        tile_sinks[tile_id] = pack_coord(get_group_id(0) * 16 + min_pos.x, get_group_id(1) * 16 + min_pos.y);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 3. Relaxation Init
    int sink_lx = min_pos.x;
    int sink_ly = min_pos.y;
    
    float start_cost = (lx == sink_lx && ly == sink_ly) ? local_h[ly][lx] : 1e10f;
    lc[ly][lx] = start_cost;
    barrier(CLK_LOCAL_MEM_FENCE);

    // 4. Synchronous Iteration (Simulated Double Buffer via Registers)
    for (int it = 0; it < 64; it++) {
        float my_h = local_h[ly][lx];
        
        // We calculate the NEXT state into private registers
        // This ensures we don't overwrite LDS before neighbors read it.
        float next_c = 1e10f;
        float next_d = 1e10f; // Track descent score locally for this step
        uint next_p = lp[ly][lx];
        
        // If we are sink, we hold our value
        if (lx == sink_lx && ly == sink_ly) {
            next_c = my_h;
            next_p = pack_coord(gx, gy);
        } else {
            // Start with current value as the baseline candidate
            // (Standard Bellman-Ford relaxation: min(old, neighbors))
            next_c = lc[ly][lx];
            // We don't store 'd' in LDS, so we can't perfectly compare against "old d".
            // However, by re-scanning neighbors, we will naturally re-find the best one.
            // So we treat 'next_c' as a candidate to beat.
        }

        // Scan Neighbors (Read from LDS)
        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue;
                
                int nx = lx + dx;
                int ny = ly + dy;
                
                if (nx >= 0 && nx < 16 && ny >= 0 && ny < 16) {
                    float dist_penalty = (dx != 0 && dy != 0) ? 1.41421f * 0.001f : 0.001f;
                    float neighbor_cost = lc[ny][nx]; // READ from shared current state
                    
                    float candidate_c = max(neighbor_cost + dist_penalty, my_h);
                    float candidate_d = neighbor_cost + dist_penalty;
                    
                    uint cand_p = pack_coord(get_group_id(0)*16 + nx, get_group_id(1)*16 + ny);

                    // Comparison Logic
                    if (candidate_c < next_c) {
                        next_c = candidate_c;
                        next_d = candidate_d;
                        next_p = cand_p;
                    } 
                    else if (candidate_c == next_c) {
                        // Tie-Breaking
                        if (candidate_d < next_d) {
                            next_d = candidate_d;
                            next_p = cand_p;
                        }
                    }
                }
            }
        }
        
        // SYNC 1: Ensure everyone has finished READING the current state
        barrier(CLK_LOCAL_MEM_FENCE);
        
        // WRITE: Update LDS with the new state computed in registers
        lc[ly][lx] = next_c;
        lp[ly][lx] = next_p;
        
        // SYNC 2: Ensure everyone has finished WRITING before next read phase
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (gx < width && gy < height) {
        cost_map[gy * width + gx] = lc[ly][lx];
        parent_map[gy * width + gx] = lp[ly][lx];
    }
}

// KERNEL 1: Calculate Steepest Descent
__kernel void calc_flow(
    __global const float* heightmap,
    __global uint* flow_map,   // Pointer to downstream neighbor
    __global int* basin_ids,   // Initial Basin IDs (-1 or self)
    const int W, const int H)
{
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    if (gx >= W || gy >= H) return;

    int idx = gy * W + gx;
    float my_h = heightmap[idx];
    
    float min_h = my_h;
    uint best_p = pack_coord(gx, gy); // Default: Flow to self (Sink)

    // Check 8 neighbors
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx == 0 && dy == 0) continue;
            int nx = gx + dx;
            int ny = gy + dy;
            
            if (nx >= 0 && nx < W && ny >= 0 && ny < H) {
                float n_h = heightmap[ny * W + nx];
                
                // Physics: Steepest drop = (HeightDiff) / Distance
                // Diagonal distance is 1.414, Cardinal is 1.0
                float dist = (dx!=0 && dy!=0) ? 1.41421f : 1.0f;
                float drop = (my_h - n_h) / dist;

                // We want to MAXIMIZE the drop (go downhill fastest)
                // Note: We only flow if n_h < my_h (drop > 0)
                if (n_h < min_h) { // Simple check first: is it lower?
                     // Use drop rate to break ties or prefer steepness?
                     // Standard watershed usually just takes absolute lowest neighbor.
                     // But taking 'Steepest Gradient' reduces kinks.
                     
                     // Let's store the neighbor that gives max gradient drop
                     // Implicitly, we track min_h for valid check, but optimized for gradient
                }
            }
        }
    }
    
    // REVISED SELECTION LOGIC FOR GRADIENT
    float max_gradient = -1.0f;
    
    for (int dy = -1; dy <= 1; dy++) {
        for (int dx = -1; dx <= 1; dx++) {
            if (dx == 0 && dy == 0) continue;
            int nx = gx + dx;
            int ny = gy + dy;
            if (nx >= 0 && nx < W && ny >= 0 && ny < H) {
                float n_h = heightmap[ny * W + nx];
                if (n_h < my_h) {
                    float dist = (dx!=0 && dy!=0) ? 1.41421f : 1.0f;
                    float gradient = (my_h - n_h) / dist;
                    
                    if (gradient > max_gradient) {
                        max_gradient = gradient;
                        best_p = pack_coord(nx, ny);
                    }
                }
            }
        }
    }

    flow_map[idx] = best_p;
    
    // If I point to myself, I am a sink. Get a unique ID (my index).
    if (best_p == pack_coord(gx, gy)) {
        basin_ids[idx] = idx;
    } else {
        basin_ids[idx] = -1; // Unknown yet
    }
}

// KERNEL 2: Pointer Jumping (Propagate Basin IDs up stream)
// In reality, we propagate the ID stored at the flow target BACK to the current pixel.
// Standard Pointer Jumping: B[i] = B[Flow[i]]
__kernel void propagate_basins(
    __global const uint* flow_map,
    __global int* basin_ids,
    __global int* changed_flag,
    const int W, const int H)
{
    int gx = get_global_id(0);
    int gy = get_global_id(1);
    int idx = gy * W + gx;
    if (gx >= W || gy >= H) return;

    int my_id = basin_ids[idx];
    
    // If I don't have an ID yet, or if my downstream neighbor has a 'better' (resolved) ID
    if (my_id == -1) {
        uint target_packed = flow_map[idx];
        int2 t = unpack_coord(target_packed);
        int target_idx = t.y * W + t.x;
        
        int target_id = basin_ids[target_idx];
        
        if (target_id != -1) {
            basin_ids[idx] = target_id;
            *changed_flag = 1;
        }
    } 
    // Pointer Jumping Optimization:
    // Even if I have an ID, check if my downstream points to something different 
    // (This helps flatten chains faster)
    else {
        uint target_packed = flow_map[idx];
        int2 t = unpack_coord(target_packed);
        int target_idx = t.y * W + t.x;
        int target_id = basin_ids[target_idx];
        
        if (target_id != -1 && target_id != my_id) {
             basin_ids[idx] = target_id;
             *changed_flag = 1;
        }
    }
}




// MurmurHash3 integer finalizer mix function for high quality randomness
uint hash_u32(uint h) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

float random_float(uint seed) {
    return (float)(hash_u32(seed) & 0x007FFFFF) / 8388607.0f; // float in [0,1]
}


// Read from buffer with clamp-to-edge
float read_tex(__global const float* data, int x, int y, int w, int h) {
    x = clamp(x, 0, w - 1);
    y = clamp(y, 0, h - 1);
    return data[y * w + x];
}


// ... (Helper functions remain the same) ...

// __kernel void accumulate_sdf(
//     __global const float* input_noise,
//     __global float* accumulator,
//     const int width,
//     const int height,
//     const float threshold,
//     const int radius)
// {
//     int gid_x = get_global_id(0);
//     int gid_y = get_global_id(1);

//     if (gid_x >= width || gid_y >= height) return;

//     // 1. Determine State
//     float center_val = read_tex(input_noise, gid_x, gid_y, width, height);
//     bool center_is_in = center_val > threshold;

//     // 2. Search neighbors (SDF)
//     float min_dist_sq = (float)(radius * radius * 2);
    
//     for (int dy = -radius; dy <= radius; dy++) {
//         for (int dx = -radius; dx <= radius; dx++) {
//             float dist_sq = (float)(dx*dx + dy*dy);
//             if (dist_sq >= min_dist_sq) continue;

//             int nx = gid_x + dx;
//             int ny = gid_y + dy;

//             float n_val = read_tex(input_noise, nx, ny, width, height);
//             bool n_is_in = n_val > threshold;

//             if (n_is_in != center_is_in) {
//                 if (dist_sq < min_dist_sq) min_dist_sq = dist_sq;
//             }
//         }
//     }

//     // 3. Process Distance
//     float dist = sqrt(min_dist_sq);
//     dist = max(0.0f, dist - 0.5f);
    
//     // Sign Flip: Inside (Land) is Positive for "Mountain" growth
//     // (Standard accumulation logic usually adds 'Land' height)
//     if (!center_is_in) dist = -dist; 
//     // Wait, let's align with the Shadertoy logic:
//     // "Inside" means closer to peak. "Outside" means closer to valley.
//     // The previous code had: if(center_is_in) dist = -dist; (Inside is negative)
//     // Then norm_dist = 0.5 - dist/... -> Inside maps to > 0.5.
    
//     // Let's stick to the mapping: 
//     // Land (Inside) -> Positive 0.5 to 1.0
//     // Sea (Outside) -> Positive 0.0 to 0.5
    
//     if (center_is_in) dist = -dist; // Inside is negative distance
    
//     // Map [-radius, +radius] to [0, 1]
//     float norm_dist = 0.5f - (dist / (float)(radius * 2.0f));
    
//     // Clamp
//     norm_dist = clamp(norm_dist, 0.0f, 1.0f);
    
//     // 4. SHARPNESS TRICK (The Fix)
//     // Apply smoothstep to soften the midpoint
//     norm_dist = smoothstep(0.0f, 1.0f, norm_dist);
    
//     // Apply Power Curve to sharpen peaks/valleys
//     // Squaring makes the 'bottom' (0.0) flatter and 'top' (1.0) pointier
//     // or vice versa depending on where the mass is.
//     // The Shadertoy used: val = val * val;
//     float contribution = norm_dist * norm_dist; 

//     // 5. Accumulate
//     int idx = gid_y * width + gid_x;
//     accumulator[idx] += contribution;
// }






// =================================================================================
// COMMON: Constants & Math Helpers (Ported from Shadertoy Common)
// =================================================================================

// #define PI 3.14159265359f

// // GLSL compatibility wrappers
// #define vec2 float2
// #define vec3 float3
// #define vec4 float4
// #define fract(x) ((x) - floor(x))
// #define mix(x, y, a) ((x) + ((y) - (x)) * (a))

// // Texture sampler helper (Bilinear interpolation wrapping)
// float sample_texture(__global const float* data, float2 uv, int2 res) {
//     float x = uv.x * res.x;
//     float y = uv.y * res.y;
    
//     // Wrap coords
//     float x_f = floor(x);
//     float y_f = floor(y);
//     int x0 = (int)x_f % res.x; if (x0 < 0) x0 += res.x;
//     int y0 = (int)y_f % res.y; if (y0 < 0) y0 += res.y;
//     int x1 = (x0 + 1) % res.x;
//     int y1 = (y0 + 1) % res.y;
    
//     float fx = x - x_f;
//     float fy = y - y_f;
    
//     float v00 = data[y0 * res.x + x0];
//     float v10 = data[y0 * res.x + x1];
//     float v01 = data[y1 * res.x + x0];
//     float v11 = data[y1 * res.x + x1];
    
//     return mix(mix(v00, v10, fx), mix(v01, v11, fx), fy);
// }

// // =================================================================================
// // NOISE FUNCTIONS (From Common Tab)
// // =================================================================================

// // Simplex 2D noise
// vec3 permute(vec3 x) { return fmod(((x*34.0f)+1.0f)*x, 289.0f); }

// float snoise(vec2 v){
//   const vec4 C = (vec4)(0.211324865405187f, 0.366025403784439f,
//            -0.577350269189626f, 0.024390243902439f);
//   vec2 i  = floor(v + dot(v, C.yy) );
//   vec2 x0 = v -   i + dot(i, C.xx);
//   vec2 i1;
//   i1 = (x0.x > x0.y) ? (vec2)(1.0f, 0.0f) : (vec2)(0.0f, 1.0f);
//   vec4 x12 = (vec4)(x0.x, x0.y, x0.x, x0.y) + (vec4)(C.x, C.x, C.z, C.z);
//   x12.xy -= i1;
//   i = fmod(i, 289.0f);
//   vec3 p = permute( permute( i.y + (vec3)(0.0f, i1.y, 1.0f ))
//   + i.x + (vec3)(0.0f, i1.x, 1.0f ));
//   vec3 m = max(0.5f - (vec3)(dot(x0,x0), dot(x12.xy,x12.xy),
//     dot(x12.zw,x12.zw)), 0.0f);
//   m = m*m ;
//   m = m*m ;
//   vec3 x = 2.0f * fract(p * C.www) - 1.0f;
//   vec3 h = fabs(x) - 0.5f;
//   vec3 ox = floor(x + 0.5f);
//   vec3 a0 = x - ox;
//   m *= 1.79284291400159f - 0.85373472095314f * ( a0*a0 + h*h );
//   vec3 g;
//   g.x  = a0.x  * x0.x  + h.x  * x0.y;
//   g.yz = a0.yz * x12.xz + h.yz * x12.yw;
//   return 0.5f + 0.5f*(130.0f * dot(m, g));
// }

// #define NUM_OCTAVES 8

// float fbm(vec2 x) {
//     float v = 0.0f;
//     float a = 0.5f;
//     vec2 shift = (vec2)(100.0f, 100.0f);
//     // Rotate to reduce axial bias
//     float c = cos(0.5f);
//     float s = sin(0.5f);
//     // mat2 rot = mat2(c, s, -s, c);
    
//     for (int i = 0; i < NUM_OCTAVES; ++i) {
//         v += a * snoise(x);
        
//         // Manual matrix mult for rotation
//         float nx = c * x.x + s * x.y;
//         float ny = -s * x.x + c * x.y;
//         x = (vec2)(nx, ny) * 2.0f + shift;
        
//         a *= 0.5f;
//     }
//     return v;
// }

// // =================================================================================
// // KERNEL: GENERATE INPUT NOISE (Simulates iChannel0 texture)
// // =================================================================================
// __kernel void generate_noise(
//     __global float* output,
//     int width,
//     int height,
//     float seed)
// {
//     int x = get_global_id(0);
//     int y = get_global_id(1);
//     if (x >= width || y >= height) return;
    
//     vec2 uv = (vec2)((float)x / width, (float)y / height);
    
//     // Scale up UV for FBM
//     float val = fbm(uv * 10.0f + (vec2)(seed, seed));
    
//     // Normalize roughly to 0..1
//     val = val * 0.5f + 0.5f;
    
//     output[y * width + x] = val;
// }

// // =================================================================================
// // KERNEL: BUFFER A (Erosion Accumulation)
// // =================================================================================
// // Helper from Buffer A
// float getOccupancy(__global const float* iChannel0, vec2 uv, int2 res) {
//     return sample_texture(iChannel0, uv, res);
// }

// bool isIn(__global const float* iChannel0, vec2 uv, int2 res, float threshold) {
//     return getOccupancy(iChannel0, uv, res) > threshold;
// }

// __kernel void buffer_a(
//     __global const float* iChannel0, // Input Noise
//     __global float* iChannel2,       // Accumulator (Self)
//     int width,
//     int height,
//     int iFrame)
// {
//     int x = get_global_id(0);
//     int y = get_global_id(1);
//     if (x >= width || y >= height) return;
    
//     int2 iResolution = (int2)(width, height);
//     vec2 fragCoord = (vec2)((float)x + 0.5f, (float)y + 0.5f);
    
//     // --- Port of Buffer A mainImage ---
    
//     // Compute Noise
//     vec3 noise = (vec3)(0.0f);
    
//     // Temporal Noise
//     vec3 temporalNoise = (vec3)((float)iFrame, (float)iFrame+1.0f, (float)iFrame+2.0f);
//     temporalNoise *= 1.618033f;
//     temporalNoise -= floor(temporalNoise);
//     noise += temporalNoise;
    
//     // Wrap values
//     noise -= floor(noise);
    
//     // Center offset
//     noise.xy -= 0.5f;
    
//     // Compute signed distance
//     float distanceToEdge;
//     {
//         vec2 samplingCenter = fragCoord + noise.xy;
//         vec2 samplingCenterUV = samplingCenter / (vec2)((float)width, (float)height);
        
//         const int iRange = 16; // Radius
//         const float range = (float)iRange;
//         const float maxSqrDist = range * range;
//         vec2 startPosition = samplingCenter;
        
//         bool fragIsIn = isIn(iChannel0, samplingCenterUV, iResolution, noise.z);
        
//         float squaredDistanceToEdge = maxSqrDist;
        
//         // Brute force search
//         for(int dx = -iRange; dx <= iRange; dx++) {
//             for(int dy = -iRange; dy <= iRange; dy++) {
//                 vec2 delta = (vec2)((float)dx, (float)dy);
//                 vec2 scanPosition = startPosition + delta;
//                 float scanDistanceSqr = dot(delta, delta);
                
//                 if(scanDistanceSqr >= maxSqrDist) continue;
//                 if(scanDistanceSqr >= squaredDistanceToEdge) continue;
                
//                 vec2 scanUV = scanPosition / (vec2)((float)width, (float)height);
//                 bool scanIsIn = isIn(iChannel0, scanUV, iResolution, noise.z);
                
//                 if (scanIsIn != fragIsIn) {
//                     squaredDistanceToEdge = scanDistanceSqr;
//                 }
//             }
//         }
        
//         distanceToEdge = sqrt(squaredDistanceToEdge);
//         distanceToEdge -= 0.5f;
//         distanceToEdge = fragIsIn ? -distanceToEdge : distanceToEdge;
        
//         // Normalize
//         distanceToEdge /= range * 2.0f;
//         distanceToEdge = 0.5f - distanceToEdge;
//     }
    
//     // Smoothstep
//     distanceToEdge = smoothstep(0.0f, 1.0f, distanceToEdge);
    
//     // Accumulate
//     // iChannel2 stores the sum. We read, add, and write back.
//     int idx = y * width + x;
//     float oldVal = iChannel2[idx];
    
//     // The shader outputs vec4(d,d,d,1). We just store float d.
//     // 'oldColor' accumulation logic: fragColor += oldColor;
//     iChannel2[idx] = oldVal + distanceToEdge;
// }

// // =================================================================================
// // KERNEL: BUFFER B (Post-Processing / Perlin Composition)
// // =================================================================================
// // Function 'perlin' from Buffer B
// float buffer_b_perlin(__global const float* iChannel0, vec2 uv, int2 res) {
//     uv += (vec2)(8.2813f, 1.42114f);
//     uv /= 4.0f;
//     vec2 occ = (vec2)(0.0f);
//     float a = 1.0f;
    
//     for(int i = 0; i < 7; i++) {
//         float val = sample_texture(iChannel0, uv, res);
//         occ += (vec2)(val, 1.0f) * a;
//         uv *= 0.5f;
//         a *= 2.0f;
//     }
    
//     float v = occ.x / occ.y;
    
//     // Contrast
//     v = v * 2.0f - 1.0f;
//     v = tanh(v * 2.0f);
//     v = v * 0.5f + 0.5f;
    
//     v *= v;
    
//     return v;
// }

// __kernel void buffer_b(
//     __global const float* iChannel0, // This is the Accumulated Buffer A
//     __global float* output,          // Heightmap Output
//     int width,
//     int height,
//     int total_frames)
// {
//     int x = get_global_id(0);
//     int y = get_global_id(1);
//     if (x >= width || y >= height) return;
    
//     int2 res = (int2)(width, height);
//     vec2 uv = (vec2)((float)x / width, (float)y / height);
    
//     // Note: The accumulator (Buffer A) contains the SUM. 
//     // Buffer B 'perlin' function expects sampled values to be roughly 0..1 range?
//     // In the shader Buffer B, 'texture(iChannel0, uv)' returns the averaged color 
//     // IF Buffer A was a float32 buffer and not cleared every frame.
//     // However, Shadertoy buffers are persistent.
//     // The perlin function sums octaves.
    
//     // WE NEED TO NORMALIZE the accumulator first, because 'sample_texture' reads raw sum.
//     // Or we can do it on the fly.
//     // Since we can't easily normalize inside sample_texture without passing factor,
//     // let's assume we pass the Normalized Buffer A as iChannel0 here.
    
//     // Logic:
//     float height_val = buffer_b_perlin(iChannel0, uv, res);
    
//     output[y * width + x] = height_val;
// }


// // =================================================================================
// // KERNEL: RENDER IMAGE (The Lighting & Coloring Shader)
// // =================================================================================

// // Helper to match Shadertoy's Skybox (used for ambient lighting)
// float3 skybox_blurry(float3 dir) {
//     float gradient = dir.y * 0.5f + 0.5f;
//     gradient = 1.0f - gradient;
//     gradient *= gradient * gradient;
//     gradient = 1.0f - gradient;
    
//     float3 gradient3 = pow((float3)(gradient), (float3)(8.0f, 1.0f, 1.0f));
//     gradient3 = (float3)(1.0f) - gradient3;
//     gradient3 = pow(gradient3, (float3)(0.4f, 0.5f, 4.0f));
//     gradient3 = (float3)(1.0f) - gradient3;
    
//     return mix((float3)(0.99f, 0.99f, 0.99f), (float3)(0.01f, 0.02f, 0.2f), gradient3) * 2.0f;
// }

// __kernel void render_map(
//     __global const float* heightmap,
//     __global float4* image_out, // Output RGBA image
//     int width,
//     int height)
// {
//     int x = get_global_id(0);
//     int y = get_global_id(1);
//     if (x >= width || y >= height) return;
    
//     int2 res = (int2)(width, height);
//     float2 uv = (float2)((float)x / width, (float)y / height);
    
//     // --- GEOMETRY CALCULATION (From Image Tab) ---
    
//     // Height Scaling (matches 'maxHeight' in shader)
//     // The shader uses 500.0 / Resolution.x. For 512, that's ~1.0.
//     float height_scale = 1.0f; 
    
//     // Precision for derivatives (matches 'precis')
//     float prec = 1.5f / width;
    
//     // Sample Neighbors for Normal/Curvature
//     // We treat the heightmap texture as the terrain function
//     float h_c = sample_texture(heightmap, uv, res);
//     float h_r = sample_texture(heightmap, uv + (float2)(prec, 0.0f), res);
//     float h_l = sample_texture(heightmap, uv - (float2)(prec, 0.0f), res);
//     float h_t = sample_texture(heightmap, uv + (float2)(0.0f, prec), res);
//     float h_b = sample_texture(heightmap, uv - (float2)(0.0f, prec), res);
    
//     float3 posC = (float3)(0.0f, h_c * height_scale, 0.0f);
//     float3 posR = (float3)(prec, h_r * height_scale, 0.0f);
//     float3 posL = (float3)(-prec, h_l * height_scale, 0.0f);
//     float3 posT = (float3)(0.0f, h_t * height_scale, prec);
//     float3 posB = (float3)(0.0f, h_b * height_scale, -prec);
    
//     float3 dx = posR - posL;
//     float3 dy = posT - posB;
    
//     // Normal
//     float3 normal = normalize(cross(dy, dx)); // Note: OpenCL cross order might differ, shader used dx, dy -> normal up?
//     // GLSL: cross(dx, dy). If x is right, z is up/forward?
//     // Let's ensure normal points UP. h_c is Y.
//     // If normal.y is negative, flip it.
//     if (normal.y < 0) normal = -normal;
    
//     // Curvature
//     float curveX = -dot(posC + posC - posR - posL, normal);
//     float curveY = -dot(posC + posC - posT - posB, normal);
//     float2 curve2 = (float2)(curveX, curveY) / prec;
    
//     // --- SHADING LOGIC (From castRayTerrain2) ---
    
//     float2 posCurve2 = max((float2)(0.0f), curve2);
//     float2 negCurve2 = -min((float2)(0.0f), curve2);
    
//     float posCurve = posCurve2.x + posCurve2.y;
//     posCurve = posCurve / (0.2f + posCurve);
    
//     float negCurve = negCurve2.x + negCurve2.y;
//     negCurve = negCurve / (0.2f + negCurve);
    
//     // Lighting Vectors
//     float3 lightDir = normalize((float3)(0.0f, 3.0f, 5.0f));
//     float3 lightColor = max(0.0f, dot(normal, lightDir)) * (float3)(1.0f, 0.75f, 0.5f) * 2.0f;
    
//     float3 ambientDir = (float3)(0.0f, 1.0f, 0.0f);
//     float ambient = dot(normal, ambientDir) * 0.5f + 0.5f;
//     ambient *= ambient;
//     ambient *= posCurve * 0.5f + 0.5f; // curve influence
//     float3 ambientColor = ambient * skybox_blurry(normal);
    
//     float3 totalDiffuse = ambientColor + lightColor;
    
//     // --- MATERIAL / ALBEDO (The Colors!) ---
    
//     float3 albedo = (float3)(1.0f);
    
//     float slopeFactor;
//     float heightFactor;
    
//     {
//         heightFactor = h_c; // Normalized height 0..1
//         heightFactor = smoothstep(0.0f, 1.0f, heightFactor);
//         heightFactor = smoothstep(0.0f, 1.0f, heightFactor);
        
//         slopeFactor = fabs(normal.y);
//         slopeFactor *= slopeFactor * slopeFactor * slopeFactor; // ^4
//         slopeFactor = 1.0f - slopeFactor;
//         slopeFactor *= slopeFactor * slopeFactor; // ^2
//         slopeFactor = 1.0f - slopeFactor;
        
//         // Palettes
//         float3 flatLow   = (float3)(0.1f, 0.4f, 0.05f);  // Greenish
//         float3 flatHigh  = (float3)(0.7f, 0.7f, 0.5f);   // Gray/Brown
//         float3 slopeLow  = (float3)(0.4f, 0.7f, 0.2f);   // Steep Green
//         float3 slopeHigh = (float3)(0.3f, 0.2f, 0.05f);  // Dark Rock
        
//         float3 flatColor = mix(flatLow, flatHigh, heightFactor);
//         float3 slopeColor = mix(slopeLow, slopeHigh, heightFactor);
        
//         albedo = mix(slopeColor, flatColor, slopeFactor);
//     }
    
//     // Add "Sediment" / Curve color
//     albedo = mix(albedo, (albedo * 1.0f + (float3)(0.75f)) * (float3)(0.7f, 0.65f, 0.3f), posCurve);
    
//     // --- RIVERS (The Blue Lines) ---
//     float rivers = max(0.0f, negCurve - posCurve);
//     rivers = pow(rivers, pow(2.0f, mix(0.5f, -0.7f, slopeFactor) + mix(-0.7f, 0.7f, heightFactor)));
//     rivers = smoothstep(0.0f, 1.0f, rivers);
//     rivers = smoothstep(0.0f, 1.0f, rivers);
    
//     albedo = mix(albedo, (float3)(0.1f, 0.2f, 0.5f), rivers);
    
//     // Specular / Final Compose
//     albedo *= albedo; // Gamma correctionish
    
//     // Top down view approximation (skip specular eye calculation for map view, or assume top-down view)
//     // To match image, let's just use the diffuse lighting
    
//     float3 final_col = albedo * totalDiffuse;
    
//     // Gamma decode for display
//     final_col = pow(final_col, (float3)(1.0f/2.2f));
    
//     image_out[y * width + x] = (float4)(final_col, 1.0f);
// }





// =================================================================================
// COMMON: Constants & Math Helpers
// =================================================================================

#define PI 3.14159265359f

// GLSL compatibility wrappers
#define vec2 float2
#define vec3 float3
#define vec4 float4
#define fract(x) ((x) - floor(x))
#define mix(x, y, a) ((x) + ((y) - (x)) * (a))

// CRITICAL FIX: GLSL style modulo for noise (handles negatives correctly)
#define mod(x, y) ((x) - (y) * floor((x) / (y)))

// Texture sampler helper (Bilinear interpolation wrapping)
float sample_texture(__global const float* data, float2 uv, int2 res) {
    float x = uv.x * res.x;
    float y = uv.y * res.y;
    
    // Wrap coords
    float x_f = floor(x);
    float y_f = floor(y);
    int x0 = (int)mod(x_f, (float)res.x);
    int y0 = (int)mod(y_f, (float)res.y);
    int x1 = (x0 + 1) % res.x;
    int y1 = (y0 + 1) % res.y;
    
    float fx = x - x_f;
    float fy = y - y_f;
    
    float v00 = data[y0 * res.x + x0];
    float v10 = data[y0 * res.x + x1];
    float v01 = data[y1 * res.x + x0];
    float v11 = data[y1 * res.x + x1];
    
    return mix(mix(v00, v10, fx), mix(v01, v11, fx), fy);
}

// =================================================================================
// NOISE FUNCTIONS (From Common Tab)
// =================================================================================

// Simplex 2D noise
vec3 permute(vec3 x) { return mod(((x*34.0f)+1.0f)*x, 289.0f); }

float snoise(vec2 v){
  const vec4 C = (vec4)(0.211324865405187f, 0.366025403784439f,
           -0.577350269189626f, 0.024390243902439f);
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);
  vec2 i1;
  i1 = (x0.x > x0.y) ? (vec2)(1.0f, 0.0f) : (vec2)(0.0f, 1.0f);
  vec4 x12 = (vec4)(x0.x, x0.y, x0.x, x0.y) + (vec4)(C.x, C.x, C.z, C.z);
  x12.xy -= i1;
  i = mod(i, 289.0f); // Uses fixed mod macro
  vec3 p = permute( permute( i.y + (vec3)(0.0f, i1.y, 1.0f ))
  + i.x + (vec3)(0.0f, i1.x, 1.0f ));
  vec3 m = max(0.5f - (vec3)(dot(x0,x0), dot(x12.xy,x12.xy),
    dot(x12.zw,x12.zw)), 0.0f);
  m = m*m ;
  m = m*m ;
  vec3 x = 2.0f * fract(p * C.www) - 1.0f;
  vec3 h = fabs(x) - 0.5f;
  vec3 ox = floor(x + 0.5f);
  vec3 a0 = x - ox;
  m *= 1.79284291400159f - 0.85373472095314f * ( a0*a0 + h*h );
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 0.5f + 0.5f*(130.0f * dot(m, g));
}

#define NUM_OCTAVES 8

float fbm(vec2 x) {
    float v = 0.0f;
    float a = 0.5f;
    vec2 shift = (vec2)(100.0f, 100.0f);
    float c = cos(0.5f);
    float s = sin(0.5f);
    
    for (int i = 0; i < NUM_OCTAVES; ++i) {
        v += a * snoise(x);
        float nx = c * x.x + s * x.y;
        float ny = -s * x.x + c * x.y;
        x = (vec2)(nx, ny) * 2.0f + shift;
        a *= 0.5f;
    }
    return v;
}

// =================================================================================
// KERNEL: GENERATE INPUT NOISE
// =================================================================================
__kernel void generate_noise(
    __global float* output,
    int width,
    int height,
    float seed)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    if (x >= width || y >= height) return;
    
    // Match Shadertoy aspect ratio scaling logic roughly if needed, 
    // but straight UV is fine for filling texture.
    vec2 uv = (vec2)((float)x / width, (float)y / height);
    
    // Scale up UV for FBM so we see pattern
    float val = fbm(uv * 5.0f + (vec2)(seed, seed));
    
    // Normalize roughly to 0..1
    val = val * 0.5f + 0.5f;
    
    output[y * width + x] = val;
}

// =================================================================================
// KERNEL: BUFFER A (Erosion Accumulation)
// =================================================================================
float getOccupancy(__global const float* iChannel0, vec2 uv, int2 res) {
    return sample_texture(iChannel0, uv, res);
}

bool isIn(__global const float* iChannel0, vec2 uv, int2 res, float threshold) {
    return getOccupancy(iChannel0, uv, res) > threshold;
}

__kernel void buffer_a(
    __global const float* iChannel0, 
    __global float* iChannel2,       
    int width,
    int height,
    int iFrame)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    if (x >= width || y >= height) return;
    
    int2 iResolution = (int2)(width, height);
    vec2 fragCoord = (vec2)((float)x + 0.5f, (float)y + 0.5f);
    
    // Temporal Noise (Threshold Jitter)
    vec3 noise = (vec3)(0.0f);
    vec3 temporalNoise = (vec3)((float)iFrame, (float)iFrame+1.0f, (float)iFrame+2.0f);
    temporalNoise *= 1.618033f;
    temporalNoise -= floor(temporalNoise);
    noise += temporalNoise;
    noise -= floor(noise);
    noise.xy -= 0.5f;
    
    float distanceToEdge;
    {
        vec2 samplingCenter = fragCoord + noise.xy;
        vec2 samplingCenterUV = samplingCenter / (vec2)((float)width, (float)height);
        
        const int iRange = 16;
        const float range = (float)iRange;
        const float maxSqrDist = range * range;
        vec2 startPosition = samplingCenter;
        
        bool fragIsIn = isIn(iChannel0, samplingCenterUV, iResolution, noise.z);
        float squaredDistanceToEdge = maxSqrDist;
        
        // Loop range: -16 to 16
        for(int dx = -iRange; dx <= iRange; dx++) {
            for(int dy = -iRange; dy <= iRange; dy++) {
                vec2 delta = (vec2)((float)dx, (float)dy);
                // Optimization: Pre-check square dist
                float scanDistanceSqr = dot(delta, delta);
                if(scanDistanceSqr >= maxSqrDist) continue;
                if(scanDistanceSqr >= squaredDistanceToEdge) continue;

                vec2 scanPosition = startPosition + delta;
                vec2 scanUV = scanPosition / (vec2)((float)width, (float)height);
                bool scanIsIn = isIn(iChannel0, scanUV, iResolution, noise.z);
                
                if (scanIsIn != fragIsIn) {
                    squaredDistanceToEdge = scanDistanceSqr;
                }
            }
        }
        
        distanceToEdge = sqrt(squaredDistanceToEdge);
        distanceToEdge -= 0.5f; // Boundary correction
        distanceToEdge = fragIsIn ? -distanceToEdge : distanceToEdge;
        
        distanceToEdge /= range * 2.0f;
        distanceToEdge = 0.5f - distanceToEdge;
    }
    
    distanceToEdge = smoothstep(0.0f, 1.0f, distanceToEdge);
    
    // Accumulate
    int idx = y * width + x;
    iChannel2[idx] += distanceToEdge;
}

// =================================================================================
// KERNEL: BUFFER B (Post-Processing)
// =================================================================================
float buffer_b_perlin(__global const float* iChannel0, vec2 uv, int2 res) {
    // Offset UVs to sample interesting part of noise
    uv += (vec2)(8.2813f, 1.42114f);
    uv /= 4.0f;
    
    vec2 occ = (vec2)(0.0f);
    float a = 1.0f;
    
    // FBM Stacking
    for(int i = 0; i < 7; i++) {
        float val = sample_texture(iChannel0, uv, res);
        occ += (vec2)(val, 1.0f) * a;
        uv *= 0.5f;
        a *= 2.0f;
    }
    
    float v = occ.x / occ.y;
    
    // Contrast Curve
    v = v * 2.0f - 1.0f;
    v = tanh(v * 2.0f);
    v = v * 0.5f + 0.5f;
    v *= v; // Peak sharpening
    
    return v;
}

__kernel void buffer_b(
    __global const float* iChannel0,
    __global float* output,
    int width,
    int height)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    if (x >= width || y >= height) return;
    
    int2 res = (int2)(width, height);
    vec2 uv = (vec2)((float)x / width, (float)y / height);
    
    float height_val = buffer_b_perlin(iChannel0, uv, res);
    output[y * width + x] = height_val;
}

// =================================================================================
// KERNEL: RENDER (Lighting & Color)
// =================================================================================
float3 skybox_blurry(float3 dir) {
    float gradient = dir.y * 0.5f + 0.5f;
    gradient = 1.0f - gradient;
    gradient *= gradient * gradient;
    gradient = 1.0f - gradient;
    
    float3 gradient3 = pow((float3)(gradient), (float3)(8.0f, 1.0f, 1.0f));
    gradient3 = (float3)(1.0f) - gradient3;
    gradient3 = pow(gradient3, (float3)(0.4f, 0.5f, 4.0f));
    gradient3 = (float3)(1.0f) - gradient3;
    
    return mix((float3)(0.99f, 0.99f, 0.99f), (float3)(0.01f, 0.02f, 0.2f), gradient3) * 2.0f;
}

__kernel void render_map(
    __global const float* heightmap,
    __global float4* image_out,
    int width,
    int height)
{
    int x = get_global_id(0);
    int y = get_global_id(1);
    if (x >= width || y >= height) return;
    
    int2 res = (int2)(width, height);
    float2 uv = (float2)((float)x / width, (float)y / height);
    
    // Render Params
    float height_scale = 1.0f; 
    float prec = 1.5f / width;
    
    // Derivatives
    float h_c = sample_texture(heightmap, uv, res);
    float h_r = sample_texture(heightmap, uv + (float2)(prec, 0.0f), res);
    float h_l = sample_texture(heightmap, uv - (float2)(prec, 0.0f), res);
    float h_t = sample_texture(heightmap, uv + (float2)(0.0f, prec), res);
    float h_b = sample_texture(heightmap, uv - (float2)(0.0f, prec), res);
    
    float3 posC = (float3)(0.0f, h_c * height_scale, 0.0f);
    float3 posR = (float3)(prec, h_r * height_scale, 0.0f);
    float3 posL = (float3)(-prec, h_l * height_scale, 0.0f);
    float3 posT = (float3)(0.0f, h_t * height_scale, prec);
    float3 posB = (float3)(0.0f, h_b * height_scale, -prec);
    
    float3 dx = posR - posL;
    float3 dy = posT - posB;
    
    float3 normal = normalize(cross(dy, dx));
    if (normal.y < 0) normal = -normal;
    
    // Curvature
    float curveX = -dot(posC + posC - posR - posL, normal);
    float curveY = -dot(posC + posC - posT - posB, normal);
    float2 curve2 = (float2)(curveX, curveY) / prec;
    
    float2 posCurve2 = max((float2)(0.0f), curve2);
    float2 negCurve2 = -min((float2)(0.0f), curve2);
    float posCurve = (posCurve2.x + posCurve2.y) / (0.2f + posCurve2.x + posCurve2.y);
    float negCurve = (negCurve2.x + negCurve2.y) / (0.2f + negCurve2.x + negCurve2.y);

    // Lighting
    float3 lightDir = normalize((float3)(0.0f, 3.0f, 5.0f));
    float3 lightColor = max(0.0f, dot(normal, lightDir)) * (float3)(1.0f, 0.75f, 0.5f) * 2.0f;
    
    float ambient = dot(normal, (float3)(0.0f, 1.0f, 0.0f)) * 0.5f + 0.5f;
    ambient *= ambient * (posCurve * 0.5f + 0.5f);
    float3 ambientColor = ambient * skybox_blurry(normal);
    
    // Material Colors
    float3 albedo;
    float slopeFactor = fabs(normal.y);
    slopeFactor = pow(slopeFactor, 4.0f);
    slopeFactor = 1.0f - slopeFactor;
    slopeFactor = pow(slopeFactor, 2.0f);
    slopeFactor = 1.0f - slopeFactor;
    
    float heightFactor = smoothstep(0.0f, 1.0f, h_c);
    heightFactor = smoothstep(0.0f, 1.0f, heightFactor);
    
    float3 flatLow = (float3)(0.1f, 0.4f, 0.05f);
    float3 flatHigh = (float3)(0.7f, 0.7f, 0.5f);
    float3 slopeLow = (float3)(0.4f, 0.7f, 0.2f);
    float3 slopeHigh = (float3)(0.3f, 0.2f, 0.05f);
    
    float3 flatColor = mix(flatLow, flatHigh, heightFactor);
    float3 slopeColor = mix(slopeLow, slopeHigh, heightFactor);
    albedo = mix(slopeColor, flatColor, slopeFactor);
    
    // Sediment
    albedo = mix(albedo, (albedo + 0.75f) * (float3)(0.7f, 0.65f, 0.3f), posCurve);
    
    // Rivers
    float rivers = max(0.0f, negCurve - posCurve);
    rivers = pow(rivers, pow(2.0f, mix(0.5f, -0.7f, slopeFactor) + mix(-0.7f, 0.7f, heightFactor)));
    rivers = smoothstep(0.0f, 1.0f, rivers);
    rivers = smoothstep(0.0f, 1.0f, rivers);
    albedo = mix(albedo, (float3)(0.1f, 0.2f, 0.5f), rivers);
    
    albedo *= albedo;
    
    float3 final_col = albedo * (ambientColor + lightColor);
    final_col = pow(final_col, (float3)(1.0f/2.2f));
    
    image_out[y * width + x] = (float4)(final_col, 1.0f);
}