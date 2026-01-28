// Pack (x, y) into uint32
uint pack_coord(int x, int y) {
    return (uint)x | ((uint)y << 16);
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