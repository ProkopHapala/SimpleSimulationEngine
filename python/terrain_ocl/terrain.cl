// Pack (x, y) into uint32
uint pack_coord(int x, int y) {
    return (uint)x | ((uint)y << 16);
}

// Local Relaxation Kernel
// Process 16x16 tiles. Each tile finds its sink and calculates internal flow.
__kernel void solve_tiles(
    __global const float* heightmap,
    __global uint* parent_map,  // Stores packed (x, y) of next hop
    __global float* cost_map,    // Cumulative cost from sink
    __global uint* tile_sinks,   // Output: Sink coord for this tile
    const int width,
    const int height,
    const float K)               // Quadratic potential weight
{
    int lx = get_local_id(0);
    int ly = get_local_id(1);
    int gx = get_group_id(0) * 16 + lx;
    int gy = get_group_id(1) * 16 + ly;
    int tile_id = get_group_id(1) * (width/16) + get_group_id(0);

    __local float local_h[16][16];
    __local float local_cost[16][16];
    __local uint local_parent[16][16];

    // 1. Load data to Local Memory
    if (gx < width && gy < height) {
        local_h[ly][lx] = heightmap[gy * width + gx];
    } else {
        local_h[ly][lx] = 1e10f;
    }
    
    // 2. Find Sink (Reduction)
    // Add quadratic potential to height to find the local 'Capital'
    float centerX = 7.5f;
    float centerY = 7.5f;
    float distSq = (lx - centerX)*(lx - centerX) + (ly - centerY)*(ly - centerY);
    float potential_h = local_h[ly][lx] + K * distSq;
    
    local_cost[ly][lx] = potential_h;
    local_parent[ly][lx] = pack_coord(gx, gy); 
    barrier(CLK_LOCAL_MEM_FENCE);

    // Parallel reduction to find min potential in tile
    __local int2 min_pos;
    if (lx == 0 && ly == 0) {
        float min_val = 1e10f;
        for(int j=0; j<16; j++) {
            for(int i=0; i<16; i++) {
                if(local_cost[j][i] < min_val) {
                    min_val = local_cost[j][i];
                    min_pos = (int2)(i, j);
                }
            }
        }
        tile_sinks[tile_id] = pack_coord(get_group_id(0) * 16 + min_pos.x, get_group_id(1) * 16 + min_pos.y);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 3. Local Relaxation (Bellman-Ford style)
    // Initialize
    int sink_lx = min_pos.x;
    int sink_ly = min_pos.y;
    local_cost[ly][lx] = 1e10f;
    local_parent[ly][lx] = pack_coord(gx, gy);
    if (lx == sink_lx && ly == sink_ly) {
        local_cost[ly][lx] = local_h[ly][lx];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Iterate to flood the 16x16 tile
    // Max distance is 32. 
    for (int it = 0; it < 32; it++) {
        float my_h = local_h[ly][lx];
        float best_c = local_cost[ly][lx];
        uint best_p = local_parent[ly][lx];

        for (int dy = -1; dy <= 1; dy++) {
            for (int dx = -1; dx <= 1; dx++) {
                if (dx == 0 && dy == 0) continue; // skip self
                int nx = lx + dx;
                int ny = ly + dy;
                if (nx >= 0 && nx < 16 && ny >= 0 && ny < 16) {
                    // Cost function: max height on path (Flooding)
                    float candidate = max(local_cost[ny][nx], my_h);
                    uint cand_p = pack_coord(get_group_id(0)*16 + nx, get_group_id(1)*16 + ny);
                    if (candidate < best_c || (candidate == best_c && cand_p < best_p)) {
                        best_c = candidate;
                        best_p = cand_p;
                    }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        local_cost[ly][lx] = best_c;
        local_parent[ly][lx] = best_p;
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // 4. Export results
    if (gx < width && gy < height) {
        cost_map[gy * width + gx] = local_cost[ly][lx];
        parent_map[gy * width + gx] = local_parent[ly][lx];
    }
}