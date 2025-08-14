__kernel void nbody_sim(__global float4* positions, __global float4* velocities, const float dt, const int particle_count) {
    int i = get_global_id(0);
    if(i==0)printf("GPU nbody_sim() dt %8.2e  particle_count %i \n", dt, particle_count);
    if (i >= particle_count) return;
    float3 my_pos = positions[i].xyz;
    float3 my_vel = velocities[i].xyz;
    float3 force = (float3)(0.0f, 0.0f, 0.0f);
    float G = 1.0f;
    for (int j = 0; j < particle_count; j++) {
        if (i == j) continue;
        float3 other_pos = positions[j].xyz;
        float3 diff = other_pos - my_pos;
        float dist_sq = dot(diff, diff) + 0.1f;
        float dist = sqrt(dist_sq);
        float3 direction = diff / dist;
        force += direction / dist_sq;
    }
    my_vel += force * G * dt;
    my_pos += my_vel * dt;
    positions[i].xyz = my_pos;
    velocities[i].xyz = my_vel;
}
