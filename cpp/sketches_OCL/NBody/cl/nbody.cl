
#define R2SAFE 1.0e-2
#define F2MAX  10.0

float2 force_LJ( float2 p1, float2 p2 ){
    float2 d   = p2 - p1;
    float  ir2 = 1/( d.x*d.x + d.y*d.y + R2SAFE );
    float  ir6 = ir2 * ir2 * ir2;
    float  fr  =  ir6 - ir6*ir6;
    return d * fr;
}

__kernel void NBody_force_naive(
    unsigned int num,
    __global float2* pos, 
    __global float2* force
){
    const int i = get_global_id (0);
    float2 p = pos[i];
    float2 f = (float2) (0.0f, 0.0f);
    for (int j=0; j < num; j++){
        float2 pj = pos[j];
        f += force_LJ( p, pj );
        //f += pj;
    }
    force[i] = f;
    //force[i] = (float2) (i, tid);
    //force[i] = (float2) (i, num);
    //force[i] = pos[i];
}

__kernel void NBody_force(
    unsigned int num,
    __global float2* pos, 
    __global float2* force
){
    const int i         = get_global_id (0);
    const int tid       = get_local_id  (0);
    const int block_dim = get_local_size(0);
    //__local  float2 pos_shared[block_dim];
    //__local  float2 pos_shared[get_local_size(0)];
    __local  float2 pos_shared[256];
    float2 p = pos[i];
    float2 f = (float2) (0.0f, 0.0f);
    int tile=0;
    //int body, j;
    for (int body=0; body <num; body += block_dim, tile++){
        pos_shared[tid] = pos[tile * block_dim + tid];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j < block_dim; j++){
            float2 pj = pos_shared[j];
           f += force_LJ( p, pj );
           //f += pj;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    force[i] = f;
    //force[i] = (float2) (i, tid);
    //force[i] = pos[i];
}


__kernel void NBody_step(
    int   num,
    __global float2* pos, 
    __global float2* vel, 
    float dt,    
    float damp, 
    int   nsteps
){
    //get our index in the array
    unsigned int i         = get_global_id(0);
    unsigned int tid       = get_local_id(0);
    unsigned int block_dim = get_local_size(0);

    for (int iter=0; iter<nsteps; iter++){

        // evaluate interaction forces
        __local  float2 pos_shared[256];
        float2 p = pos[i];
        float2 f = (float2) (0.0f, 0.0f);
        int tile=0;
        //int body, j;
        for (int body=0; body <num; body += block_dim, tile++){
            pos_shared[tid] = pos[tile * block_dim + tid];
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int j=0; j < block_dim; j++){
                float2 pj = pos_shared[j];
               f += force_LJ( p, pj );
               //f += pj;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        //force[i] = f;
        
        float2 v = vel[i];
        
        v  = v*damp +  f * dt;
        p += v * dt;

        pos[i] = p;
        vel[i] = v;

        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}


/*
__kernel void NBody_step(
    __global float4* pos, 
    __global float4* vel, 
    __local  float4* pos_shared, 
    float dt, 
    unsigned int num
){
    //get our index in the array
    unsigned int i         = get_global_id(0);
    unsigned int tid       = get_local_id(0);
    unsigned int block_dim = get_local_size(0);

    for (int iter=0; iter<10; iter++){
        float4 p = pos[i];
        float4 v = vel[i];

        float4 acc = (float4) (0.0f, 0.0f, 0.0f, 0.0f);

        int tile=0;
        int body, j;
        for (body=0; body <num; body += block_dim, tile++){
            pos_shared[tid] = pos_gen[tile * block_dim + tid];
            barrier(CLK_LOCAL_MEM_FENCE);
            for (j=0; j < block_dim; j++){
                float4 pj = pos_shared[j];
                float4 r =  pj - p;
                float dist = sqrt(r.x * r.x + r.y * r.y + r.z * r.z);
                dist += 0.1;
                float inv_dist = 1.0f / (dist);
                float inv_dist_cubed = inv_dist * inv_dist * inv_dist;
                float s = pj.w * inv_dist_cubed;
                acc -= r * s;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        v -= acc * dt;
        p += v   * dt;

        pos[i] = p;
        vel[i] = v;

        barrier(CLK_GLOBAL_MEM_FENCE);
    }

}

*/

