
#define K_bend     5.0f 
#define param_k  -10.0f
#define param_l0   1.0f
#define param_kX  -5.0f
#define param_l0X  1.41421356237f

float3 stickForce( float3 p1, float3 p2, float k, float l0 ){
    float3 d = p2 - p1;
    float  l = sqrt(dot(d,d));
    d*=  k*(l-l0)/(l+1e-8f);
    return d;
}

float3 harmonicForce( float3 p1, float3 p2, float k ){
    float3 d = p2 - p1;
    return d * k;
}

__kernel void cloth_force(
    __global float4* pos, 
    __global float4* force
){    

    const int i = get_global_id(1)*get_global_size(0) + get_global_id(0);
    //printf("global_id (%li,%li) %i \n", get_global_id(0), get_global_id(1), i );
    float3 p = pos[i].xyz;
    float3 f = (float3) (0.0f,0.0f,0.0f);
        
    if( get_global_id(1)>0                   ){ f += stickForce( pos[i-get_global_size(0)].xyz, p, param_k, param_l0 ); }
    if( get_global_id(1)<get_global_size(1)-1){ f += stickForce( pos[i+get_global_size(0)].xyz, p, param_k, param_l0 ); }
    if( get_global_id(0)>0                   ){ f += stickForce( pos[i-1                 ].xyz, p, param_k, param_l0 ); }
    if( get_global_id(0)<get_global_size(0)-1){ f += stickForce( pos[i+1                 ].xyz, p, param_k, param_l0 ); }
    f.y += -0.05f;
    
    force[i] = (float4)(f,0.0f);
    //force[i] = 1545.0f;
}

__kernel void cloth_force_bend(
    __global float4* pos, 
    __global float4* force
){    

    const int i = get_global_id(1)*get_global_size(0) + get_global_id(0);
    //printf("global_id (%li,%li) %i \n", get_global_id(0), get_global_id(1), i );
    float3 p = pos[i].xyz;
    float3 f = (float3) (0.0f,0.0f,0.0f);
    float3 c;    
    if( get_global_id(1)>0                   ){ float3 pj =pos[i-get_global_size(0)].xyz; f += stickForce( pj, p, param_k, param_l0 ); c=pj; }
    if( get_global_id(1)<get_global_size(1)-1){ float3 pj =pos[i+get_global_size(0)].xyz; f += stickForce( pj, p, param_k, param_l0 ); 
        if( get_global_id(1)>0               ){ f+=harmonicForce( p, (c+pj)*0.5f, K_bend );}   }
    if( get_global_id(0)>0                   ){ float3 pj =pos[i-1                 ].xyz; f += stickForce( pj, p, param_k, param_l0 ); c=pj; }
    if( get_global_id(0)<get_global_size(0)-1){ float3 pj =pos[i+1                 ].xyz; f += stickForce( pj, p, param_k, param_l0 ); 
        if( get_global_id(0)>0               ){ f+=harmonicForce( p, (c+pj)*0.5f, K_bend ); }   }
    f.y += -0.05f;
    
    force[i] = (float4)(f,0.0f);
    //force[i] = 1545.0f;
}

__kernel void sheet_force(
    __global float4* pos, 
    __global float4* force
){    
    const int i = get_global_id(1)*get_global_size(0) + get_global_id(0);
    //printf("global_id (%li,%li) %i \n", get_global_id(0), get_global_id(1), i );
    float3 p = pos[i].xyz;
    float3 f = (float3) (0.0f,0.0f,0.0f);
    float3 pj,c;    
    
    if(     get_global_id(1)>0                    ){ 
        if( get_global_id(0)>0 ){ pj =pos[i-get_global_size(0)-1 ].xyz; f += stickForce( pj, p, param_kX, param_l0X ); }
                                  pj =pos[i-get_global_size(0)   ].xyz; f += stickForce( pj, p, param_k,  param_l0 ); c=pj; 
        if( get_global_id(0)<get_global_size(0)-1 ){ pj =pos[i-get_global_size(0)+1 ].xyz; f += stickForce( pj, p, param_kX, param_l0X ); }
    }
    if(     get_global_id(1)<get_global_size(1)-1){ 
        if( get_global_id(0)>0                  ){ pj =pos[i+get_global_size(0)-1 ].xyz; f += stickForce( pj, p, param_kX, param_l0X ); }
                                                     pj =pos[i+get_global_size(0)].xyz;    f += stickForce( pj, p, param_k,  param_l0  ); 
        if( get_global_id(1)>0                    ){ f+=harmonicForce( p, (c+pj)*0.5f, K_bend ); }   
        if( get_global_id(0)<get_global_size(0)-1 ){ pj =pos[i+get_global_size(0)+1 ].xyz; f += stickForce( pj, p, param_kX, param_l0X ); }
    }
    if(     get_global_id(0)>0                    ){ pj =pos[i-1      ].xyz; f += stickForce( pj, p, param_k, param_l0 ); c=pj; }
    if(     get_global_id(0)<get_global_size(0)-1 ){ pj =pos[i+1      ].xyz; f += stickForce( pj, p, param_k, param_l0 ); 
        if( get_global_id(0)>0                    ){ f+=harmonicForce( p, (c+pj)*0.5f, K_bend ); }   
    }
    f.y += -0.05f;
    
    force[i] = (float4)(f,0.0f);
    //force[i] = 1545.0f;
}

__kernel void move_leapfrog(
    __global float4* pos,
    __global float4* vel, 
    __global float4* force,
    float dt, float damp
){    
    const int i = get_global_id(0);
    float3 f = force[i].xyz;
    float3 v = vel  [i].xyz;
    float3 p = pos  [i].xyz;
    v        = v*damp + f*dt;
    p       += v*dt;
    vel[i]   = (float4)(v,0.0f);
    pos[i]   = (float4)(p,0.0f);
}


__kernel void cloth_dynamics(
    __global float4* pos,
    __global float4* vel, 
    float dt, float damp
){    
    const int i = get_global_id(1)*get_global_size(0) + get_global_id(0);
    float3 p = pos[i].xyz;
    float3 f = (float3) (0.0f,0.0f,0.0f);
    
    if( get_global_id(1)>0                   ){ f += stickForce( pos[i-get_global_size(0)].xyz, p, param_k, param_l0 ); }
    if( get_global_id(1)<get_global_size(1)-1){ f += stickForce( pos[i+get_global_size(0)].xyz, p, param_k, param_l0 ); }
    if( get_global_id(0)>0                   ){ f += stickForce( pos[i-1                 ].xyz, p, param_k, param_l0 ); }
    if( get_global_id(0)<get_global_size(0)-1){ f += stickForce( pos[i+1                 ].xyz, p, param_k, param_l0 ); }
    f.y += -0.05f;
    
    if( get_global_id(1)<get_global_size(1)-1 ){
        float3 v = vel[i].xyz;        
        v        = v*damp + f*dt;
        p       += v*dt;
        vel[i]   = (float4)(v,0.0f);
        pos[i]   = (float4)(p,0.0f);
    }
}

__kernel void harmonic_constr(
    __global int   * iconstrains,
    __global float4* constrains,
    __global float4* pos,
    __global float4* force 
){    
    const int ii   = get_global_id(0);
    float4 constr  = constrains[ii];
    const int  i   = iconstrains[ii];
    force[i]      += (float4)( harmonicForce( pos[i].xyz, constr.xyz, constr.w )  , 0.0f  );
}









