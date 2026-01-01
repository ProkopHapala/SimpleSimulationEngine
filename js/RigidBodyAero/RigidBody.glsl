#version 300 es
precision highp float;

// ================= INPUTS (Previous Frame) =================
uniform sampler2D u_tex_pos;
uniform sampler2D u_tex_vel;
uniform sampler2D u_tex_quat;
uniform sampler2D u_tex_ang_vel;
uniform sampler2D u_tex_force;  // Stored Force from prev step
uniform sampler2D u_tex_torque; // Stored Torque from prev step

// ================= SIMULATION PARAMETERS =================
uniform float u_dt;
uniform vec4  u_gravity; // e.g., (0, -9.81, 0, 0)
uniform bool  u_wind_tunnel_mode;

// ================= BODY DEFINITION (Wing Segments) =================
#define MAX_POINTS 16 
uniform int  u_point_count;
uniform vec3 u_points_pos[MAX_POINTS]; // Local positions
uniform vec3 u_points_t[MAX_POINTS];   // Chord/Tangent (World basis)
uniform vec3 u_points_n[MAX_POINTS];   // Normal/Up (World basis)

struct AeroParams {
    float area;
    float CD0;
    float dCD;
    float dCDS;
    float dCL;
    float dCLS;
    float sStall;
    float wStall;
};
uniform AeroParams u_aero_params[MAX_POINTS];

uniform vec3 u_inertia_inv;        // Diagonal of Body-Space Inverse Inertia

// ================= OUTPUTS (Next Frame) =================
// We write to multiple textures simultaneously (MRT)
layout(location = 0) out vec4 out_pos;
layout(location = 1) out vec4 out_vel;
layout(location = 2) out vec4 out_quat;
layout(location = 3) out vec4 out_ang_vel;
layout(location = 4) out vec4 out_force;
layout(location = 5) out vec4 out_torque;

// ================= MATH HELPERS =================

// Quaternion to Rotation Matrix
mat3 quat_to_mat3(vec4 q) {
    float x2 = q.x + q.x; float y2 = q.y + q.y; float z2 = q.z + q.z;
    float xx = q.x * x2;  float xy = q.x * y2;  float xz = q.x * z2;
    float yy = q.y * y2;  float yz = q.y * z2;  float zz = q.z * z2;
    float wx = q.w * x2;  float wy = q.w * y2;  float wz = q.w * z2;

    // Standard Body-to-World (Column-Major for GLSL)
    return mat3(
        1.0 - (yy + zz), xy + wz,          xz - wy,
        xy - wz,         1.0 - (xx + zz),  yz + wx,
        xz + wy,         yz - wx,          1.0 - (xx + yy)
    );
}

// Quaternion Multiplication
vec4 quat_mul(vec4 q1, vec4 q2) {
    return vec4(
        q1.w * q2.xyz + q2.w * q1.xyz + cross(q1.xyz, q2.xyz),
        q1.w * q2.w - dot(q1.xyz, q2.xyz)
    );
}

void main() {
    // 1. Get current particle index based on pixel coordinate
    ivec2 uv = ivec2(gl_FragCoord.xy);

    // 2. Fetch Previous State
    vec4  p_data = texelFetch(u_tex_pos,     uv, 0);
    vec4  v_data = texelFetch(u_tex_vel,     uv, 0);
    vec4  q_data = texelFetch(u_tex_quat,    uv, 0);
    vec4  w_data = texelFetch(u_tex_ang_vel, uv, 0);
    vec3  F_old  = texelFetch(u_tex_force,   uv, 0).xyz;
    vec3  T_old  = texelFetch(u_tex_torque,  uv, 0).xyz;

    vec3  pos  = p_data.xyz;
    float mass = p_data.w;
    vec3  vel  = v_data.xyz;
    vec4  quat = q_data;
    vec3  ang  = w_data.xyz;


    // Safety check for empty pixels
    if (mass <= 0.0) { discard; }

    // =========================================================
    // STEP 1: FIRST KICK (Velocity Verlet)
    // Advance velocities to t + 0.5*dt
    // =========================================================

    // Linear
    // Note: F_old includes gravity from previous step accumulation
    vec3 v_half = vel + (F_old / mass) * (0.5 * u_dt);

    // Angular
    // We need to convert Torque (World) -> Alpha (World)
    // alpha = R * I_inv_body * R^T * T
    mat3 R = quat_to_mat3(quat);
    mat3 R_t = transpose(R);
    
    // Transform Torque to Body Space
    vec3 T_body_old = R_t * T_old;
    // Apply Inertia
    vec3 alpha_body_old = u_inertia_inv * T_body_old;
    vec3 alpha_world_old = R * alpha_body_old;

    vec3 w_half = ang + alpha_world_old * (0.5 * u_dt);

    // =========================================================
    // STEP 2: DRIFT
    // Advance Position to t + dt using v_half
    // =========================================================
    
    if (!u_wind_tunnel_mode) {
        pos += v_half * u_dt;

        // Update Quaternion
        // dq/dt = 0.5 * w * q
        vec4 w_quat = vec4(w_half, 0.0);
        vec4 dq = quat_mul(w_quat, quat);
        quat += dq * (0.5 * u_dt);
        quat = normalize(quat); // Important to combat drift!
    } else {
        pos = vec3(0.0);
    }

    // Recompute Rotation Matrix for the NEW orientation
    // We need this to transform points for Force evaluation
    R = quat_to_mat3(quat); 

    // =========================================================
    // STEP 3: EVAL FORCE (at t+1, using v_half)
    // =========================================================
    
    vec3 F_new = vec3(0.0);
    vec3 T_new = vec3(0.0);

    // A. Gravity
    F_new += u_gravity.xyz * mass;

    // B. Aerodynamics (Loop over wing segments)
    for(int i = 0; i < u_point_count; i++) {
        vec3 r_local = u_points_pos[i];
        vec3 r_world = R * r_local;
        
        // Wind Tunnel: constant incident wind
        vec3 v_wind = u_wind_tunnel_mode ? vec3(0.0, 0.0, -5.0) : vec3(0.0);
        vec3 v_point = v_half + cross(w_half, r_world);
        vec3 v_air = v_wind - v_point;
        
        float v2 = dot(v_air, v_air);
        if(v2 > 0.0001) {
            vec3 u_air = normalize(v_air);
            
            // Local frame in world space
            vec3 t_world = R * u_points_t[i];
            vec3 n_world = R * u_points_n[i];
            
            float ca = dot(u_air, n_world);
            float sa = dot(u_air, t_world);
            
            float abs_sa = abs(sa);
            AeroParams p = u_aero_params[i];
            
            // Stall blending
            float wS = 0.0;
            if (abs_sa > p.sStall) {
                float x = (abs_sa - p.sStall) / p.wStall;
                wS = (x > 1.0) ? 1.0 : (3.0*x*x - 2.0*x*x*x);
            }
            float mS = 1.0 - wS;
            
            // Drag
            float CD = p.CD0 + (mS * p.dCD * abs_sa + wS * p.dCDS) * abs_sa;
            
            // Lift
            float CL_ca = ca;
            float CL_sa = sa;
            if (CL_ca < 0.0) { CL_ca = -CL_ca; CL_sa = -CL_sa; }
            float CL = (mS * p.dCL + wS * p.dCLS * CL_ca) * CL_sa;
            
            float prefactor = v2 * p.area;
            
            // Lift direction
            vec3 liftDir = normalize(t_world - sa * u_air);
            
            vec3 f_total_pt = (liftDir * CL + u_air * CD) * prefactor;

            F_new += f_total_pt;
            T_new += cross(r_world, f_total_pt);
        }
    }

    // =========================================================
    // STEP 4: SECOND KICK
    // Advance velocities from t+0.5 to t+1
    // =========================================================
    
    vec3 v_cm = v_half;
    vec3 w_cm = w_half;

    if (!u_wind_tunnel_mode) {
        // a. Linear kick
        v_cm = v_half + (F_new / mass) * (0.5 * u_dt);

        // b. Angular kick (body frame)
        mat3 Rt = transpose(R);
        vec3 T_body = Rt * T_new;
        vec3 alpha_body_new = u_inertia_inv * T_body;
        
        vec3 alpha_world_new = R * alpha_body_new;
        w_cm = w_half + alpha_world_new * (0.5 * u_dt);
    } else {
        v_cm = vec3(0.0);
        w_cm = vec3(0.0);
    }

    // Final pack
    out_pos     = vec4(pos, mass);
    out_vel     = vec4(v_cm, 0.0);
    out_quat    = quat;
    out_ang_vel = vec4(w_cm, 0.0);
    out_force   = vec4(F_new, 0.0);
    out_torque  = vec4(T_new, 0.0);
}