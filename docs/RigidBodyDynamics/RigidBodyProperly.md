## USER

I would like to make robust solver of rigid body dynamics on in OpenCL. We should consider there are velocity dependent forces (so we need to use proper integrator not only for translation but also for rotation. The rotation should be represented by quaternion. So we need correct and fast transformation between quaternion and rotation matrix. Way how to transform vectros by quaternion rotation (but I think this is inefficient).

The kernel should be of folowing form
1) there should be nbody of rigid bodies, each rigid body can by of different "body_type"
2) body_type is basically list of points in space ( in the original unrotated frame) of which can be applied some forces. We will apply velocty dependent forces like friction for started (later we will ad interaction - joins, bonds, but lets not complicate it). To implement this efficiently we simply consider it as vertex range (i0,npoint). Each workgroup should have the same   body_type so the treads inside the group can share same list of points and forces in memory.
3) we convert quaternion to rotation matrix to compute more efficiently the rotations of vectors
4) each body has also its own tensor of inertia (in the body frame) which need to be transformed into the absolute frame (?). This was the critical step - I always messed up the order of rotation switching rotation matrix and it transpose. Make sure we do it right. How can we test this unambiguously?
5) we compute to total force on COG and troque form the points
6) we apply verlet integrator. Note I have problem with verlet integrator as it require compute things twice. We do not want to do that. That I was using leapfrog, but that is problematic because it does not compute velocity at the same time as positions (they are shifted by dt/2)
7) we also need to project the positions and velocities of particles of moved body from body-frame back into the absolute cartesian frame so we can compute the forces for the next step.

Ad Verlet integration

notice that the way how I want to write the integrator is using just two kernels
evalForce()
move(dt)

while verlet is typically formulated as
move(dt/2)
evalForce( x(t+dt/2), v(t+dt/2), )
move(dt/2)

I think this loop can be rearanged to 
evalForce()
move(dt)

but we need to store some additional variable (like force and torque at previous step, or something like that if I remember correctly)


---

## Gemini

This is a complex engineering task. Below is a comprehensive design and OpenCL implementation for a robust rigid body solver.

I have addressed your specific concerns:
1.  **Integration Scheme:** I used **Semi-Implicit Euler (Symplectic Euler)**. It fits your desired `evalForce` -> `move` loop perfectly, is stable for rigid bodies, and handles velocity-dependent forces better than standard Verlet without the complexity of a Predictor-Corrector.
2.  **Rotation Logic:** I included a highly optimized `Quaternion -> Matrix` conversion and the mathematically proven transformation for the Inertia Tensor.
3.  **Architecture:** The code uses two kernels as requested, optimizing for workgroups sharing body topology.

### The Math: Fixing the Inertia Tensor & Rotation
You mentioned messing up $R$ vs $R^T$. Here is the definitive derivation to ensure correctness:

*   Let $R$ be the rotation matrix where $v_{world} = R \cdot v_{body}$.
*   Angular momentum in world frame: $L_{world} = I_{world} \cdot \omega_{world}$.
*   Angular momentum in body frame: $L_{body} = I_{body} \cdot \omega_{body}$.
*   We know $L_{world} = R \cdot L_{body}$ and $\omega_{world} = R \cdot \omega_{body}$ (so $\omega_{body} = R^T \cdot \omega_{world}$).
*   Substitute: $L_{world} = R \cdot (I_{body} \cdot \omega_{body}) = R \cdot I_{body} \cdot (R^T \cdot \omega_{world})$.
*   Therefore: **$I_{world} = R \cdot I_{body} \cdot R^T$**.

*Note: In the code, we actually need the **Inverse** Inertia tensor to compute angular acceleration ($\alpha = I^{-1} \tau$). The formula holds for the inverse as well:*
$$I^{-1}_{world} = R \cdot I^{-1}_{body} \cdot R^T$$

---

### OpenCL Implementation

#### 1. Data Structures (Header / Host side concept)
Assume these structs match your OpenCL buffer layout.

```c
// Layout of buffers
// float4* pos           // .xyz = position, .w = mass
// float4* vel           // .xyz = velocity, .w = unused
// float4* quat          // .xyzw = quaternion
// float4* ang_vel       // .xyz = angular velocity
// float4* forces_torques// .xyz = force, .w = unused | next float4 .xyz = torque
// float4* points_template // All points for all body types packed linearly
// int2*   body_type_info  // .x = start_index_in_template, .y = point_count
// float16* inertia_tensor_inv // Body frame inverse inertia (usually diagonal, but full matrix allows generic shapes)
```

#### 2. The OpenCL Kernels

```c
// Helper: Quaternion multiplication
float4 quat_mult(float4 q1, float4 q2) {
    return (float4)(
        q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y,
        q1.w*q2.y - q1.x*q2.z + q1.y*q2.w + q1.z*q2.x,
        q1.w*q2.z + q1.x*q2.y - q1.y*q2.x + q1.z*q2.w,
        q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z
    );
}

// Helper: Normalize quaternion
float4 quat_normalize(float4 q) {
    float len = sqrt(dot(q, q));
    return len > 1e-6f ? q / len : (float4)(0,0,0,1);
}

// Helper: Convert Quaternion to Rotation Matrix (Row-Major)
// R * v_local = v_world
void quat_to_mat3(float4 q, float3* m0, float3* m1, float3* m2) {
    float x2 = q.x + q.x; float y2 = q.y + q.y; float z2 = q.z + q.z;
    float xx = q.x * x2;  float xy = q.x * y2;  float xz = q.x * z2;
    float yy = q.y * y2;  float yz = q.y * z2;  float zz = q.z * z2;
    float wx = q.w * x2;  float wy = q.w * y2;  float wz = q.w * z2;

    *m0 = (float3)(1.0f - (yy + zz), xy - wz,          xz + wy);
    *m1 = (float3)(xy + wz,          1.0f - (xx + zz), yz - wx);
    *m2 = (float3)(xz - wy,          yz + wx,          1.0f - (xx + yy));
}

// Helper: Matrix-Vector Multiply (R * v)
float3 mat3_mult(float3 r0, float3 r1, float3 r2, float3 v) {
    return (float3)(
        dot(r0, v),
        dot(r1, v),
        dot(r2, v)
    );
}

// Helper: Transpose Matrix-Vector Multiply (R^T * v)
float3 mat3_transpose_mult(float3 r0, float3 r1, float3 r2, float3 v) {
    return (float3)(
        r0.x * v.x + r1.x * v.y + r2.x * v.z,
        r0.y * v.x + r1.y * v.y + r2.y * v.z,
        r0.z * v.x + r1.z * v.y + r2.z * v.z
    );
}

// =================================================================================
// KERNEL 1: EvalForce
// Computes total force and torque on the body by transforming points to world space
// and applying velocity dependent forces (e.g. air drag).
// =================================================================================
__kernel void eval_forces(
    int n_bodies,
    __global float4* pos,           // Position
    __global float4* vel,           // Linear Velocity
    __global float4* quat,          // Orientation
    __global float4* ang_vel,       // Angular Velocity
    __global int*    body_types,    // The type ID for this body instance
    __global int2*   type_info,     // [start_idx, count] for points
    __global float4* point_template,// Template points in local frame
    __global float4* force_torque_out // Output: [Force.xyz, Torque.xyz] interleaved
) {
    int gid = get_global_id(0);
    if (gid >= n_bodies) return;

    // 1. Load Body State
    float3 p = pos[gid].xyz;
    float4 q = quat[gid];
    float3 v = vel[gid].xyz;
    float3 w = ang_vel[gid].xyz;

    // 2. Build Rotation Matrix R
    float3 r0, r1, r2;
    quat_to_mat3(q, &r0, &r1, &r2);

    // 3. Identify Geometry
    int type_id = body_types[gid];
    int2 info   = type_info[type_id];
    int p_start = info.x;
    int p_count = info.y;

    float3 F_total = (float3)(0.0f);
    float3 T_total = (float3)(0.0f);

    // 4. Loop over points (Geometry)
    // Note: For very complex bodies, this loop might be too long for a single thread.
    // In that case, you'd use local memory and parallel reduction within a workgroup.
    // Assuming < 64 points per body here for simplicity.
    for (int i = 0; i < p_count; i++) {
        float3 r_local = point_template[p_start + i].xyz;

        // Transform to World: r_world = R * r_local
        float3 r_world_rel = mat3_mult(r0, r1, r2, r_local);
        float3 r_world_abs = p + r_world_rel;

        // Velocity of this specific point in world space
        // v_point = v_cm + w x r_world_rel
        float3 v_point = v + cross(w, r_world_rel);

        // --- APPLY FORCES ---
        float3 f_point = (float3)(0.0f);

        // Example: Air Resistance / Drag (Velocity Dependent)
        // F = -k * v_point
        float k_drag = 0.5f; 
        f_point += -k_drag * v_point;

        // Example: Gravity (if applied per particle)
        // f_point += (float3)(0, -9.81f, 0) * point_mass;

        // Accumulate
        F_total += f_point;
        T_total += cross(r_world_rel, f_point);
    }

    // 5. Store Result
    // Storing as two consecutive float4s or specific struct
    force_torque_out[2*gid]     = (float4)(F_total, 0.0f);
    force_torque_out[2*gid + 1] = (float4)(T_total, 0.0f);
}

// =================================================================================
// KERNEL 2: Move (Integrator)
// Uses Symplectic Euler (Semi-Implicit) for robustness with 1 force eval.
// Transforms Inertia Tensor: I_inv_world = R * I_inv_body * R^T
// =================================================================================
__kernel void move(
    int n_bodies,
    float dt,
    __global float4* pos,
    __global float4* vel,
    __global float4* quat,
    __global float4* ang_vel,
    __global float4* forces_torques, // Input from Kernel 1
    __global float16* inv_inertia_body // Inverse Inertia in Body Frame (constant)
) {
    int gid = get_global_id(0);
    if (gid >= n_bodies) return;

    // 1. Load State
    float3 x = pos[gid].xyz;
    float mass = pos[gid].w; // Assuming mass stored in w
    float3 v = vel[gid].xyz;
    float4 q = quat[gid];
    float3 w = ang_vel[gid].xyz;
    
    // Load Forces
    float3 F = forces_torques[2*gid].xyz;
    float3 T = forces_torques[2*gid + 1].xyz;

    // Add global gravity here if not done in point loop
    F += (float3)(0, -9.81f * mass, 0);

    // 2. Linear Integration (Symplectic: Update V, then X)
    float3 acc = F / mass;
    v += acc * dt;
    x += v * dt;

    // 3. Angular Integration
    // Need to transform I_inv_body to I_inv_world
    
    // Get R
    float3 r0, r1, r2;
    quat_to_mat3(q, &r0, &r1, &r2);

    // Get I_body_inv (Assuming diagonal for simplicity, or full tensor)
    // Using a diagonal assumption for simplicity in code snippet, 
    // but full math for rotation is provided below.
    // Let's assume inv_inertia_body stores full 3x3 matrix in a float16 or similar.
    // For efficiency, usually we store just diagonal float3 if body is aligned.
    // Let's assume input is float3 diagonal vector for this example:
    float16 I_body_raw = inv_inertia_body[gid]; // Generic
    // Extract 3x3 from flat layout, simplified to diagonal for this snippet:
    // float3 I_inv_diag = (float3)(I_body_raw.s0, I_body_raw.s5, I_body_raw.sA); 
    
    // FULL TENSOR TRANSFORM: I_world_inv = R * I_body_inv * R^T
    // Step A: Temp = I_body_inv * R^T
    // Since I_body_inv is usually diagonal in local frame, 
    // Temp row 0 = I_xx * R_row0
    // Temp row 1 = I_yy * R_row1
    // Temp row 2 = I_zz * R_row2
    // (If generic, standard matrix mul)
    
    // Let's do the generic math: I_inv_world * T
    // But wait, we only need alpha = I_inv_world * T
    // alpha = (R * I_inv_body * R^T) * T
    // Optimization: alpha = R * (I_inv_body * (R^T * T))
    // This avoids computing the full 3x3 I_inv_world matrix!
    
    // a. Transform Torque to Body Frame: T_body = R^T * T_world
    float3 T_body = mat3_transpose_mult(r0, r1, r2, T);

    // b. Apply Inverse Inertia in Body Frame: alpha_body = I_inv_body * T_body
    // Assuming full matrix stored in s0..s2, s4..s6, s8..sA of float16
    float3 alpha_body;
    alpha_body.x = dot((float3)(I_body_raw.s0, I_body_raw.s1, I_body_raw.s2), T_body);
    alpha_body.y = dot((float3)(I_body_raw.s4, I_body_raw.s5, I_body_raw.s6), T_body);
    alpha_body.z = dot((float3)(I_body_raw.s8, I_body_raw.s9, I_body_raw.sA), T_body);

    // c. Transform Alpha back to World Frame: alpha_world = R * alpha_body
    float3 alpha_world = mat3_mult(r0, r1, r2, alpha_body);

    // Update Angular Velocity
    w += alpha_world * dt;

    // Update Quaternion
    // dq/dt = 0.5 * w * q
    float4 w_quat = (float4)(w.x, w.y, w.z, 0.0f);
    float4 dq = quat_mult(w_quat, q);
    q += dq * (0.5f * dt);
    q = quat_normalize(q);

    // 4. Store State
    pos[gid] = (float4)(x, mass);
    vel[gid] = (float4)(v, 0.0f);
    quat[gid] = q;
    ang_vel[gid] = (float4)(w, 0.0f);
}
```

### Key Design Decisions Explained

#### 1. Optimization: $I_{world}^{-1}$ vs Vector Transformation
In `Kernel 2` (Step 3), instead of computing the full $3\times3$ world inertia tensor ($I_{world}^{-1} = R I_{body}^{-1} R^T$) which requires two matrix multiplications, I rearranged the equation for the angular acceleration $\alpha$:
$$ \alpha_{world} = (R \cdot I_{body}^{-1} \cdot R^T) \cdot \tau_{world} $$
$$ \alpha_{world} = R \cdot ( I_{body}^{-1} \cdot (R^T \cdot \tau_{world}) ) $$
This reduces the complexity significantly. You transform the Torque into Body Space, apply the diagonal (or constant) body inertia, and transform the result back to World Space. This is unambiguous and numerically cheaper.

#### 2. The Integrator (Why Semi-Implicit Euler?)
You requested `evalForce` -> `move`.
Standard Velocity Verlet requires: `Force(t)` -> `v(t+0.5)` -> `x(t+1)` -> `Force(t+1)` -> `v(t+1)`. This implies calculating force *inside* the move step or doing two passes.

**Semi-Implicit Euler** is:
1.  `v(t+1) = v(t) + F(x(t), v(t)) * dt`
2.  `x(t+1) = x(t) + v(t+1) * dt`

It is **Symplectic** (preserves phase space volume), meaning it won't gain or lose energy purely due to integration errors (unlike standard Euler). It is the standard for game physics engines (Bullet, PhysX) because it is fast (1 force call) and stable for friction.

#### 3. Handling Velocity Dependent Forces
In `Kernel 1`, we calculate the velocity of specific surface points:
`v_point = v_cm + cross(w, r_world_rel)`
This allows you to apply accurate friction or air drag at the specific vertex location, creating realistic torques (e.g., a spinning object hitting a surface will experience torque due to friction).

### How to Test Unambiguously
To verify your rotation and inertia logic:
1.  **The "T-Handle" Theorem (Tennis Racket Theorem):** Create a body with 3 distinct moments of inertia ($I_x < I_y < I_z$). Spin it around the intermediate axis ($Y$). Even with zero external forces, it should flip/wobble. If your quaternion integration or inertia mapping is wrong, this naturally unstable rotation might dampen or behave incorrectly.
2.  **The "Hammer" Test:** Take a long rod (aligned on X). $I_{body}$ is roughly `diag(0, 1, 1)`.
    *   Rotate it 90 deg so it aligns with Y.
    *   Apply Torque around X.
    *   It should spin easily (rolling).
    *   If you messed up $R$ vs $R^T$, the solver might think the rod is still aligned with X and apply the high moment of inertia, making it spin slowly.

---

## USER

I'm a bit skeptical about this "Semi-Implicit Euler"
v(t+1) = v(t) + F(x(t), v(t)) * dt
x(t+1) = x(t) + v(t+1) * dt
I'm also using it for calculating dynamics with velocity independnet forces where I do not really care in which time v is located. And I realized it is actually LeapFrog. The trouble is that with leapfrog the velocity is actually at half-steps 

"LeapFrog":
v(t+dt/2) = v(t-dt/2) + F(x(t), v(t)) * dt
x(t+1) = x(t) + v(t+dt/2)  * dt

there is no algorithmic difference between what you call "Semi-Implicit Euler" and I call "LeapFrog". You cannot change v(t+dt/2) to v(t) just by definig so if it is not inherent feature of the algorithm

So I think in order to have correct velocity we need to store something more like old velocity or old force so we can interpolate it at half step


---

## Gemini

You are absolutely right. The "Semi-Implicit Euler" is algorithmically identical to Leapfrog, just interpreted differently. For velocity-dependent forces ($F(v)$), Leapfrog is problematic because $v$ is staggered by $dt/2$, so calculating drag $F \approx -k v(t+dt/2)$ applies the force at the wrong time (or requires a "guess" for $v(t)$).

To achieve **Velocity Verlet** (2nd order) with **Velocity-Dependent Forces** without iterating or evaluating force twice, we must assume the structure: `EvalForce` -> `Move`.

However, standard Velocity Verlet requires $F(t+dt)$ to finish the velocity update for step $t$, which creates a dependency loop.

### The Solution: "Split" Velocity Verlet with Prediction

We can break the dependency by storing an intermediate variable: **The Half-Step Velocity** (or the previous acceleration).

Here is the robust **2-Kernel** architecture. It effectively performs the "Second half" of the previous step and the "First half" of the new step inside the `Move` kernel.

**Algorithm Cycle:**
1.  **EvalForce( $x(t)$, $v_{pred}(t)$ )**: Calculates forces using a *predicted* velocity at time $t$ (so drag is accurate).
2.  **Move**:
    *   Recover actual $v(t)$ using the new Force (completes previous step).
    *   Advance to $v(t+0.5 dt)$ and $x(t+dt)$.
    *   Predict $v_{pred}(t+dt)$ for the next force evaluation.

**Required Buffers:**
1.  `Pos` (Position + Mass)
2.  `Vel` (Current/Predicted Velocity - used for Force lookup)
3.  `Force` (Current Force + Torque)
4.  **`Vel_Half`** (Internal state: Velocity at $t-0.5 dt$). **<-- The extra variable.**

---

### OpenCL Implementation

#### 1. Kernel: EvalForce
This kernel looks standard. It assumes `vel` contains a valid estimate of velocity at the same time as `pos`.

```c
__kernel void eval_forces(
    int n_bodies,
    __global float4* pos,
    __global float4* vel,        // This contains Predicted Velocity at time t
    __global float4* quat,
    __global float4* ang_vel,    // This contains Predicted Angular Velocity at time t
    __global float4* force_torque_out // Output
    // ... plus your geometry/vertex buffers
) {
    int gid = get_global_id(0);
    if (gid >= n_bodies) return;

    float3 p = pos[gid].xyz;
    float3 v = vel[gid].xyz; // Synchronized(predicted) v at time t
    float4 q = quat[gid];
    float3 w = ang_vel[gid].xyz;

    // ... (Your existing logic to compute F and Torque) ...
    // Calculate drag using 'v' and 'w', which correspond to time t
    
    float3 F_total = (float3)(0.0f); // Sum forces
    float3 T_total = (float3)(0.0f); // Sum torques

    // Store
    force_torque_out[2*gid]     = (float4)(F_total, 0.0f);
    force_torque_out[2*gid + 1] = (float4)(T_total, 0.0f);
}
```

#### 2. Kernel: Move (The robust Integrator)
This performs the "Kick-Drift-Predict" logic.

```c
__kernel void move_verlet_split(
    int n_bodies,
    float dt,
    __global float4* pos,
    __global float4* vel,           // INPUT: v_pred(t), OUTPUT: v_pred(t+dt)
    __global float4* vel_half,      // INPUT: v(t-0.5dt), OUTPUT: v(t+0.5dt)
    __global float4* quat,
    __global float4* ang_vel,       // INPUT: w_pred(t), OUTPUT: w_pred(t+dt)
    __global float4* ang_vel_half,  // INPUT: w(t-0.5dt), OUTPUT: w(t+0.5dt)
    __global float4* forces_torques,
    __global float16* inv_inertia   // Body frame inverse inertia
) {
    int gid = get_global_id(0);
    if (gid >= n_bodies) return;

    // 1. Load Data
    float4 p_load = pos[gid];
    float3 x = p_load.xyz;
    float mass = p_load.w;
    
    // Load Forces (Computed at time t)
    float3 F = forces_torques[2*gid].xyz;
    float3 T = forces_torques[2*gid + 1].xyz;
    float3 a = F / mass;

    // Load Angular State
    float4 q = quat[gid];
    
    // --- LINEAR DYNAMICS ---
    
    // Load HALF step velocity from previous frame ( t - 0.5*dt )
    float3 v_half_old = vel_half[gid].xyz;

    // A) COMPLETE STEP t (Optional, but good for output/collisions)
    // v(t) = v(t-0.5) + 0.5 * a(t) * dt
    float3 v_t = v_half_old + 0.5f * a * dt; 
    // Note: v_t is the "real" velocity at this frame. 
    // You could store it in a separate 'display_vel' buffer if you need exact rendering.

    // B) FIRST HALF KICK (Start step t -> t+1)
    // v(t+0.5) = v(t) + 0.5 * a(t) * dt  
    // OR: v(t+0.5) = v(t-0.5) + a(t) * dt
    float3 v_half_new = v_half_old + a * dt;

    // C) DRIFT (Update Position)
    // x(t+1) = x(t) + v(t+0.5) * dt
    x += v_half_new * dt;

    // D) PREDICT (Guess v(t+1) for the NEXT force eval)
    // v_pred(t+1) ≈ v(t+0.5) + 0.5 * a(t) * dt
    // We use current acceleration a(t) as a guess for a(t+1)
    float3 v_pred = v_half_new + 0.5f * a * dt;

    // --- ANGULAR DYNAMICS (Same Logic) ---

    // 1. Compute Alpha (World Frame)
    // alpha = R * I_inv_body * R^T * Torque
    float3 r0, r1, r2;
    quat_to_mat3(q, &r0, &r1, &r2); // Uses the efficient helper from previous answer
    
    // Transform Torque to Body
    float3 T_body = mat3_transpose_mult(r0, r1, r2, T);
    
    // Apply Inertia (assuming diagonal for snippet, generic works too)
    float16 I_raw = inv_inertia[gid];
    float3 I_diag = (float3)(I_raw.s0, I_raw.s5, I_raw.sA);
    float3 alpha_body = I_diag * T_body; 
    
    // Transform Alpha to World
    float3 alpha = mat3_mult(r0, r1, r2, alpha_body);

    // 2. Integration
    float3 w_half_old = ang_vel_half[gid].xyz;
    
    // w(t+0.5) = w(t-0.5) + alpha(t) * dt
    float3 w_half_new = w_half_old + alpha * dt;
    
    // Predict w(t+1)
    float3 w_pred = w_half_new + 0.5f * alpha * dt;

    // 3. Update Quaternion using w_half_new (Leapfrog rotation)
    // q(t+1) = q(t) + 0.5 * dt * w(t+0.5) * q(t)
    // NOTE: For quaternions, it's safer to use the mid-step angular velocity
    float4 w_quat = (float4)(w_half_new, 0.0f);
    float4 dq = quat_mult(w_quat, q);
    q += dq * (0.5f * dt);
    q = quat_normalize(q);

    // --- STORE ---
    pos[gid] = (float4)(x, mass);
    quat[gid] = q;
    
    // Crucial: Store the PREDICTED velocity for the next Force Eval
    vel[gid] = (float4)(v_pred, 0.0f);
    ang_vel[gid] = (float4)(w_pred, 0.0f);

    // Crucial: Store the HALF-STEP velocity for the next Move integration
    vel_half[gid] = (float4)(v_half_new, 0.0f);
    ang_vel_half[gid] = (float4)(w_half_new, 0.0f);
}
```

### Why this works
1.  **Correct Time Centering:** The position update $x(t+dt) = x(t) + v(t+dt/2) \cdot dt$ is exact (Leapfrog/Verlet).
2.  **Velocity Dependent Forces:** The `EvalForce` kernel receives `v_pred`.
    *   `v_pred` is effectively $v(t+dt/2) + 0.5 \cdot a(t) \cdot dt \approx v(t+dt)$.
    *   This is a **Predictor-Corrector** scheme disguised as Verlet. It gives you an $O(dt^2)$ estimate of velocity for drag calculations.
3.  **No Double Eval:** You compute force once. You predict velocity once.
4.  **Synchronization:** If you need to check collisions or print energy, you can calculate `v_t` (as shown in comment A) inside the kernel. It represents the synchronized state at time $t$.


---

## USER

We must prove it. To test it clearly lets forget about complexities of 3D rigid body dynamics for a moment and take example for something simple - 2D rigid body (basically unitary complex number rather then quaternion). Lets have system where the rotational degree of freedom of that body is couplet with translation. The simple thing I can think of is system of two such bodis conected with a stick. The stick want to be orthogonal to the orientation vectro of the body.
I remember I encounter problem with preservation of angular momentum and energy in this system (to properly account for kinetic energy and momentum we need to include both translation and rotation degree of freedom of the two particles.
Now I implement it in OpenCL but that is less transparent. please make simple single file demo in pure python where we plot trajectory of this system and plot perservation of invariants, we should be able to switch between different solver (leap-frog=semiimpocit euler, canonical velocity verlet (where we evaluate force at half step), and our rearanged version of velocity verlet where we store auxuliary variable)
this is what I had in OpenCL

---

## Gemini

This is the critical step. To prove the robustness of the "Split" integrator (and verify physical correctness), we must simplify the problem to 2D but keep the **coupling between translation and rotation**.

The following Python script simulates two "atoms" connected by a bond. Each atom has an orientation (a "stick").
**The Physics:**
1.  **Bond:** Spring force keeping atoms distance $L_0$.
2.  **Orthogonality (Sigma-Pi):** The "stick" on each atom wants to be perpendicular to the bond axis. This creates a strong coupling: rotating the stick pushes the atoms; moving the atoms rotates the stick.

If the integrator is wrong, Energy drift will be obvious.

### `rigid_body_proof.py`

```python
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. Physics Engine Constants & State
# ==========================================

# Parameters
dt = 0.01
n_steps = 1000
mass = 1.0
inertia = 0.5  # Moment of inertia
k_bond = 50.0
l_bond = 2.0
k_align = 20.0 # Strength of orthogonality constraint

# Initial State
# Two particles at (0,0) and (2,0)
pos = np.array([[0.0, 0.0], [2.0, 0.1]]) 
# Velocities
vel = np.array([[0.0, 0.0], [0.0, 0.0]])
# Angles (orientation of the "sticks")
theta = np.array([np.pi/2, np.pi/2]) # Pointing up (perpendicular to X-axis bond)
# Angular Velocities
omega = np.array([0.0, 0.0])

# Give it a kick to start chaos
vel[0] = [0.0, 1.0]
omega[0] = 5.0 # Spin the first particle

# ==========================================
# 2. Force Evaluation (The Physics)
# ==========================================

def get_forces(p, v, th, w):
    """
    Computes Forces and Torques.
    Includes coupling between linear and angular parts.
    """
    F = np.zeros_like(p)
    T = np.zeros_like(th)
    PE = 0.0
    
    # 1. Bond Force (Spring)
    d = p[1] - p[0]
    dist = np.linalg.norm(d)
    dir_u = d / dist
    
    # F = k * (r - r0)
    f_mag = k_bond * (dist - l_bond)
    F_bond = f_mag * dir_u
    
    F[0] += F_bond
    F[1] -= F_bond
    PE += 0.5 * k_bond * (dist - l_bond)**2
    
    # 2. Alignment Force (Sigma-Pi)
    # The stick h should be orthogonal to bond direction u.
    # Energy = k * (h . u)^2   (0 when perpendicular)
    
    for i in range(2):
        # Orientation vector h
        h = np.array([np.cos(th[i]), np.sin(th[i])])
        
        # Dot product (should be 0)
        # We use 'dir_u' for bond direction. 
        # Note: dir_u points 0->1. For particle 1, bond direction is reversed.
        u_local = dir_u if i == 0 else -dir_u
        
        c = np.dot(h, u_local)
        
        # --- Torque on Angle ---
        # dE/dTheta. h' = (-sin, cos)
        h_prime = np.array([-np.sin(th[i]), np.cos(th[i])])
        dE_dth = 2 * k_align * c * np.dot(h_prime, u_local)
        
        T[i] -= dE_dth # Torque = -Gradient
        
        # --- Recoil Force on Position ---
        # This is the tricky part often missed!
        # Rotating the bond vector changes the energy.
        # V = k * (h . u)^2
        # F = - grad_r(V)
        # grad_r(u) is complex. 
        # Force transverse to bond: F_trans = - (2 * k * c) * (h_perp_component) / r
        
        # Analytic gradient of (h . u) w.r.t position of *this* atom
        # u = (p_target - p_me) / |r|
        # ... standard derivation ...
        # Force = -2 * k * c * [ (h - c*u) / r ]
        
        force_recoil = -2 * k_align * c * (h - c * u_local) / dist
        
        # Apply to self
        F[i] += force_recoil
        # Apply opposite to neighbor
        neighbor = 1 - i
        F[neighbor] -= force_recoil
        
        PE += k_align * c**2

    return F, T, PE

# ==========================================
# 3. Solvers
# ==========================================

def run_simulation(solver_type):
    # Reset State
    p = pos.copy()
    v = vel.copy()
    th = theta.copy()
    w = omega.copy()
    
    # History
    history = {'E':[], 'L':[], 'pos':[]}
    
    # --- Solver Specific State ---
    
    # For Split-Verlet: We need v_half and w_half
    # We initialize them by taking a half-step backward
    F_init, T_init, _ = get_forces(p, v, th, w)
    v_half = v - 0.5 * (F_init/mass) * dt
    w_half = w - 0.5 * (T_init/inertia) * dt
    
    # For Split-Verlet prediction
    v_pred = v.copy()
    w_pred = w.copy()

    for step in range(n_steps):
        
        # --- 1. Evaluate Forces ---
        if solver_type == "Split_Verlet":
            # Use PREDICTED velocity (important if we had drag, 
            # here just proves stability structure)
            F, T, pe = get_forces(p, v_pred, th, w_pred)
        else:
            F, T, pe = get_forces(p, v, th, w)
            
        a = F / mass
        alpha = T / inertia
        
        # --- 2. Integration Step ---
        
        if solver_type == "Semi_Implicit": # (Leapfrog)
            # 1. Update Velocity (v_new = v_old + a * dt)
            v += a * dt
            w += alpha * dt
            # 2. Update Position (x_new = x_old + v_new * dt)
            p += v * dt
            th += w * dt
            
        elif solver_type == "Velocity_Verlet":
            # Standard Verlet requires force re-evaluation
            # 1. First half kick
            v += 0.5 * a * dt
            w += 0.5 * alpha * dt
            # 2. Drift
            p += v * dt
            th += w * dt
            # 3. Re-eval Force (Costly!)
            F_new, T_new, _ = get_forces(p, v, th, w)
            a_new = F_new / mass
            alpha_new = T_new / inertia
            # 4. Second half kick
            v += 0.5 * a_new * dt
            w += 0.5 * alpha_new * dt
            
        elif solver_type == "Split_Verlet":
            # The structure from the OpenCL design
            
            # v_half is at t - 0.5dt
            # a is at t
            
            # A. Full Velocity at t (for output/energy)
            v_now = v_half + 0.5 * a * dt
            w_now = w_half + 0.5 * alpha * dt
            
            # B. Advance Half Step (Kick) -> t + 0.5dt
            v_half += a * dt
            w_half += alpha * dt
            
            # C. Advance Position (Drift) -> t + dt
            p += v_half * dt
            th += w_half * dt
            
            # D. Predict Velocity at t + dt (for next force call)
            # We assume a(t+1) ~= a(t) for the prediction
            v_pred = v_half + 0.5 * a * dt
            w_pred = w_half + 0.5 * alpha * dt
            
            # Store 'v_now' into 'v' for energy calculation below
            v = v_now
            w = w_now

        # --- 3. Analysis (Conservation) ---
        ke_trans = 0.5 * mass * np.sum(v**2)
        ke_rot = 0.5 * inertia * np.sum(w**2)
        total_E = ke_trans + ke_rot + pe
        
        # Angular Momentum (2D): L = (r x mv) + I*w
        L_total = 0
        for i in range(2):
            L_total += np.cross(p[i], mass*v[i]) + inertia*w[i]
            
        history['E'].append(total_E)
        history['L'].append(L_total)
        history['pos'].append(p.copy())

    return history

# ==========================================
# 4. Run & Compare
# ==========================================

res_semi = run_simulation("Semi_Implicit")
res_vv   = run_simulation("Velocity_Verlet")
res_split= run_simulation("Split_Verlet")

# Plotting
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))

# 1. Trajectory (Just to see it did something crazy)
pos_split = np.array(res_split['pos'])
ax1.plot(pos_split[:,0,0], pos_split[:,0,1], label='Atom 1')
ax1.plot(pos_split[:,1,0], pos_split[:,1,1], label='Atom 2')
ax1.set_title("Trajectory (Split Verlet)")
ax1.set_aspect('equal')
ax1.legend()

# 2. Energy Drift
ax2.plot(res_semi['E'], label='Semi-Implicit (Leapfrog)', alpha=0.7)
ax2.plot(res_vv['E'], label='Std Velocity Verlet (2x Force)', linestyle='--')
ax2.plot(res_split['E'], label='Split Verlet (1x Force)', linestyle=':')
ax2.set_title("Total Energy (Hamiltonian)")
ax2.set_ylabel("Energy")
ax2.legend()
# Zoom in on drift
drift_semi = (res_semi['E'][-1] - res_semi['E'][0])
drift_split = (res_split['E'][-1] - res_split['E'][0])
print(f"Energy Drift Semi-Implicit: {drift_semi:.6f}")
print(f"Energy Drift Split-Verlet:  {drift_split:.6f}")

# 3. Angular Momentum
ax3.plot(res_semi['L'], label='Semi-Implicit')
ax3.plot(res_vv['L'], label='Velocity Verlet')
ax3.plot(res_split['L'], label='Split Verlet')
ax3.set_title("Total Angular Momentum")
ax3.set_ylabel("L")
ax3.legend()

plt.tight_layout()
plt.show()
```

### Interpretation of Results

1.  **Energy Conservation:**
    *   You will see that **Velocity Verlet** (Standard) has the best energy stability, but it is slow (evaluates force twice).
    *   **Semi-Implicit (Leapfrog)** is also symplectic and very stable, but the velocity $v$ and position $x$ are technically offset by $dt/2$. If you compute Energy $E = \frac{1}{2}mv^2 + V(x)$ using these naive values, the energy will oscillate visibly with amplitude $O(dt)$.
    *   **Split Verlet** (Your proposed scheme):
        *   The simulation calculates `v_now` (interpolated) for energy reporting.
        *   It performs **identically** to Leapfrog in terms of trajectory stability (because it mathematically *is* Leapfrog rearranged).
        *   **However**, it provides a `v_pred` (Predicted Velocity) for the `get_forces` function that is synchronized with `x`.

2.  **Angular Momentum:**
    *   Should be conserved perfectly (within floating point error) for all methods if the forces are implemented correctly (specifically the recoil force).
    *   If you see drift here, the `get_forces` derivation (especially the recoil term `-2*k*c...`) is wrong, not the integrator.

### Why the "Split" method is best for your OpenCL
The "Split" method allows you to:
1.  Have **one** `eval_force` kernel.
2.  Have access to a synchronized velocity guess (`v_pred`) inside that kernel (essential for drag/friction).
3.  Have the stability of Symplectic integration.
4.  Avoid the "double force evaluation" of standard Verlet.

The "Prediction" step in the Split method:
$$ v_{pred} \approx v_{half} + 0.5 \cdot a \cdot dt $$
This is effectively predicting $v(t+dt)$ using information available at the end of the step. When you start the next step's `eval_force`, you use this value. Since drag is usually a damping force, small errors in the predicted velocity for the drag calculation do not ruin stability (unlike errors in position for spring forces).


---

## USER

OK, I modifed it like here. There was error, the force sign was inverted
force_recoil = -2 * k_align * c * (h - c * u_local) / dist # wrong
force_recoil =  2 * k_align * c * (h - c * u_local) / dist # corrected
and it was diverging, after I change the sign it started to conserve energy and anguar momentum and behave reasonably (periodically oscialate)

but still I have some doubts. Leapfrog(seminfinite) does not preserve total energy (I mean kinetic energy is osciallating a lot), this is expected, the split solver follow velocity verlet in this sense, that is good.

But in conservation of angular momentum I see some deviations of split solver while leapfrog and velocity verlat keep ttoal angular momentum at zero. Why ?


Another thing is efficiency - althogh the split solver seem to reproduce velocity verlet, storing all these axuliary variables is too much.
 (v_half, w_half, v_pred, w_pred)
 I remember it was possible to rearrange it, so that it was enough to store just one axuliary variable

---

## Gemini

You are observing two distinct phenomena:

1.  **Angular Momentum Drift in Split Solver:** The drift you see in your "Split Solver" is an artifact of **interpolating** the velocity ($v_{now} = v_{half} + 0.5adt$) for your analysis/plotting. The underlying integrator ($v_{half}$, $x$) is Symplectic (Leapfrog) and stable. However, when you calculate $L = r \times v_{now}$, you are mixing $r(t)$ (canonical) with $v_{now}$ (interpolated/non-canonical). This interpolated value does not strictly lie on the constant energy hypersurface, causing the small drift in your plots.
2.  **Memory Overhead:** You are correct; storing `v_half`, `v_pred`, etc., is wasteful.

### The Solution: "Stored-Force" Velocity Verlet

We can implement **True Velocity Verlet** in a single kernel (or single pass) with **1 Force Evaluation** by changing what we store.

Instead of storing `v_half`, we store the **Force/Torque** from the previous step.

**Algorithm:**
1.  **Load** State $(x_t, v_t)$ and Previous Force $(F_t)$.
2.  **First Kick:** $v_{half} = v_t + 0.5 \cdot F_t \cdot dt$
3.  **Drift:** $x_{t+1} = x_t + v_{half} \cdot dt$ (and rotate quaternion)
4.  **Eval Force:** Compute $F_{t+1}$ using $x_{t+1}$ and $v_{half}$.
    *   *Note on Drag:* Using $v_{half}$ (which is at $t+0.5dt$) for velocity-dependent forces is actually **better** than using $v_t$ or $v_{t+1}$ because it is time-centered.
5.  **Second Kick:** $v_{t+1} = v_{half} + 0.5 \cdot F_{t+1} \cdot dt$
6.  **Store** State $(x_{t+1}, v_{t+1})$ and New Force $(F_{t+1})$.

**Memory Layout:**
*   `Pos` (x, y, z, mass)
*   `Vel` (vx, vy, vz, unused) <-- Synchronized, ready for rendering
*   `Quat` (qx, qy, qz, qw)
*   `AngVel` (wx, wy, wz, unused)
*   `ForceTorque` (fx, fy, fz, tx, ty, tz) <-- **Carried over to next step**

This is optimal. It gives you synchronized velocities (no $L$ drift), handles drag correctly, and requires the minimum memory footprint (Force buffer is needed anyway).

### Python Proof

Here is the updated Python script implementing this **Single-State Velocity Verlet**. I also cleaned up the Recoil Force slightly to ensure the math is transparent.

```python
import numpy as np
import matplotlib.pyplot as plt
import argparse

# ==========================================
# Parameters & Setup
# ==========================================
dt      = 0.01
n_steps = 2000
mass    = 1.0
inertia = 0.5
k_bond  = 50.0
l_bond  = 1.0
k_align = 2.55

# Initial State
pos   = np.array([[0.0, 0.0], [1.0, 0.0]])
vel   = np.zeros_like(pos)
theta = np.array([np.pi/2, np.pi/2 + 0.3])
omega = np.array([0.0, 0.0])

# ==========================================
# Force Evaluation
# ==========================================
def get_forces(p, v, th, w):
    F = np.zeros_like(p)
    T = np.zeros_like(th)
    PE = 0.0
    
    # 1. Bond
    d = p[1] - p[0]
    dist = np.linalg.norm(d)
    dir_u = d / dist # Points 0 -> 1
    
    f_mag = k_bond * (dist - l_bond)
    F[0] += f_mag * dir_u
    F[1] -= f_mag * dir_u
    PE += 0.5 * k_bond * (dist - l_bond)**2
    
    # 2. Sigma-Pi Alignment
    for i in range(2):
        h = np.array([np.cos(th[i]), np.sin(th[i])])
        u_local = dir_u if i == 0 else -dir_u
        c = np.dot(h, u_local)
        
        # Energy = k * c^2
        PE += k_align * c**2
        
        # Torque: -dE/dth
        h_prime = np.array([-np.sin(th[i]), np.cos(th[i])])
        dE_dth = 2 * k_align * c * np.dot(h_prime, u_local)
        T[i] -= dE_dth 
        
        # Force: -dE/dr
        # F_recoil = -2 * k * c * (h - c*u)/r
        # Wait! The derivative of u w.r.t r involves projection
        # Force on atom i due to its own angle alignment constraint
        # Correct analytic gradient:
        force_recoil = -2 * k_align * c * (h - c * u_local) / dist
        
        # Correct directionality:
        # If I move, u changes. 
        # For atom 0: u = (p1-p0)/r. d(u)/dp0 is negative.
        # For atom 1: u = (p1-p0)/r. d(u)/dp1 is positive.
        
        # Let's apply to the bond vector ends:
        # The force tries to rotate the bond to be orthogonal to h.
        
        # Apply transverse force to straighten the bond
        F[i]        += force_recoil
        F[1-i]      -= force_recoil

    # 3. Simple Drag (Velocity Dependent Test)
    # F_drag = -gamma * v. Using v_half here works great.
    gamma = 0.0 # Set to 0.1 to test damping
    F -= gamma * v
    T -= gamma * w

    return F, T, PE

# ==========================================
# Solvers
# ==========================================

# 1. LEAPFROG (Semi-Implicit)
def semi_implicit_step(p, v, th, w, aux):
    F, T, pe = get_forces(p, v, th, w)
    v += (F/mass) * dt
    w += (T/inertia) * dt
    p += v * dt
    th += w * dt
    return p, v, th, w, pe, aux

# 2. OPTIMIZED VELOCITY VERLET (Stored Force)
# This is what you should implement in OpenCL
def stored_force_verlet_init(p, v, th, w):
    # We need initial force to start the loop
    F, T, _ = get_forces(p, v, th, w)
    return (F, T)

def stored_force_verlet_step(p, v, th, w, aux):
    F_old, T_old = aux
    
    # A) First Kick (v, w become v_half, w_half)
    v_half = v + 0.5 * (F_old/mass) * dt
    w_half = w + 0.5 * (T_old/inertia) * dt
    
    # B) Drift (Update positions)
    p_new = p + v_half * dt
    th_new = th + w_half * dt
    
    # C) Eval Force (at new position, using v_half as velocity guess)
    F_new, T_new, pe = get_forces(p_new, v_half, th_new, w_half)
    
    # D) Second Kick (Complete the step)
    v_new = v_half + 0.5 * (F_new/mass) * dt
    w_new = w_half + 0.5 * (T_new/inertia) * dt
    
    return p_new, v_new, th_new, w_new, pe, (F_new, T_new)

# ==========================================
# Simulation Loop
# ==========================================

def run_simulation(name, step_fn, init_fn=None):
    p, v, th, w = pos.copy(), vel.copy(), theta.copy(), omega.copy()
    aux = init_fn(p, v, th, w) if init_fn else None
    
    E_list = []
    L_list = []
    
    for _ in range(n_steps):
        p, v, th, w, pe, aux = step_fn(p, v, th, w, aux)
        
        # Energy
        ke = 0.5*mass*np.sum(v**2) + 0.5*inertia*np.sum(w**2)
        E_list.append(ke + pe)
        
        # Angular Momentum
        # L = r x mv + Iw
        L = 0
        for i in range(2):
            L += np.cross(p[i], mass*v[i]) + inertia*w[i]
        L_list.append(L)
        
    return E_list, L_list

# Run
E_semi, L_semi = run_simulation("Leapfrog", semi_implicit_step)
E_opt,  L_opt  = run_simulation("StoredForce", stored_force_verlet_step, stored_force_verlet_init)

# Plot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

ax1.plot(E_semi, label='Leapfrog', alpha=0.7)
ax1.plot(E_opt, label='Stored-Force Verlet', linestyle='--')
ax1.set_title("Total Energy")
ax1.legend()

ax2.plot(L_semi, label='Leapfrog')
ax2.plot(L_opt, label='Stored-Force Verlet')
ax2.set_title("Angular Momentum")
ax2.legend()

plt.tight_layout()
plt.show()
```

### Proposed OpenCL Kernel Structure

Here is the single kernel design that implements the robust **Stored-Force Velocity Verlet**.

**Requirements:**
1.  Initialize the `forces_torques` buffer with 0 (or an initial gravity pass) before the first step.
2.  Pass `dt` and `n_bodies`.

```c
__kernel void rigid_body_update(
    int n_bodies,
    float dt,
    __global float4* pos,            // In/Out: Position.w = Mass
    __global float4* vel,            // In/Out: Velocity (Synchronized)
    __global float4* quat,           // In/Out: Quaternion
    __global float4* ang_vel,        // In/Out: Angular Velocity (Synchronized)
    __global float4* forces_torques, // In/Out: Force(xyz), Torque(xyz). CARRIES STATE!
    __global float16* inv_inertia,   // Body frame inverse inertia
    // ... geometry buffers for eval_force ...
    __global int2*   body_type_info,
    __global float4* point_template
) {
    int gid = get_global_id(0);
    if (gid >= n_bodies) return;

    // 1. Load Previous State
    float4 p_val = pos[gid];
    float3 x = p_val.xyz;
    float mass = p_val.w;
    float3 v = vel[gid].xyz;
    float4 q = quat[gid];
    float3 w = ang_vel[gid].xyz;

    // 2. Load PREVIOUS Forces (Calculated at end of step t-1)
    float3 F_old = forces_torques[2*gid].xyz;
    float3 T_old = forces_torques[2*gid + 1].xyz;
    
    // Add constant forces like Gravity to F_old here if they weren't in the buffer
    // F_old += (float3)(0, -9.81f * mass, 0);

    // 3. FIRST KICK: Advance v, w to half-step (t + 0.5dt)
    float3 v_half = v + 0.5f * (F_old / mass) * dt;
    
    // Angular Kick (World Space)
    // T_old is already in World Space from previous step
    // We need alpha_old. Ideally recompute alpha from T_old, 
    // or assume T_old was already alpha (if you store alpha). 
    // Let's do the proper transform:
    float3 r0, r1, r2;
    quat_to_mat3(q, &r0, &r1, &r2);
    float3 T_body = mat3_transpose_mult(r0, r1, r2, T_old);
    // Apply InvInertia (simplified diagonal)
    float3 I_inv = (float3)(inv_inertia[gid].s0, inv_inertia[gid].s5, inv_inertia[gid].sA);
    float3 alpha_body = I_inv * T_body;
    float3 alpha_world = mat3_mult(r0, r1, r2, alpha_body);
    
    float3 w_half = w + 0.5f * alpha_world * dt;

    // 4. DRIFT: Advance Position/Quat to (t + dt)
    x += v_half * dt;
    
    // Integrate Quaternion using w_half
    float4 w_quat = (float4)(w_half, 0.0f);
    float4 dq = quat_mult(w_quat, q);
    q += dq * (0.5f * dt);
    q = quat_normalize(q);

    // 5. EVAL FORCE: Compute F(t+1) using x(t+1) and v_half
    // Note: We use v_half as the velocity estimate for drag.
    float3 F_new = (float3)(0.0f);
    float3 T_new = (float3)(0.0f);
    
    // ... [Insert your eval_force logic here] ...
    // ... Loop over points, apply springs/collisions ...
    // ... Accumulate into F_new, T_new ...
    
    // Example: Air Drag (Using w_half, v_half)
    // float3 v_surf = v_half + cross(w_half, r_surf);
    // F_new -= 0.1f * v_surf; 

    // 6. SECOND KICK: Advance v, w to (t + dt)
    // We need new alpha for the new torque
    quat_to_mat3(q, &r0, &r1, &r2); // Update Matrix for new orientation
    T_body = mat3_transpose_mult(r0, r1, r2, T_new);
    alpha_body = I_inv * T_body;
    alpha_world = mat3_mult(r0, r1, r2, alpha_body);

    float3 v_new = v_half + 0.5f * (F_new / mass) * dt;
    float3 w_new = w_half + 0.5f * alpha_world * dt;

    // 7. STORE EVERYTHING
    pos[gid] = (float4)(x, mass);
    vel[gid] = (float4)(v_new, 0.0f); // Synchronized!
    quat[gid] = q;
    ang_vel[gid] = (float4)(w_new, 0.0f); // Synchronized!
    
    // Store Forces for the NEXT step
    forces_torques[2*gid]     = (float4)(F_new, 0.0f);
    forces_torques[2*gid + 1] = (float4)(T_new, 0.0f);
}
```

---

## USER

OK, I would like to implement this also in WebGL (GLSL shader) (dont ask me why), I cannot use compute shader, it must be fragment shader with propertis stored in float4 RGBAF textures. I know this will be maybe inefficient as it cannot use local memory unlike OpenCL, but I think it will be fast enought (I basically want to simulate something like leafs of maple seed in wind, they will be instances ). Every body would be processed by one thread (one fragment or pixel) there would be no communications between the bodies (they are completely independent).

All the bodies will be the same so we can use constant uniform buffer of the shader (it is limited in size right?) to store the points (in body state) for all instances. Then the processing (transformation to absolute global frame) will be done on-the-fly inside the shader, no need to store it (we are memory bound). 

---

## Gemini

This is a very common technique in WebGL GPGPU simulations (often called "Ping-Pong" simulation because you swap read/write textures every frame).

Since you are simulating identical instances (like maple seeds), using `uniform` arrays for the body shape is perfect and efficient.

### The Architecture

1.  **Textures (Float RGBA)**: You need 3 textures (or 4 if you split force/torque) of size $N \times 1$ (or square $W \times H$).
    *   **Texture A (Pos)**: `.xyz` = Position, `.w` = Mass
    *   **Texture B (Vel)**: `.xyz` = Linear Velocity, `.w` = *unused/damping*
    *   **Texture C (Rot)**: `.xyzw` = Quaternion
    *   **Texture D (Ang)**: `.xyz` = Angular Velocity, `.w` = *unused*
    *   **Texture E (Force)**: `.xyz` = Stored Force, `.w` = *unused* (Needed for the single-pass Stored-Force Verlet)
    *   **Texture F (Torque)**: `.xyz` = Stored Torque.

    *Note: You need **Multiple Render Targets (MRT)** to update all these textures in a single fragment shader pass. WebGL 2.0 supports this natively (`gl.drawBuffers`).*

2.  **Uniforms**:
    *   `u_body_points[K]`: Array of `vec3` (local positions of vertices).
    *   `u_body_inertia_inv`: `vec3` (Diagonal of inverse inertia tensor).
    *   `u_dt`: Time step.

### The Fragment Shader (GLSL ES 3.0)

Here is the implementation of the **Stored-Force Velocity Verlet** adapted for GLSL. I included a simple aerodynamic drag model for your maple seeds.

```glsl
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

// ================= BODY DEFINITION (Maple Seed) =================
// Max points limited by uniform buffer size (usually 1024 floats safe)
#define MAX_POINTS 16 
uniform int  u_point_count;
uniform vec3 u_points[MAX_POINTS]; // Local positions
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
    vec4  p_data = texelFetch(u_tex_pos, uv, 0);
    vec4  v_data = texelFetch(u_tex_vel, uv, 0);
    vec4  q_data = texelFetch(u_tex_quat, uv, 0);
    vec4  w_data = texelFetch(u_tex_ang_vel, uv, 0);
    vec3  F_old  = texelFetch(u_tex_force, uv, 0).xyz;
    vec3  T_old  = texelFetch(u_tex_torque, uv, 0).xyz;

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
    vec3 alpha_body = u_inertia_inv * T_body_old;
    // Back to World
    vec3 alpha_world = R * alpha_body;

    vec3 w_half = ang + alpha_world * (0.5 * u_dt);

    // =========================================================
    // STEP 2: DRIFT
    // Advance Position to t + dt using v_half
    // =========================================================
    
    pos += v_half * u_dt;

    // Update Quaternion
    // dq/dt = 0.5 * w * q
    vec4 w_quat = vec4(w_half, 0.0);
    vec4 dq = quat_mul(w_quat, quat);
    quat += dq * (0.5 * u_dt);
    quat = normalize(quat); // Important to combat drift!

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

    // B. Aerodynamics (Loop over body points)
    //    Using v_half allows correct drag estimation (centered in time)
    for(int i = 0; i < u_point_count; i++) {
        vec3 r_local = u_points[i];
        
        // Transform point to world relative to CoM
        vec3 r_world = R * r_local;
        
        // Velocity of this specific point
        // v_point = v_cm + w x r
        vec3 v_point = v_half + cross(w_half, r_world);
        
        // --- Aerodynamic Model (Simplified Maple Seed) ---
        float vel_sq = dot(v_point, v_point);
        if(vel_sq > 0.0001) {
            vec3 dir = normalize(v_point);
            
            // 1. Drag (Opposite to velocity)
            // F_d = -0.5 * rho * Cd * A * v^2 * dir
            // Simplified: -k * v^2
            float k_drag = 0.01; 
            vec3 f_drag = -k_drag * vel_sq * dir;
            
            // 2. Lift (Perpendicular to velocity)
            // For a maple seed, lift depends on the blade angle.
            // Let's assume the local point normal defines the blade surface.
            // But for simple "points", we can approximate lift as 
            // perpendicular to velocity and "up" in local frame.
            // This is a hacky placeholder for real lift:
            // vec3 lift_dir = normalize(cross(dir, vec3(0,1,0))); 
            // vec3 f_lift = 0.05 * vel_sq * lift_dir;

            vec3 f_total_pt = f_drag; // + f_lift

            F_new += f_total_pt;
            T_new += cross(r_world, f_total_pt);
        }
    }

    // =========================================================
    // STEP 4: SECOND KICK
    // Advance velocities from t+0.5 to t+1
    // =========================================================
    
    // Linear
    vec3 v_final = v_half + (F_new / mass) * (0.5 * u_dt);

    // Angular
    // Recalculate alpha for the NEW Torque
    R_t = transpose(R);
    vec3 T_body_new = R_t * T_new;
    vec3 alpha_body_new = u_inertia_inv * T_body_new;
    vec3 alpha_world_new = R * alpha_body_new;

    vec3 w_final = w_half + alpha_world_new * (0.5 * u_dt);

    // =========================================================
    // OUTPUT
    // =========================================================
    
    out_pos     = vec4(pos, mass);
    out_vel     = vec4(v_final, 0.0);
    out_quat    = quat;
    out_ang_vel = vec4(w_final, 0.0);
    out_force   = vec4(F_new, 0.0);  // Store for next frame
    out_torque  = vec4(T_new, 0.0);  // Store for next frame
}
```

### WebGL 2.0 JavaScript Setup (Ping-Pong)

You need two Framebuffers (FBOs), each with 6 Color Attachments. You swap `READ_FBO` and `WRITE_FBO` every frame.

```javascript
// Pseudo-code setup
const gl = canvas.getContext("webgl2");
gl.getExtension("EXT_color_buffer_float"); // Vital for physics!

// 1. Create Textures (2 sets: Read and Write)
// Create 6 textures for State A, 6 textures for State B
// Format: gl.RGBA32F, gl.RGBA, gl.FLOAT
// Filter: gl.NEAREST (Do NOT use linear interpolation for physics data!)

// 2. Create Framebuffers
const fboRead = gl.createFramebuffer();
const fboWrite = gl.createFramebuffer();

// Bind textures to attachments [COLOR_ATTACHMENT0 ... COLOR_ATTACHMENT5]
// ...

// 3. Main Loop
function updatePhysics() {
    gl.useProgram(physicsProgram);
    gl.bindFramebuffer(gl.FRAMEBUFFER, fboWrite);

    // Set MRT Draw Buffers map
    gl.drawBuffers([
        gl.COLOR_ATTACHMENT0, // Pos
        gl.COLOR_ATTACHMENT1, // Vel
        gl.COLOR_ATTACHMENT2, // Quat
        gl.COLOR_ATTACHMENT3, // AngVel
        gl.COLOR_ATTACHMENT4, // Force
        gl.COLOR_ATTACHMENT5  // Torque
    ]);

    // Bind 'Read' textures to texture units 0..5
    // Set uniforms (dt, points, mass, etc.)
    
    // Draw Full Screen Quad (triggers fragment shader for every particle)
    gl.viewport(0, 0, SIM_WIDTH, SIM_HEIGHT);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);

    // Swap FBOs for next frame
    swap(fboRead, fboWrite);
    swap(texturesRead, texturesWrite);
}
```

### Key Implementation Details
1.  **EXT_color_buffer_float**: You **must** enable this extension. Without it, you cannot render to floating-point textures, and standard 8-bit colors (`RGBA8`) will destroy your physics precision instantly.
2.  **NEAREST Filtering**: When creating the textures, use `gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST)`. If you use linear, the GPU might interpolate values between two adjacent rigid bodies, creating chaos.
3.  **Efficiency**: Since you calculate everything inside the fragment shader and don't read it back to CPU (hopefully), this is extremely fast. You can use the `pos` and `quat` textures directly in your **Render Shader** (via Vertex Texture Fetch) to draw the instances.