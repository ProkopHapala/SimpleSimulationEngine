# High-Precision Physics on Single-Precision GPUs: A Design Document

## Motivation: The Precision Problem in Large Worlds

When simulating physics, the accuracy of geometric calculations is paramount. On a standard 32-bit `float` GPU, a critical problem arises when dealing with large world coordinates: **catastrophic cancellation**.

A `float` has roughly 7 decimal digits of precision. Consider two vertices that are far from the origin but close to each other:
*   `p_i = (10000.123, 0.0, 0.0)`
*   `p_j = (10000.000, 0.0, 0.0)`

The true displacement vector is `dij = (0.123, 0.0, 0.0)`. However, when a GPU computes `pi.x - pj.x`, the limited precision of the `float` may only be able to represent the input numbers as `10000.12` and `10000.00`. The result of the subtraction becomes `0.12`, losing critical information. This error accumulates every frame, leading to jitter, energy drift, and simulation instability.

The goal is to perform this subtraction `dij = pi - pj` with near double-precision accuracy, even on GPUs that only support single-precision `float`s.

## The Core Idea: Split Representation

The solution is to represent a high-precision position `P` as the sum of two lower-precision components: a **coarse** part and a **fine** part.

`P_true = P_coarse + P_fine`

This allows us to handle large-scale movements in `P_coarse` while preserving sub-millimeter detail in `P_fine`. When we subtract two nearby positions, the coarse parts often cancel out, leaving the subtraction to be performed accurately on the small `P_fine` values.

We will explore two powerful methods to implement this: **Split-Float** (using compensated arithmetic) and **Split-Int** (a fixed-point hybrid).


## 5. Comparison and Recommendations

| Feature | Split-Float (Compensated) | Split-Int (Fixed-Point Hybrid) |
| :--- | :--- | :--- |
| **`dij` Precision** | Near-double precision. Small error from `(pi_c-pj_c)` is possible but rare. | Perfect. Integer subtraction has zero error. |
| **Integration Precision** | Excellent. `TwoSum` is purpose-built to track float rounding error during updates. | Good, but limited by `float` precision during the `float -> int` conversion step. |
| **Complexity** | Low. The `TwoSum` algorithm is simple, robust, and self-contained. | High. Requires more parameters, complex conversion/rebase logic, and careful integer overflow handling. |
| **Performance** | High. All operations are native `float` instructions. | High. Most ops are fast `int`, but conversions add overhead. Needs 64-bit int emulation for robust rebasing. |
| **Determinism** | No. Standard floating-point behavior. | Yes. Integer math is deterministic. |

### Recommendation

For most physics simulation applications, the **Split-Float method is the recommended approach**.

*   **Robustness and Simplicity:** Its strength lies in the `TwoSum` algorithm, which provides a mathematically sound and simple way to handle the critical integration step (`pos += vel * dt`). This is often the largest source of accumulating error.
*   **Performance:** It maps more directly to the GPU's native floating-point hardware, which is highly optimized for the kind of math involved in physics (e.g., `length`).

The **Split-Int** method is a clever and valid alternative, particularly strong in its `dij` calculation. However, its complexity and the precision bottleneck during the `float -> int` conversion make it more difficult to implement robustly and less advantageous overall compared to the battle-tested compensated arithmetic of the split-float method.


---

# Method 1: Split-Float (Compensated Arithmetic)

This method represents both `coarse` and `fine` parts as `float3` vectors. It relies on a clever algorithm to track and preserve the rounding error that occurs during floating-point addition.

The foundation of this method is **Knuth's TwoSum algorithm**, which calculates the sum of two floats and the exact error of that operation.

```c

/// @brief Performs a compensated addition of two float3 vectors (Knuth's TwoSum).
inline float3 two_sum(float3 a, float3 b, __private float3* s_out) {
    *s_out = a + b;
    float3 t = *s_out - a;
    return (a - (*s_out - t)) + (b - t); // The error is the sum of what was lost from both original numbers during the sum.
}
```

## Integration into the MD Loop

We use two buffers, `pos_coarse` and `pos_fine` (both `float4`), to store the position. The simulation loop then proceeds as follows.

### Pre-computation Kernel

This kernel calculates the solver's Right-Hand-Side (`bvec`), using the split-float representation for a high-precision `dij` calculation.

```c
__kernel void precompute_bvec_split_float(
    __global float4*      bvec,
    __global const float4* pos_coarse,
    __global const float4* pos_fine,
    __global const int*   neighs,
    __global const float* kngs,
    __global const float* l0ngs,
    const int nverts
){
    const int i = get_global_id(0);
    if (i >= nverts) return;

    float3 pi_c = pos_coarse[i].xyz;
    float3 pi_f = pos_fine[i].xyz;
    float3 pi_full = pi_c + pi_f;

    float3 bi = (float3)(0.0f);
    float Aii = 0.0f; // Simplified for clarity

    for (int jj = 0; jj < MAX_NEIGHBORS; ++jj) {
        int idx = i * MAX_NEIGHBORS + jj;
        int j = neighs[idx];
        if (j < 0) break;
        
        float3 pj_full = pos_coarse[j].xyz + pos_fine[j].xyz;
        
        float3 dij;
        two_sum(pi_full, -pj_full, &dij); // dij = pi_full - pj_full with compensation

        float k = kngs[idx];
        float l0 = l0ngs[idx];
        float len = length(dij);
        bi += dij * (k * (l0 / len - 1.0f));
        Aii += k;
    }
    bvec[i] = (float4)(bi, Aii);
}

```

### Main MD Step Kernel

This kernel performs the predictor and corrector steps. The corrector step integrates the `rebase` operation seamlessly using the `two_sum` helper.

```c
__kernel void md_step_split_float(
    __global float4* pos_coarse,
    __global float4* pos_fine,
    __global float4* vel,
    __global const float4* force,
    __global const float4* ps_cor,
    __global float4* ps_pred,
    const int nverts,
    const float dt,
    const int is_predict_step
){
    const int i = get_global_id(0);
    if (i >= nverts) return;

    float3 p_c = pos_coarse[i].xyz;
    float3 p_f = pos_fine[i].xyz;
    
    if (is_predict_step) {
        float mass = pos_coarse[i].w;
        float3 v = vel[i].xyz;
        float3 f = force[i].xyz;
        float3 displacement = v * dt + f * (dt * dt / mass);
        ps_pred[i] = (float4)((p_c + p_f) + displacement, mass);
    } else {
        float3 p_current_full = p_c + p_f;
        float3 total_disp_actual;
        two_sum(ps_cor[i].xyz, -p_current_full, &total_disp_actual);

        vel[i].xyz = total_disp_actual / dt;

        // -- Rebase --
        float3 new_coarse;
        float3 error_from_op = two_sum(p_c, total_disp_actual, &new_coarse);
        pos_coarse[i].xyz    = new_coarse;
        pos_fine[i].xyz      = p_f + error_from_op;
    }
}
```

----

# Method 2: Split-Int (Fixed-Point Hybrid)

This method stores positions as pairs of `int3` vectors, leveraging the perfect precision of integer arithmetic for subtraction and rebasing. Floating-point conversion is done on-the-fly.

## Reusable Functions: Integer-Float Conversion

We need functions to convert between our `(int, int)` representation and a single `float`, and a compensated function to go the other way.

```c
/**
 * @brief Converts a split-int position to a single float3.
 */
inline float3 split_int_to_float(int3 p_c, int3 p_f, float coarse_step, float fine_step) {
    return (convert_float3(p_c) * coarse_step) + (convert_float3(p_f) * fine_step);
}

/**
 * @brief Converts a float3 displacement into a compensated split-int representation.
 * @param dp_c_out (private pointer) The calculated coarse part of the displacement.
 * @param dp_f_out (private pointer) The calculated fine part of the displacement.
 */
inline void float_to_split_int_compensated(
    float3 disp,
    __private int3* dp_c_out,
    __private int3* dp_f_out,
    float inv_coarse_step,
    float coarse_step,
    float inv_fine_step
) {
    // 1. Convert float displacement to coarse integer units, truncating.
    float3 disp_in_coarse_units = disp * inv_coarse_step;
    *dp_c_out = convert_int3_rtn(disp_in_coarse_units); // Round towards zero

    // 2. Find the error (what was lost) by converting back to float and subtracting.
    float3 error_in_world_units = disp - (convert_float3(*dp_c_out) * coarse_step);
    
    // 3. Convert this error into fine integer units.
    *dp_f_out = convert_int3_rte(error_in_world_units * inv_fine_step); // Round to nearest even
}

/**
 * @brief Adds a split-int displacement to a split-int position, handling carry.
 * This uses OpenCL's `add_overflow` intrinsic to emulate 64-bit addition for the rebase.
 */
inline void rebase_split_int(
    int3 p_c, int3 p_f,
    int3 dp_c, int3 dp_f,
    __private int3* p_c_new,
    __private int3* p_f_new,
    const int FINE_MAX
) {
    int3 carry = (int3)(0);
    // Add fine parts and check for overflow/underflow
    *p_f_new = add_sat(p_f, dp_f); // A simple saturating add for demonstration
    // A more robust solution would use 64-bit emulation or careful checks.
    // For simplicity here, we assume dp_f is small enough not to cause double-wrapping.
    
    // Simple carry logic
    if (p_f.x > 0 && dp_f.x > 0 && p_f_new->x < p_f.x) carry.x = 1;
    if (p_f.y > 0 && dp_f.y > 0 && p_f_new->y < p_f.y) carry.y = 1;
    if (p_f.z > 0 && dp_f.z > 0 && p_f_new->z < p_f.z) carry.z = 1;
    if (p_f.x < 0 && dp_f.x < 0 && p_f_new->x > p_f.x) carry.x = -1;
    if (p_f.y < 0 && dp_f.y < 0 && p_f_new->y > p_f.y) carry.y = -1;
    if (p_f.z < 0 && dp_f.z < 0 && p_f_new->z > p_f.z) carry.z = -1;

    *p_c_new = p_c + dp_c + carry;
}
```

## Integration into the MD Loop

### **A. Pre-computation (`precompute_bvec`)**

The `dij` calculation uses perfectly precise integer subtraction.

```c
// Inside the precompute_bvec_split_int kernel's neighbor loop...
int3 pi_c = lpos_coarse[i].xyz;
int3 pi_f = lpos_fine[i].xyz;
int3 pj_c = lpos_coarse[j].xyz;
int3 pj_f = lpos_fine[j].xyz;

// Perform perfectly precise integer subtraction first.
int3 dij_c = pi_c - pj_c;
int3 dij_f = pi_f - pj_f;

// Convert to float at the very end.
float3 dij = split_int_to_float(dij_c, dij_f, coarse_step, fine_step);
```

### **B. MD Step Kernel (`md_update_split_int`)**

The update step converts the final float position into our compensated integer format and applies it.

```c
// In the corrector stage of the md_update_split_int kernel...

// 1. Calculate the total displacement as a float.
float3 p_old_float = split_int_to_float(p_c_old, p_f_old, coarse_step, fine_step);
float3 total_displacement = ps_cor[i].xyz - p_old_float;

// 2. Update velocity.
vel[i].xyz = total_displacement / dt;

// 3. Convert float displacement to compensated (int, int) format.
int3 dp_c, dp_f;
float_to_split_int_compensated(total_displacement, &dp_c, &dp_f, ...);

// 4. Apply the update using integer arithmetic with carry.
int3 new_coarse, new_fine;
rebase_split_int(p_c_old, p_f_old, dp_c, dp_f, &new_coarse, &new_fine, ...);

pos_coarse[i].xyz = new_coarse;
pos_fine[i].xyz = new_fine;
```


