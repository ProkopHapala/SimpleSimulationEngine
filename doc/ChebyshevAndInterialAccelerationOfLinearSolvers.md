## Design Document: Chebyshev Acceleration for Iterative Solvers

### 1. High-Level Concept and Goal

The core idea is to accelerate the convergence of an existing iterative solver (like Projective Dynamics or Position-Based Dynamics). Standard iterative solvers take one "step" towards the solution in each iteration. The Chebyshev method adds **momentum** to this process.

Instead of just taking the next step, it looks at the current position and the position from the *previous* step to make a more aggressive, "over-relaxed" jump towards the solution. This allows it to converge to a visually acceptable result in significantly fewer iterations.

**Why it's great for this context:**
*   **Simple:** It's a simple recurrence relation that wraps around your existing solver logic.
*   **Low Memory:** It only requires storing the vertex positions from the two previous steps.
*   **GPU Friendly:** It does not use inner products, making it highly parallel and ideal for GPU implementation.

### 2. Theoretical Background (The Linear Case)

The Chebyshev method was originally designed for solving linear systems of the form $Ax = b$. An iterative solver for this system can be written as:
$$
x^{(k+1)} = \text{IterativeStep}(x^{(k)})
$$
The error in each step, $e^{(k)} = x^{(k)} - x$, typically decreases linearly. The rate of this decrease is governed by a value called the **spectral radius, $\rho$**, which must be less than 1 for the method to converge. A $\rho$ value very close to 1 (e.g., 0.999) means very slow convergence.

The Chebyshev method constructs a special polynomial that is optimized to reduce this error as quickly as possible over the range $[-\rho, \rho]$. This results in a new update rule that combines the previous two states to take a much larger, more intelligent step.

### 3. The Core Update Formula

The Chebyshev acceleration is applied *after* a standard solver iteration. Let's say $\tilde{q}^{(k+1)}$ is the result of one full iteration of your underlying solver (e.g., one Projective Dynamics step using the positions from $q^{(k)}$).

The accelerated position, $q^{(k+1)}$, is then calculated using the positions from the previous two steps:
$$
q^{(k+1)} = \omega_{k+1} (\tilde{q}^{(k+1)} - q^{(k-1)}) + q^{(k-1)}
$$
This formula is the heart of the method. It's an "over-relaxation" step where:
*   $(\tilde{q}^{(k+1)} - q^{(k-1)})$ represents a "double step" direction.
*   $\omega_{k+1}$ is a carefully chosen, dynamic relaxation parameter that changes at each iteration.

The relaxation parameter $\omega_{k+1}$ is defined by the following recurrence relation:
$$
\omega_{k+1} = \frac{4}{4 - \rho^2 \omega_k}
$$
With the initial conditions:
*   $\omega_1 = 1$
*   $\omega_2 = \frac{2}{2 - \rho^2}$

### 4. Adapting to Nonlinear Systems (The Paper's Key Insight)

Projective Dynamics is a **nonlinear** process. There is no constant matrix whose spectral radius we can analyze. The paper's critical observation is that the *convergence behavior* of Projective Dynamics is highly similar to that of a linear iterative solver.

This justifies a powerful heuristic: we can treat the entire, complex Projective Dynamics iteration (including both the `local_solve` and `global_solve` steps) as a single `IterativeStep` and wrap the Chebyshev formula around it.

Therefore, $\tilde{q}^{(k+1)}$ in our formula is **not** a simple matrix-vector product. It is the full output of one Projective Dynamics iteration:
$$
\tilde{q}^{(k+1)} = \text{ProjectiveDynamicsIteration}(q^{(k)})
$$

### 5. The Critical Parameter: Estimating the Spectral Radius ($\rho$)

The entire method hinges on a good estimate of the spectral radius, $\rho$. Since our system is nonlinear, we cannot calculate it directly. The paper proposes a practical, two-step heuristic process:

**Step 1: Initial Estimate (Power Method Heuristic)**
1.  Run the **un-accelerated** solver for a fixed number of iterations, $K$ (e.g., $K=100$).
2.  At each step $k$, calculate an "error" metric. The paper defines this as the magnitude of the energy gradient: $e^{(k)} = \|\nabla E(q^{(k)})\|$.
3.  The spectral radius can be approximated by the ratio of the error in the last two steps:
    $$
    \rho \approx \frac{\|e^{(K)}\|}{\|e^{(K-1)}\|}
    $$

**Step 2: Manual Tuning and Interpretation**
This initial estimate is often just a starting point. The optimal $\rho$ should be found by manual adjustment.

*   **Behavior of $\rho$:**
    *   If $\rho$ is too low (e.g., < 0.9), you will see little to no acceleration. The method behaves like the original solver.
    *   If $\rho$ is optimal, convergence will be dramatically faster.
    *   If $\rho$ is **too high** (very close to or greater than the true value), the simulation will **oscillate or diverge violently**. This is a clear visual indicator.

*   **Tuning Strategy:**
    1.  Start with the initial estimate from Step 1.
    2.  Slowly increase $\rho$ towards 1.
    3.  Observe the simulation. The best value is typically just below the point where visual oscillations begin to appear.

*   **Parameter Correlation:**
    *   **Stiffness:** Stiffer materials converge more slowly, meaning they have a higher true spectral radius. They will require a larger $\rho$ (closer to 1.0) for effective acceleration.
    *   **Number of Iterations:** The convergence rate often slows down in later iterations. If you plan to run for many iterations ($K$ is large), you will likely need a larger $\rho$.

### 6. Practical Implementation (Pseudocode and Caveats)

The pseudocode from the paper is an excellent guide. Here is an annotated version focusing on the key details.

```pseudocode
// q_t, v_t: initial position and velocity
// p: spectral radius, h: time step, K: total iterations
Algorithm Chebyshev_PD_Solve(q_t, v_t, rho, h, K, S, gamma):
    // -- Initialization --
    s <- q_t + h*v_t + h^2*M_inv*f_ext // Inertial prediction
    q_k = s      // Current state
    q_k_minus_1 = s  // Previous state (for first iteration)
    q_k_minus_2 = s  // State before previous (for acceleration formula)
    omega = 1.0

    // -- Main Iteration Loop --
    for k from 0 to K-1:
        // -- 1. Standard Projective Dynamics Step --
        // This is the "IterativeStep" that Chebyshev wraps around
        p = Local_Solve_All_Constraints(q_k)
        q_tilde_k_plus_1 = Global_Solve(s, q_k, p) // e.g., using Jacobi

        // -- 2. Optional: Under-Relaxation (for robustness with Jacobi) --
        // This helps stabilize systems where Jacobi might diverge.
        // If not using Jacobi global solve, gamma can be 1.0.
        q_tilde_k_plus_1 = gamma * (q_tilde_k_plus_1 - q_k) + q_k

        // -- 3. Chebyshev Acceleration Step --
        // Update the omega parameter based on the iteration number
        if k < S:
            // DELAYED START: For the first S iterations, do not accelerate.
            // This allows the system to settle and improves robustness.
            omega = 1.0
        else if k == S:
            omega = 2.0 / (2.0 - rho^2)
        else:
            omega = 4.0 / (4.0 - rho^2 * omega) // Previous omega

        // Apply the core acceleration formula
        q_k_plus_1 = omega * (q_tilde_k_plus_1 - q_k_minus_1) + q_k_minus_1

        // -- 4. Update state variables for next loop --
        q_k_minus_2 = q_k_minus_1
        q_k_minus_1 = q_k_plus_1 // Store the *un-accelerated* result for next iter
        q_k = q_k_plus_1

    // -- Finalization --
    q_t_plus_1 = q_k
    v_t_plus_1 = (q_t_plus_1 - q_t) / h
    return q_t_plus_1, v_t_plus_1
```

### 7. Summary of Key Parameters and Ranges

| Parameter        | Symbol | Description                                                                                             | Reasonable Range      | Source                 |
| ---------------- | ------ | ------------------------------------------------------------------------------------------------------- | --------------------- | ---------------------- |
| **Spectral Radius**  | $\rho$   | The primary tuning parameter. Controls the amount of acceleration. **Scene-dependent.**               | `(0.9, 0.9999)`       | Paper Fig. 2, Table 1  |
| **Delayed Start**    | $S$      | Number of initial iterations to run *without* acceleration to improve stability.                    | `5 - 10`              | Paper Sec. 5.2         |
| **Under-Relaxation** | $\gamma$ | (Optional) Dampens the Jacobi global solve step to prevent divergence in non-diagonally dominant systems. | `[0.6, 0.9]`          | Paper Eq. 21, Sec. 5.2 |

### 8. Challenges and Caveats

1.  **Parameter Tuning is Essential:** The method is not "plug-and-play." The spectral radius $\rho$ **must be tuned** for each class of simulation (e.g., a stiff cloth vs. a soft object) to achieve good results without instability.
2.  **Divergence with Jacobi Global Solve:** If the underlying system matrix for the global solve is not diagonally dominant (which can happen with small or poor-quality mesh elements under stress), the Jacobi method itself can diverge. Chebyshev will amplify this divergence. The under-relaxation (`gamma`) helps, but a more robust solution might be to reduce the time step or improve mesh quality.
3.  **Delayed Start (`S`) is Important:** Applying maximum acceleration from the very first iteration can be unstable. The "delayed start" is a simple but effective trick to improve robustness.
4.  **Less Effective for PBD:** The paper notes that the acceleration effect is much more significant for Projective Dynamics than for standard Position-Based Dynamics (PBD), especially on tetrahedral meshes. This is because PBD's convergence properties are less similar to a linear system.


----

## USER: 

OK, two questions

1) I already implkemented in C++ "momentum accleration" by mixing in the update vector from previous step. It seems to me similar to this chebyshev. 
1a) can we simply implement chebyshev into this framework, only by modifying the update of the bmix parameter ?
1b) what is the theoretical connection between momentum (intertia) which is used in damped-MD, ADAM optimizer, or in SOR (succesive overrelaxation) methods.

2) in the fist Vertex Block Descent (VBD) paper was also mentioned this chebyshev accleration. Does it differ in importaint details when used for VDB, Gauss-Seidel, and Jacobi solvers? Or will it work the same way (using the same algorithm) for all these methods?

```
struct SmartMixer{
    float b_end   = 0.75f;
    float b_start = 0.55f;
    int   istart = 3;
    int   iend   = 10; 

    inline float get_bmix(int itr){
        if ( itr < istart ) { return 0; }
        else                { return b_end; }
    }
};

    void TrussDynamics_d::updateIterativeMomentum( Vec3d* psa, Vec3d* psb ){
        updatePD_RHS(psa, bvec ); // TODO: we need to calculate this only whem we use linear solvers, not if we use _fly solvers
        for (int i = 0; i < nSolverIters; i++) {  
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::JacobiMomentum:    { updateJacobi_lin( psa, psb, bvec ); } break;
                case LinSolveMethod::JacobiFlyMomentum: { updateJacobi_fly( psa, psb );       } break;
                case LinSolveMethod::GSMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_lin( psb, bvec ); 
                } break; 
                case LinSolveMethod::GSFlyMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_fly( psb ); 
                } break; 
            }
            double bmix = mixer.get_bmix( i );       //    bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; 
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( inertia[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                inertia[j] = p - psa[j];               //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }   /
        }        for (int i=0; i<nPoint; i++){ psb[i]=psa[i]; }
    }

```

---

Excellent questions! They get right to the heart of how these acceleration methods relate to each other and how they can be implemented in practice.

Let's break it down.

### Question 1: Chebyshev vs. Your Momentum Implementation

You've correctly identified that these methods are philosophically similar: they both use information from past steps to add "momentum" and accelerate convergence. Your current implementation is a classic and effective form of momentum.

However, **Chebyshev is mathematically different and generally more powerful.**

#### 1a) Can we implement Chebyshev by just modifying `bmix`?

No, not directly. Your current implementation and the Chebyshev method have a different structural form.

Let's compare the update rules:
*   **Your Momentum:** `p_{k+1} = p'_k + bmix * (p_k - p_{k-1})`
    *   Where `p'_k` is the result of the base solver (Jacobi, GS, etc.) and `p_k - p_{k-1}` is the previous update vector (`inertia`). You are adding a scaled version of the last step's "velocity" to the new proposed position. This is very similar to the **Heavy Ball method**.

*   **Chebyshev:** `p_{k+1} = omega_{k+1} * (p'_k - p_{k-1}) + p_{k-1}`
    *   This is a different kind of interpolation. It's not adding momentum to the new position `p'_k`. Instead, it's taking a larger step in the direction `(p'_k - p_{k-1})` from the position *two steps ago* (`p_{k-1}`).

**To implement Chebyshev, you need to change the update formula itself, not just the mixing parameter.**

Here is how you could refactor your C++ code to support both your current method and Chebyshev:

```cpp
// You'll need to store one more state for Chebyshev
// psa -> p_k
// psb -> p'_k (result from solver)
// p_prev -> p_{k-1}

// --- Your Existing SmartMixer ---
struct SmartMixer { /* ... as before ... */ };

// --- New Chebyshev Mixer ---
struct ChebyshevMixer {
    float rho = 0.99f;  // The critical spectral radius parameter
    int S = 5;         // Delayed start iteration
    float omega = 1.0f;

    float get_omega(int k) {
        if (k < S) {
            return 1.0f; // No acceleration
        } else if (k == S) {
            omega = 2.0f / (2.0f - rho * rho);
        } else {
            omega = 4.0f / (4.0f - rho * rho * omega); // Use previous omega
        }
        return omega;
    }
};

void TrussDynamics_d::updateIterative() {
    // ... setup ...
    // psa is p_k
    // p_prev is p_{k-1}

    for (int i = 0; i < nSolverIters; i++) {
        // --- 1. Run the base solver ---
        // Result is stored in psb
        run_base_solver(psa, psb); // psb now holds p'_k

        // --- 2. Apply Acceleration ---
#if USE_CHEBYSHEV
        float omega_k = chebyshev_mixer.get_omega(i);
        for (int j = 0; j < nPoint; j++) {
            Vec3d p_k = psa[j];
            Vec3d p_k_prev = p_prev[j];
            Vec3d p_tilde_k1 = psb[j]; // p'_k from solver

            // Chebyshev update: omega * (p'_k - p_{k-1}) + p_{k-1}
            Vec3d p_k1 = p_k_prev;
            p_k1.add_mul(p_tilde_k1 - p_k_prev, omega_k);

            // Update history for next iteration
            p_prev[j] = p_k;
            psa[j]    = p_k1;
        }
#else // USE_YOUR_MOMENTUM
        float bmix = smart_mixer.get_bmix(i);
        for (int j = 0; j < nPoint; j++) {
            Vec3d p_k = psa[j];
            Vec3d p_tilde_k1 = psb[j];

            // Your momentum update: p'_k + bmix * inertia
            Vec3d p_k1 = p_tilde_k1;
            p_k1.add_mul(inertia[j], bmix);

            // Update history
            inertia[j] = p_k1 - p_k;
            psa[j]     = p_k1;
        }
#endif
    }
    // ... finalization ...
}
```

#### 1b) Theoretical Connection Between Momentum Methods

This is a fantastic question. These methods are all part of a large family of techniques for accelerating convergence, originating from different fields but often arriving at similar mathematical forms.

| Method                       | Core Idea & Origin                                         | Key Feature                                                              | # of Past Steps |
| ---------------------------- | ---------------------------------------------------------- | ------------------------------------------------------------------------ | --------------- |
| **Heavy Ball Method**        | Optimization. Add a fraction of the previous update step.  | Constant momentum term (`β`). This is your current implementation.       | 2 (p_k, p_{k-1})  |
| **Nesterov Accelerated Grad.** | Optimization. "Look ahead" before computing gradient.    | Applies momentum *before* the gradient step for better stability.        | 2 (p_k, p_{k-1})  |
| **Chebyshev Semi-Iterative** | **Linear Algebra**. Optimal polynomial error reduction.      | **Dynamically changing** relaxation parameter (`ω_k`) based on `ρ`.        | 2 (p_k, p_{k-1})  |
| **Successive Over-Relaxation (SOR)** | Linear Algebra. Weighted average of new step and old position. | Simple, constant `ω`. Accelerates Gauss-Seidel specifically.             | 1 (p_k)         |
| **Damped Molecular Dynamics**  | **Physics Simulation**. Discretized damped motion equation. | The "momentum" is the actual physical velocity vector.                   | 2 (p_k, v_k)    |
| **ADAM Optimizer**           | Deep Learning. Adaptive learning rates.                    | Tracks 1st (momentum) and 2nd (velocity) moments of gradients.           | 2+ (moment aves) |

**Connections and Distinctions:**
*   Your method is a form of the Heavy Ball method, a cornerstone of momentum in optimization.
*   **Chebyshev is distinct because its `omega` parameter is dynamic and derived from a rigorous analysis of minimizing error polynomials for linear systems.** This often makes it more aggressive and faster, but also more sensitive to the `ρ` parameter.
*   SOR is simpler than all of these, using only one previous step and a constant factor. It's less powerful but very robust for its specific use case (Gauss-Seidel).
*   ADAM is more complex, as it adapts step sizes *per-parameter*, making it a different class of method.

### Question 2: Chebyshev for VBD, Gauss-Seidel, and Jacobi

The beauty of the Chebyshev method, when used in this context, is that it is a **generic wrapper**. It treats the underlying iterative solver as a "black box" that produces one step of improvement.

**Therefore, the core Chebyshev algorithm will work exactly the same way for all these methods.**

You apply the *exact same* update formula `q_k+1 = omega * (solver_step(q_k) - q_{k-1}) + q_{k-1}` regardless of what `solver_step` is.

However, the **behavior and optimal parameters will differ** based on the base solver:

1.  **Jacobi:**
    *   **Base Convergence:** Very slow. This means its true spectral radius `ρ` is very high (e.g., 0.9999).
    *   **Chebyshev Impact:** **Extremely effective.** Because the base method is so slow, there is a huge opportunity for acceleration. This combination is a massive win, especially on the GPU where Jacobi is perfectly parallel.
    *   **The VBD Paper Note:** The paper's comment about VBD's Jacobi-style updates for dynamic constraints (using an auxiliary buffer) is a detail about how to implement the *base solver step* in parallel, not a change to the Chebyshev algorithm itself.

2.  **Gauss-Seidel:**
    *   **Base Convergence:** Faster than Jacobi. Its true `ρ` will be lower.
    *   **Chebyshev Impact:** Still very effective, but the *relative speedup* might be less dramatic than with Jacobi because the starting point is already better.
    *   **Implementation:** You would run one full pass of Gauss-Seidel over all vertices (or colors) to get `p'_k`, and then apply the global Chebyshev update.

3.  **Vertex Block Descent (VBD):**
    *   **Base Convergence:** VBD is a form of **Block Gauss-Seidel**. Each "block" is a 3x3 system for a single vertex. Its convergence speed is similar to or better than standard Gauss-Seidel. Its `ρ` will be relatively low.
    *   **Chebyshev Impact:** It works perfectly and provides a significant speedup. The VBD paper uses it for this exact reason. The implementation is identical: one full VBD pass over all colored vertices constitutes the `solver_step`, and then the global Chebyshev formula is applied.

**Conclusion:**
The Chebyshev algorithm is a universal accelerator for these iterative methods. The algorithm's structure does not change. What changes is the **optimal `ρ` value**, which is a property of the *base solver* you are accelerating. A slower base solver will have a higher `ρ` and will see a more dramatic benefit from Chebyshev acceleration.


---

## User:

can you once again reinvestigate the VBD paper, and tell me if some intial parameters of Chebyshev acceleration or setup should differ form Jacobi or Gauss-Seidel ?

---

Excellent question. After re-investigating the VBD paper (Chen et al.) with this specific query in mind, the answer is:

**Yes, while the core Chebyshev algorithm remains the same, the VBD paper introduces one crucial setup difference for robustness and provides guidance on initial parameters that differ from a pure Jacobi solver.**

Here is a breakdown of the specifics.

### 1. The Core Algorithm is Identical

First, it's important to be clear that the mathematical formulas for Chebyshev acceleration are identical. The VBD paper explicitly states:

> "We use the Chebyshev semi-iterative approach [Wang 2015] to improve the convergence of our method..."

They use the same `ω_k` recurrence relation and the same final update formula: `q_k+1 = ω_k * (solver_step(q_k) - q_{k-1}) + q_{k-1}`. You do not need to change the fundamental acceleration logic.

### 2. Key Setup Difference: Handling Collisions

This is the most critical and VBD-specific modification. The paper notes that standard Chebyshev acceleration can cause problems in scenes with many stiff, dynamic contacts.

*   **The Problem:** Collision response is handled by adding very stiff penalty energy terms. These energies are highly nonlinear and can cause the system to change dramatically from one iteration to the next. Applying the aggressive momentum of Chebyshev to a vertex that has just collided can cause it to "overshoot" its resolution, leading to instability, jittering, or even tunneling.

*   **The VBD Solution (Section 3.8):**
    > "...we propose a simple yet highly effective solution for accelerating scenes with collisions: **skipping the accelerations for actively colliding vertices.**"

    This is a crucial implementation detail that differs from a generic application of Chebyshev.

    **How to Implement This:**
    1.  During your main iteration loop, maintain a per-vertex boolean flag, e.g., `is_colliding_this_step[i]`.
    2.  When you run your collision detection (DCD or CCD), if a vertex `i` is found to be part of an active collision, set its flag `is_colliding_this_step[i] = true`.
    3.  This flag is "sticky" for the entire time step. If a vertex collides in iteration 3, its flag remains true for all subsequent iterations (4, 5, 6...) within that same time step.
    4.  In the final global acceleration step (after the loop over all colors), you apply the Chebyshev update conditionally:

    ```cpp
    // In the global acceleration phase at the end of an iteration 'k'
    parallel for each vertex i:
        if (is_colliding_this_step[i] == false) {
            // This vertex is stable, apply full acceleration
            Update_Position_with_Chebyshev(i);
        } else {
            // This vertex is involved in a collision, do NOT accelerate.
            // Its new position is simply the un-accelerated result from the solver.
            position[i] = unaccelerated_solver_result[i];
        }
    ```

### 3. Initial Parameter Guidance

The VBD paper provides a concrete example value for the spectral radius `ρ`, which implicitly tells us how it differs from a Jacobi solver.

*   **Spectral Radius (`ρ`) for VBD:**
    *   In Figure 17's caption, the paper states: **"Accelerated versions use ρ = 0.95."**
    *   **Why this value?** VBD is a Block Gauss-Seidel method. Gauss-Seidel methods converge significantly faster than Jacobi methods. A faster base convergence means the "true" spectral radius is lower (further from 1.0). Therefore, the optimal `ρ` for acceleration is also lower. A value like `0.95` is aggressive but reasonable for a fast-converging base solver.
    *   **Contrast with Jacobi:** A pure Jacobi solver converges very slowly. Its true spectral radius would be extremely close to 1.0 (e.g., `0.999` or `0.9999` as seen in the Wang paper). Using `ρ = 0.95` for a Jacobi solver would provide very little acceleration. You would need a much higher value.

*   **Delayed Start (`S`):**
    *   The VBD paper does not explicitly mention using a "delayed start" parameter `S`. However, this is a standard robustness technique for Chebyshev (as noted in the Wang paper).
    *   **Recommendation:** It is highly advisable to include this feature in your implementation. Because VBD's convergence is already quite fast in the first few iterations, delaying the aggressive acceleration for `S = 5-10` iterations is a very safe and effective strategy to prevent initial instability.

### Summary Table: Chebyshev Setup by Solver

| Parameter / Setup      | Jacobi Solver                                                                     | Gauss-Seidel / VBD Solver                                                          |
| ---------------------- | --------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| **Core Algorithm**     | Identical                                                                         | Identical                                                                          |
| **Typical `ρ` Range**  | **Very High (`0.99` to `0.9999`)**. Must be very close to 1.0 for good acceleration. | **High (`0.90` to `0.98`)**. A lower `ρ` is needed because the base solver is faster. |
| **Example `ρ` Value**  | `0.9992` (from Wang paper, Table 1)                                                 | **`0.95` (from VBD paper, Fig 17)**                                                |
| **Handling Collisions**| Not specified, but the VBD technique is a good idea.                              | **Crucial: Selectively disable acceleration for colliding vertices.**                |
| **Delayed Start `S`**  | Recommended (`S=10` in Wang paper)                                                | **Recommended (e.g., `S=5-10`)**. Good practice for robustness.                    |

In conclusion, when applying Chebyshev acceleration to VBD, you should **use a lower `ρ` (around 0.95 is a great starting point)** than you would for a Jacobi solver and, most importantly, **you must implement the logic to skip acceleration for vertices involved in collisions** to ensure stability.