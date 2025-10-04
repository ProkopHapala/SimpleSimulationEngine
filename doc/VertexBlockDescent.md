# Vertex Block Descent (VBD)

## Purpose

This document rewrites article [Vertex Block Descent, A.Chen,(2024)](https://dl.acm.org/doi/10.1145/3658179) focusing on the main VBD optimization equations and Algorithm 1 into Markdown format (with inline LaTeX). It collects implementation-oriented details for a C++ / OpenCL implementation, focusing on:

* Global optimization formulation (Sec. 3.1)
* Local system solver and Newton step (Sec. 3.2)
* Collisions (Sec. 3.5)
* Initialization (Sec. 3.7)
* Accelerated iterations (Sec. 3.8)
* Parallelization (Sec. 3.9)
* Algorithm 1 (pseudocode for one time step)

---

## Design notes: modular implicit dynamics in `Truss`

- **`Truss.run_dynamics()` outer loop**
  - Encapsulate implicit Euler stepping over `nsteps`, updating `(x, v)` using shared routines for inertia/anchor handling and force accumulation.
  - Accepts solver handle (callable or object) operating on current state; solver returns updated positions after solving local system.
  - Maintains diagnostics/trajectory logging flags; keeps plotting separate per global repo rules.

- **Solver interface** (Python prototype)
  ```python
  def solve_step(truss: Truss, dt: float, gravity: np.ndarray, *, solver_state: dict) -> np.ndarray:
      """Return new positions for current step given truss.points, velocities."""
  ```
  - `solver_state` allows caching matrices/preconditioners between steps.
  - Velocity update stays in `run_dynamics()` (`v = (x_new - x_old) / dt`).

- **Existing solvers to plug in**
  - Projective Dynamics dense/iterative Python versions (`python/pyTruss/projective_dynamics.py`, `python/pyTruss/projective_dynamics_iterative.py`).
  - Julia PD reference (`Julia/try_julia_plot/ProjectiveDynamics.jl`) mirrors the same structure: outer implicit loop + solver selection (`solve_linear_system!`, Gauss-Seidel, Jacobi, etc.).
  - Vertex Block Descent (`solve_vbd_numpy` and OpenCL kernel) should implement the same interface, returning relaxed positions per step.

- **System assembly reuse**
  - `Truss` exposes geometry/mass/rest-length queries; add helpers to build inertia term (`M/h^2`), gather constraint projections shared by PD and VBD.
  - Solver receives `y = x + dt*v + dt^2*a_ext`, fixed mask, bond data; avoids duplicating setup code seen today in `solve_vbd_numpy()`.

- **Next steps before coding**
  - Refactor CPU PD/VBD solvers to adopt common interface (`SolverContext` objects storing pre-factorizations).
  - Implement `Truss.run_dynamics()` orchestrating: compute `y`, call solver, update `(x, v)`, append trajectories if requested.
  - Update CLI (`run_vbd_cloth.py`) to select solver via flag (`--solver pd|vbd|gs|jacobi`), pass options through to `run_dynamics()`.

---

## Notation

* Number of vertices: $N$
* Positions: $x \in \mathbb{R}^{3N}$, with vertex $i$ at $x_i \in \mathbb{R}^3$
* Time step size: $h$
* Vertex mass: $m_i$
* External acceleration: $a^{\mathrm{ext}}$
* Force elements for vertex $i$: $F_i$
* Total potential energy: $E(x) = \sum_j E_j(x)$

---

## 3.1 Global optimization

$$
x^{t+1} = \arg\min_x G(x)
$$

$$
G(x) = \tfrac{1}{2h^2}|x-y|_M^2 + E(x), \quad y = x^t + h v^t + h^2 a^{\mathrm{ext}}
$$

Velocity:
$$
v^{t+1} = \tfrac{1}{h}(x^{t+1} - x^t)
$$

**Vertex-local energy:**
$$
G_i(x) = \tfrac{m_i}{2h^2}||x_i - y_i||_M^2 + \sum_{j \in F_i} E_j(x)
$$

---

## 3.2 Local system solver

Newton step:
$$
H_i \Delta x_i = f_i
$$

with
$$
f_i = -\tfrac{m_i}{h^2}(x_i-y_i) - \sum_{j \in F_i}\tfrac{\partial E_j(x)}{\partial x_i}
$$

$$
H_i = \tfrac{m_i}{h^2}I_{3\times3} + \sum_{j \in F_i}\tfrac{\partial^2 E_j(x)}{\partial x_i^2}
$$

* Solve $\Delta x_i = H_i^{-1} f_i$ analytically.
* Skip if $|\det H_i| \le \varepsilon$.

Solver of 3x3 matrix
```OpenCL
inline int solve3x3(__private const float* H, __private const float* f, __private float* dx, float eps) {
    float c00 = H[4]*H[8] - H[5]*H[7];
    float c01 = H[2]*H[7] - H[1]*H[8];
    float c02 = H[1]*H[5] - H[2]*H[4];
    float c10 = H[5]*H[6] - H[3]*H[8];
    float c11 = H[0]*H[8] - H[2]*H[6];
    float c12 = H[2]*H[3] - H[0]*H[5];
    float c20 = H[3]*H[7] - H[4]*H[6];
    float c21 = H[1]*H[6] - H[0]*H[7];
    float c22 = H[0]*H[4] - H[1]*H[3];

    float det = H[0]*c00 + H[1]*c10 + H[2]*c20;
    if(fabs(det) < eps) return 0;

    float invDet = 1.0f/det;

    dx[0] = (c00*f[0] + c01*f[1] + c02*f[2]) * invDet;
    dx[1] = (c10*f[0] + c11*f[1] + c12*f[2]) * invDet;
    dx[2] = (c20*f[0] + c21*f[1] + c22*f[2]) * invDet;
    return 1;
}
```

**Optional line search:**
$$
x_i \leftarrow x_i + \alpha \Delta x_i, \quad \alpha \in (0,1]
$$

**Damping:**
$$
H_i \leftarrow H_i + \tfrac{k_d}{h}\sum_{j \in F_i}\tfrac{\partial^2 E_j}{\partial x_i^2}
$$
$$
f_i \leftarrow f_i - \Big(\sum_{j \in F_i} \tfrac{k_d}{h}\tfrac{\partial^2 E_j}{\partial x_i^2}\Big)(x_i - x_i^t)
$$



---

## 3.5 Collisions

**Energy:**
$$
E_c(x) = \tfrac{1}{2}k_c d^2, \quad d = \max(0,(x_b-x_a)\cdot\hat{n})
$$

**Gradient:**
$$
\frac{\partial E_c}{\partial x_i} = k_c d \Big( \frac{\partial (x_b - x_a)}{\partial x_i}\cdot\hat{n} \Big)
$$

**Hessian (approx):**
$$
\frac{\partial^2 E_c}{\partial x_i^2} \approx k_c \frac{\partial d}{\partial x_i}\frac{\partial d}{\partial x_i}^\top + k_c d \frac{\partial^2 d}{\partial x_i^2}
$$

---

## Friction

Relative increment:
$$
\delta x_c = (x_a - x_a^t) - (x_b - x_b^t)
$$

Projection:
$$
u_c = T_c^\top \delta x_c
$$

Friction force:
$$
f_{c,i} = -\mu_c \lambda_{c,i},\frac{\partial \delta x_c}{\partial x_i}T_c f_1(|u_c|)\frac{u_c}{|u_c|}
$$

with
$$
f_1(u) =
\begin{cases}
2(u/\epsilon_v h) - (u/\epsilon_v h)^2, & 0<u<\epsilon_v h \\
1, & u \ge \epsilon_v h
\end{cases}
$$

---

## 3.7 Initialization

Choices:

1. $x = x^t$
2. $x = x^t + h v^t$
3. $x = x^t + h v^t + h^2 a^{\mathrm{ext}}$

**Adaptive:**
$$
a^t = \tfrac{v^t-v^{t-1}}{h}, \quad a^t_{\mathrm{ext}} = a^t\cdot\hat{a}^{\mathrm{ext}}, \quad \hat{a}^{\mathrm{ext}} = \tfrac{a^{\mathrm{ext}}}{|a^{\mathrm{ext}}|}
$$

$$
\alpha =
\begin{cases}
1, & a^t_{\mathrm{ext}}>|a^{\mathrm{ext}}| \\
0, & a^t_{\mathrm{ext}}<0 \\
\tfrac{a^t_{\mathrm{ext}}}{|a^{\mathrm{ext}}|}, & \text{otherwise}
\end{cases}
$$

$$
\tilde a = \alpha a^{\mathrm{ext}}, \quad x = x^t + h v^t + h^2 \tilde a
$$

---

## 3.8 Accelerated iterations (Chebyshev)

$$
x^{(n)} = \omega_n(\bar{x}^{(n)} - x^{(n-2)}) + x^{(n-2)}
$$

$$
\omega_1=1,\quad \omega_2=\tfrac{2}{2-\rho^2},\quad \omega_n=\tfrac{4}{4-\rho^2\omega_{n-1}}
$$

Skip for colliding vertices.

---

## 3.9 Parallelization

* Precolor vertices (material forces).
* Collisions: no recoloring, use buffer $x_{new}$.
* Steps per color:

  1. Compute $f_i,H_i$
  2. Solve $\Delta x_i = H_i^{-1} f_i$
  3. Store to $x_{new}$, then copy back

---

## Algorithm 1: One VBD time step

* nInput: $x^t, v^t, a_{ext}$
* Output: $x^{t+1}, v^{t+1}$

1. $y ← x^t + h v^t + h^2 a_{ext}$
2. Run DCD on $x^t$ → collisions
3. x ← initialization (adaptive)
4. For n=1..n_max:
   a. If n mod n_col=1: run CCD, update collisions
   b. For each color c:
      Parallel over vertices i in color c:
        - $f_i = -∂G_i/∂x_i$
        - $H_i = ∂^2 G_i/∂x_i^2 + (m_i/h^2)I$
        - $Δx_i = H_i^{-1} f_i$  (skip if det small)
        - $x_{new}[i] = x_i + Δx_i$
      Copy $x_{new}→x$ for color c
   c. Chebyshev acceleration (non-colliding)
5. $v^{t+1} = (x-x^t)/h$
6. Return x, v

---

## Implementation notes

* Per-vertex 3x3 inverse (analytic).
* Skip small-determinant Hessians.
* Use buffer $x_{new}$ to avoid races.
* GPU: one block per vertex, threads sum $F_i$.
* Tune $\rho$ (0.75–0.95).
* CCD every $n_{col}$ iterations.

---

## Gogle Gemini version

Of course. Here is the content from the LaTeX document converted to Markdown, with equations formatted as requested. This version is structured to be a clear and practical guide for implementation.

***

# Vertex Block Descent (VBD) Algorithm for Implementation

This document provides a detailed explanation of the Vertex Block Descent (VBD) algorithm, focusing on the core equations and logic necessary for implementation in C++ and OpenCL.

## Global Optimization

The VBD method is based on the optimization form of implicit Euler time integration. The goal is to find the next set of vertex positions $x^{t+1}$ that minimizes a global variational energy function $G(x)$.

The optimization problem is:
$$
x^{t+1} = \underset{x}{\text{argmin}} \ G(x)
$$

The variational energy $G(x)$ consists of an inertia potential and a potential energy term $E(x)$ (e.g., elastic, collision energies):
$$
G(x) = \frac{1}{2} \|x - y\|_M^2 + E(x)
$$
Where:
*   $\| \cdot \|_M$ is the mass-weighted norm.
*   $y$ is the predicted position based on the previous state, incorporating velocity and external acceleration ($a_{\text{ext}}$):
    $$
    y = x^t + h v^t + h^2 a_{\text{ext}}
    $$

VBD uses a coordinate descent approach, optimizing one vertex at a time. The local variational energy $G_i(x)$ for a single vertex $i$ is:
$$
G_i(x) = \frac{m_i}{2h^2} \|x_i - y_i\|^2 + \sum_{j \in \mathcal{F}_i} E_j(x)
$$
Where:
*   $m_i$ is the mass of vertex $i$.
*   $\mathcal{F}_i$ is the set of all force/energy elements (like tetrahedra, springs, collision constraints) acting on vertex $i$.
*   $E_j(x)$ is the potential energy of element $j$.

The position of each vertex is updated by minimizing its local energy, and these updates are performed iteratively in a Gauss-Seidel manner.
$$
x_i \leftarrow \underset{x_i}{\text{argmin}} \ G_i(x)
$$

After all iterations, the final velocity is computed:
$$
v^{t+1} = \frac{1}{h} (x^{t+1} - x^t)
$$

## Local System Solver

To minimize the local energy $G_i$ for vertex $i$, we use a single Newton's method step. This involves solving the 3x3 linear system:
$$
H_i \Delta x_i = -f_i
$$
Where:
*   $\Delta x_i$ is the change in position for vertex $i$.
*   $f_i$ is the **gradient** of the local energy $G_i$ with respect to $x_i$:
    $$
    f_i = \frac{\partial G_i(x)}{\partial x_i} = \frac{m_i}{h^2} (x_i - y_i) + \sum_{j \in \mathcal{F}_i} \frac{\partial E_j(x)}{\partial x_i}
    $$
*   $H_i$ is the **Hessian** (second derivative) of the local energy $G_i$:
    $$
    H_i = \frac{\partial^2 G_i(x)}{\partial x_i \partial x_i^T} = \frac{m_i}{h^2} I + \sum_{j \in \mathcal{F}_i} \frac{\partial^2 E_j(x)}{\partial x_i \partial x_i^T}
    $$

The position update is then $\Delta x_i = -H_i^{-1} f_i$. For a 3x3 system, the inverse $H_i^{-1}$ can be computed analytically and efficiently.

**Stability Check:** To prevent issues with singular matrices, if $|\text{det}(H_i)|$ is below a small threshold $\epsilon$, the update for vertex $i$ is skipped for that iteration.

## Collisions

Collisions are handled by adding a quadratic penalty potential energy $E_c$ for each contact constraint.
$$
E_c(x) = \frac{1}{2} k_c d^2, \quad \text{with} \quad d = \text{max}(0, (x_b - x_a) \cdot \hat{n})
$$
Where:
*   $k_c$ is the collision stiffness.
*   $d$ is the penetration depth.
*   $x_a$ and $x_b$ are the contact points on the two colliding primitives.
*   $\hat{n}$ is the contact normal.

The gradient and Hessian of this energy term are added to the corresponding local systems ($f_i$ and $H_i$) for all vertices involved in the contact. The contact normal $\hat{n}$ is assumed to be constant during differentiation.

## Friction

Friction is modeled as a force that is added to the system. To fit this into the optimization framework, the friction force and its derivative are added to the gradient and Hessian, respectively.

1.  **Relative Motion:** First, calculate the tangential relative motion at the contact point.
    $$
    \delta x_c = (x_a - x_a^t) - (x_b - x_b^t)
    $$
    Project this motion onto the tangent plane to get the 2D tangential translation $u_c = T_c^T \delta x_c$.

2.  **Friction Force:** The friction force vector for vertex $i$ is:
    $$
    F^{\text{fric}}_{c,i} = -\mu_c \lambda_{c,i} \left( \frac{\partial \delta x_c}{\partial x_i} \right)^T T_c f_1(\|u_c\|) \frac{u_c}{\|u_c\|}
    $$
    Where:
    *   $\mu_c$ is the friction coefficient.
    *   $\lambda_{c,i}$ is the magnitude of the normal force on vertex $i$, derived from the collision energy gradient.
    *   $f_1(u)$ is a function that smoothly transitions from static to dynamic friction:
        $$
        f_1(u) =
        \begin{cases}
            \frac{u}{h e_v} - \frac{u^2}{2 (h e_v)^2}, & \text{if } u \in (0, he_v) \\
            1, & \text{if } u \ge he_v
        \end{cases}
        $$

3.  **Gradient and Hessian Update:**
    *   Add the friction force contribution to the gradient: $f_i \leftarrow f_i - F^{\text{fric}}_{c,i}$.
    *   Add the friction Hessian contribution: $H_i \leftarrow H_i - \frac{\partial F^{\text{fric}}_{c,i}}{\partial x_i}$. The paper provides the following approximation for the friction Hessian:
        $$
        \frac{\partial F^{\text{fric}}_{c,i}}{\partial x_i} \approx -\mu_c \lambda_{c,i} \left( \frac{\partial \delta x_c}{\partial x_i} \right)^T T_c \frac{f_1(\|u_c\|)}{\|u_c\|} T_c^T \frac{\partial \delta x_c}{\partial x_i}
        $$

## Initialization

A good initial guess for $x$ at the start of the iterative process is crucial for fast convergence. The paper proposes an adaptive scheme:
$$
x_{\text{initial}} = x^t + h v^t + h^2 \tilde{a}
$$
Where $\tilde{a}$ is a carefully scaled version of the external acceleration $a_{\text{ext}}$.
$$
\tilde{a} = \alpha a_{\text{ext}}
$$
The scaling factor $\alpha$ is determined by comparing the previous frame's acceleration $a^t = (v^t - v^{t-1})/h$ to the external acceleration:
*   Let $\tilde{a}_{\text{ext}} = a^t \cdot \hat{a}_{\text{ext}}$ be the component of the previous acceleration along the external acceleration direction.
*   Then $\alpha$ is calculated as:
    $$
    \alpha =
    \begin{cases}
        1, & \text{if } \tilde{a}_{\text{ext}} > \|a_{\text{ext}}\| \\
        0, & \text{if } \tilde{a}_{\text{ext}} < 0 \\
        \tilde{a}_{\text{ext}} / \|a_{\text{ext}}\|, & \text{otherwise}
    \end{cases}
    $$
This scheme helps objects in free-fall to start closer to their final position, while preventing objects at rest from being initialized in a penetrating state.

## Accelerated Iterations

Convergence can be improved using the Chebyshev semi-iterative approach. After each full Gauss-Seidel pass (over all colors), the positions are updated globally.

Let $\tilde{x}^{(n)}$ be the vertex positions after the Gauss-Seidel pass in iteration $n$. The accelerated position $x^{(n)}$ is then computed as:
$$
x^{(n)} = \omega_n ( \tilde{x}^{(n)} - x^{(n-2)} ) + x^{(n-2)}
$$
The acceleration ratio $\omega_n$ changes at each iteration:
$$
\omega_n = \frac{4}{4 - \rho^2 \omega_{n-1}}
$$
With initial conditions $\omega_1 = 1$, $\omega_2 = \frac{2}{2-\rho^2}$, and $\rho \in (0, 1)$ being the estimated spectral radius (e.g., a constant like 0.95).

## Parallelization

*   **Static Forces:** Vertices are pre-colored based on the mesh connectivity (e.g., tetrahedra). All vertices of the same color can be processed in parallel because they don't share any static force elements. The updates are performed serially over the color groups.

*   **Dynamic Forces (Collisions):** To handle race conditions from dynamically created collision constraints, a temporary position buffer is used.
    *   During a color pass, updates for vertices are written to an auxiliary buffer `x_new`.
    *   After all vertices of a color are processed, the contents of `x_new` are copied back to the main position buffer `x`.
    This ensures that vertices within the same color group use position data from the *start* of the color pass, effectively mixing Jacobi-style updates for dynamic forces with Gauss-Seidel updates for static forces.

*   **GPU workgroup scheduling (TODO):** Process each pre-colored set sequentially on the device by issuing one kernel launch per color. Within a launch, split the color into chunks that map to workgroups of size 32. Each chunk must include all vertices assigned to the workgroup (up to 32 indices) and the complete list of their neighbors (up to 128 indices) so that no `neighs` entry is omitted. The Python host code must build these packed groups before dispatch, ensuring every neighbor reference is resident in the chunk buffers.

***

# Algorithm 1: VBD simulation for one time step

```pseudocode
function VBD_TimeStep(x_t, v_t, a_ext):
    // --- Initialization ---
    y <- x_t + h*v_t + h^2*a_ext
    // Perform initial Discrete Collision Detection (DCD) to find penetrating vertices
    InitialDCD(x_t)
    // Set the initial guess for the iterative solver
    x <- AdaptiveInitialization(x_t, v_t, v_{t-1}, a_ext)

    // --- Main Iteration Loop ---
    for n from 1 to n_max:
        // Perform Continuous Collision Detection (CCD) periodically
        if n % n_col == 1:
            CCD(x)

        // Iterate through each color group sequentially
        for each color c:
            // --- Block-level parallelism (e.g., one GPU block per vertex) ---
            parallel for each vertex i in color c:
                // --- Thread-level parallelism (e.g., threads in a block sum forces) ---
                // 1. Compute sums of gradients (f_i) and Hessians (H_i)
                f_i <- (m_i / h^2) * (x_i - y_i)
                H_i <- (m_i / h^2) * I

                // Aggregate from all force elements (elastic, collision, friction)
                parallel for each element j connected to vertex i:
                    f_i <- f_i + gradient(E_j, x_i)
                    H_i <- H_i + hessian(E_j, x_i)

                // 2. Solve the 3x3 local system
                delta_x_i <- -inverse(H_i) * f_i

                // 3. Store result in temporary buffer to avoid race conditions
                x_new[i] <- x_i + delta_x_i

            // Copy updated positions back to the main buffer after the color pass
            parallel for each vertex i in color c:
                x[i] <- x_new[i]

        // --- Optional: Acceleration Step (after all colors are processed) ---
        if use_acceleration:
            parallel for each vertex i:
                Update x[i] using the Chebyshev formula

    // --- Finalization ---
    v_t_plus_1 <- (x - x_t) / h
    return x, v_t_plus_1
```
