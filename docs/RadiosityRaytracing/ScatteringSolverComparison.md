# Scattering Solver Comparison: `Scatterer.h` vs `Scatterer2.h`

This document compares two C++ headers, `Scatterer.h` and `Scatterer2.h`, which implement different approaches to solving light or particle scattering problems.

## 1. `Scatterer.h` (Radiosity-like / Matrix-based)

*   **Core Idea:** This class, `Scattering`, is designed to solve a global linear system that describes the exchange of "flux" (e.g., light energy) between discrete surface elements. This approach is highly analogous to the classic Radiosity algorithm used in computer graphics.
*   **Key Components & Functionality:**
    *   **`SurfElement`:** Represents a small patch on a surface, defined by its position, normal, area, and surface ID.
    *   **`makeCouplingMatrix()`:** This is the central method. It explicitly constructs a dense, symmetric `M` matrix. Each element `M[i][j]` represents the "coupling" (or form factor) between `SurfElement i` and `SurfElement j`. This coupling considers:
        *   **Geometric Factor:** Calculated based on the relative orientation and distance between elements (`eli.geomCoupling(elj)`).
        *   **Occlusion:** It uses `getOcclusion()` (from the `TriangleRayTracer` base class) to determine if any obstacles block the line of sight between elements, ensuring realistic light transport.
    *   **`dotFunc()`:** This method overrides a virtual function from the `LinSolver` base class. It implements the matrix-vector multiplication `Ax` using the pre-computed `M` matrix. This indicates that the class is intended to be used with an iterative linear solver (such as Conjugate Gradient, typically supported by `LinSolver`) to find the unknown flux values.
    *   **`step_Direct()`:** This method performs a single, direct update step: `vals[i] = M[i*n+j] * (vals[j] + sources[j])`. In the context of a linear solver, this would typically represent one iteration of a Jacobi or Gauss-Seidel method if used iteratively, or a direct propagation step if the system is simple enough.
*   **Analogy to Radiosity:** The explicit construction of a global coupling matrix, which incorporates geometric factors and occlusion, is a defining characteristic of the Radiosity method. The objective is to determine the equilibrium distribution of flux by solving a system of linear equations.

## 2. `Scatterer2.h` (Direct Iterative Simulation / Deterministic Monte Carlo-like)

*   **Core Idea:** This class, `Scattering2`, simulates the transport of "flux" through a network of discrete "channels" that connect scattering elements. Unlike the matrix-based approach, it directly propagates flux through these channels in an iterative fashion.
*   **Key Components & Functionality:**
    *   **`HalfChannel`:** Represents one end of a connection, storing `fluxIn` and `fluxOut` values.
    *   **`Channel`:** Represents a bidirectional connection between two scattering elements, defined by a direction (`dir`) and two `HalfChannel` ends.
    *   **`ScatterElem`:** Similar to `SurfElement` but also includes indices to its associated `Channel`s, indicating its connectivity within the network.
    *   **`makeChannles()`:** This method establishes the `Channel` connections between `ScatterElem`s, also taking `occlusion` into account.
    *   **`scatterAmp()`:** Calculates the scattering amplitude between two directions for a `ScatterElem`, similar to the scattering kernel in `Scatterer.h`.
    *   **`scatterBiDir()`:** This is a crucial method. It takes two `Channel`s and simulates a bidirectional scattering event between them, transferring `fluxOut` from one to `fluxIn` of the other based on the calculated `scatterAmp`.
    *   **`dotFunc()`:** This method is explicitly empty (`override {}`). This is a strong indicator that `Scattering2` *does not* utilize the `LinSolver`'s matrix-solving capabilities. It does not form or solve a global linear system.
    *   **`step_Direct()`:** This is the primary simulation loop. It operates in two distinct phases:
        *   **Scattering Phase:** Iterates through all `ScatterElem`s and their associated `Channel`s, calling `scatterBiDir()` to simulate local scattering events.
        *   **Transfer Phase:** Iterates through all `Channel`s, updating `fluxOut` based on `fluxIn` for both ends. This simulates the propagation of flux along the channels.
*   **Analogy to Monte Carlo / Discrete Ordinates:** This approach is more akin to a deterministic Monte Carlo simulation or a discrete ordinates method. Instead of solving a global system, it directly simulates the "movement" and "scattering" of flux packets (or discrete flux values) between predefined angular and spatial bins (the channels). It's an explicit, iterative propagation.

## Key Differences and Relationship to Monte Carlo/Radiosity

| Feature             | `Scatterer.h` (Older)                               | `Scatterer2.h` (Newer)                                     |
| :------------------ | :-------------------------------------------------- | :--------------------------------------------------------- |
| **Primary Approach** | **Radiosity-like (Matrix-based)**                   | **Direct Iterative Simulation (Deterministic Monte Carlo-like)** |
| **Global System**   | Builds and intends to solve a dense `M` matrix.     | Does NOT build a global matrix.                            |
| **`LinSolver` Use** | Fully utilizes `LinSolver` for matrix operations.   | Inherits `LinSolver` but `dotFunc()` is empty; does not use matrix solving. |
| **Flux Calculation**| Solves for global equilibrium flux distribution.    | Explicitly propagates flux through discrete channels.      |
| **Simulation Step** | `step_Direct()` is one iteration of a linear solver. | `step_Direct()` is a full simulation cycle of scattering and transfer events. |
| **Complexity**      | Higher memory usage for `M` matrix (N^2). Computationally intensive for large N (matrix inversion/iterative solution). | Lower memory usage (stores channels, not full matrix). Computationally scales more linearly with N. |
| **Nature**          | Global, implicit, equilibrium-seeking.              | Local, explicit, step-by-step propagation.                 |

## Usage in `test_Scatterer.cpp`

The `test_Scatterer.cpp` file (and its identical counterpart `test_Scatterer-new.cpp`) exclusively uses `Scatterer2.h`. The lines related to `makeCouplingMatrix()` and `prepare()` (which would be necessary for the matrix-based approach of `Scatterer.h`) are commented out. This confirms that the direct iterative simulation approach of `Scatterer2.h` is the one actively implemented and tested in this specific application.

In essence, `Scatterer.h` represents a more traditional, global approach to solving transport problems via linear algebra, while `Scatterer2.h` offers a more direct, event-driven simulation that can be more flexible and scalable for certain types of problems.

### 4. Comparison of Element Types (`SurfElement` vs. `ScatterElem*`)

The simulation uses different C++ structs to represent elements, each with a specific role. Understanding them is key to understanding the different solvers.

#### `SurfElement` (`TriangleRayTracer.h`)
*   **Purpose:** A purely **geometric** entity. It represents a small, flat patch (a panel) on a larger surface.
*   **Key Properties:**
    *   `pos`: Position in 3D space.
    *   `normal`: The surface normal vector, defining its orientation.
    *   `area`: The surface area of the patch.
    *   `isurf`: An ID for the original triangle it belongs to, used to prevent self-occlusion during ray tracing.
*   **Role:** `SurfElement` is the fundamental building block for discretizing geometry. It is used by the `TriangleRayTracer` to define surfaces that can cause occlusion and by the `Scattering` (Radiosity) solver as the basis for exchanging flux. It is essentially a passive surface.

#### `ScatterElem2` (`Scatterer2.h`)
*   **Purpose:** A **physical** entity representing an active scattering center or node in the simulation network.
*   **Key Properties:**
    *   `pos`, `isurf`: Shares these geometric properties with `SurfElement`.
    *   `rot`, `thicks`, `areas`, `beta`: These define the anisotropic scattering physics of the element. It does not have a single `normal` or `area`; its interaction depends on the direction of incoming and outgoing flux.
    *   `ichan`, `nchan`, etc.: These are connectivity indices, linking the element to the `channels` that transport flux.
*   **Role:** `ScatterElem2` is the core component of the `Scattering2` direct simulation. It's not just a passive surface; it's an active node that receives flux from channels, re-distributes it according to physical scattering laws, and sends it back out into other channels.

#### `ScatterElement` (`Scatterer.h`)
*   This is another, more complex, physical entity from the older `Scatterer.h`. It also represents an active scattering point but focuses on discretizing the outgoing solid angle into a grid of `rays` to manage scattered flux distribution.

### Should `ScatterElem2` inherit from `SurfElement`?

While they share properties like `pos` and `isurf`, **inheritance is not recommended**. The "is-a" relationship does not hold well. A `ScatterElem2` is not a special type of `SurfElement`; they represent different concepts.

*   **Conceptual Mismatch:** `SurfElement` is a 2D surface patch. `ScatterElem2` is a point-like object with volumetric, anisotropic scattering properties. Forcing inheritance would mean `ScatterElem2` inherits a `normal` and `area` that don't accurately describe its physics.
*   **Better Approach (Composition/Factory):** A better design is to keep them separate. The geometry discretization process (e.g., `trinagleToElements2`) creates a list of `SurfElement`s. Then, a factory or setup function can iterate through these `SurfElement`s and create a corresponding list of `ScatterElem2`s, copying over the `pos` and `isurf` and initializing the physical properties (`rot`, `thicks`, etc.) as needed.

This separation of concerns (geometry vs. physics) leads to a cleaner and more flexible design. The current implementation in `test_Scatterer.cpp` has a logical bug where it creates `SurfElement`s but the solver expects `ScatterElem2`s, highlighting the need for this explicit creation step.
