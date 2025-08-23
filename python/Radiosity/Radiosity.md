# Radiosity Heat Transfer

I want to create [Radiosity](https://en.wikipedia.org/wiki/Radiosity_(computer_graphics)) solver for radiative heat transfer. The input should be set of elements represented by point in space (but for debuging I would have them in 2D so that it is easier to visualize), each point has oriented surface with normal n and surface aread A and we assume [Lambert cosine law](https://en.wikipedia.org/wiki/Lambert%27s_cosine_law) for diffuse reflection reflection but also for thermal radiation. There should be some "albedo" or "absorbace/emmisivity" coeeficient \epsilon which control the coupling between heat energy in the elements and thermal radiation. This coupling should be given by well know [Stefan–Boltzmann law](https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law)

$$
P = A \epsilon \sigma T^{4}
$$

In the end we want to expres it in standard matrix form used for radiosity calculations

$$
B_i = E_i + \rho_i \sum_{j=1}^n F_{ij} B_j
$$

$$
B = (I - \rho F)^{-1} E 
$$

$$
A_i B_i = A_i E_i + \rho_i \sum_{j=1}^n A_j B_j F_{ji}
$$

However there are some chaveates:
1) To express it as linear system we need to store $T^4$ rather then T at the elements
2) the surfaces may be two-sided, where from oposite site they have different emmisivity \epsilon

Note that all surface cna both absorb and emit light. Most light is lost in outer space never hitting other surface. Most of surfaces are passive (i.e. they are in equlibirum where they emit as much energy as they absorb), but some surfaces may be out of equlibirum in stacionary state where certain constant  inflow or outplow of thermal energy comes/leave it by non-radiative transfer. This should be some right-hand side of the equation.

---

Now please analyse the problem, write it in systematic manner (both as a tutorial for student and design document for software developer).

Then implement it using python and numpy. Use efficient numpy operations (array indexing, masking). Try to avoid slow explicity loop interations and branching in python as they are much slower than numpy array operations. Finally visualize the results.

The taks comprise of building of the matrix. The input should be list of lines in 2D plane (just few for the moment) which define our surfaces. For example each surface can be defined like this ((x1,y1),(x2,y2), nelements, epsilon, heatingPower) and it is cut to nelements.
You compute matrix of optical couplings between them
Then you solve that matrix using numpy solve,
then you visualize the temperature (T) of each element by scaterpolot colors


---

This document provides a comprehensive guide to understanding and implementing a 2D radiosity solver for radiative heat transfer. It's structured to serve as both a tutorial for students learning the concepts and a design document for software developers implementing the solution.

### Part 1: Analysis and Theoretical Formulation

#### **Tutorial for Students**

##### **1. Introduction to Radiative Heat Transfer and Radiosity**

Imagine a collection of surfaces, some hot and some cold. The hot surfaces emit thermal energy in the form of radiation, and all surfaces can absorb and reflect radiation that falls on them. Radiative heat transfer is the process by which this energy moves between surfaces.

The **radiosity method** is a powerful technique to calculate the steady-state temperature distribution in such a system. It's widely used in computer graphics to create realistic lighting and in engineering to solve thermal problems.

The core concept is **Radiosity (B)**, which is the total energy leaving a surface per unit area. This includes two components:
*   **Emitted Energy (E):** Energy the surface emits on its own due to its temperature. This is like a light bulb glowing.
*   **Reflected Energy:** Energy that hits the surface from other sources and bounces off.

So, for any patch of a surface, we can write:
`Radiosity = Emitted Energy + Reflected Energy`

##### **2. The Key Ingredients**

*   **Emissivity (ε) and Absorptivity (α):** Emissivity is a value between 0 and 1 that describes how efficiently a surface radiates energy compared to a perfect "black body". For grey bodies, which we assume here, the ability to absorb energy (absorptivity α) is equal to its emissivity (ε).
*   **Reflectivity (ρ):** This is the fraction of incoming radiation that a surface reflects. For an opaque surface that doesn't transmit energy, all energy is either absorbed or reflected. Therefore, `α + ρ = 1`, which means `ρ = 1 - ε`.
*   **View Factor (Fij):** This is a crucial geometric factor. The view factor `F_ij` is the fraction of radiation leaving surface `i` that directly strikes surface `j`. It depends only on the geometry—how the two surfaces are positioned and oriented relative to each other.
*   **Stefan-Boltzmann Law:** This law connects temperature to emitted energy. The power emitted by a surface is given by `P = A * ε * σ * T^4`, where `A` is the surface area, `ε` is the emissivity, `σ` is the Stefan-Boltzmann constant, and `T` is the absolute temperature. The emitted energy per unit area, `E`, is therefore `ε * σ * T^4`.

##### **3. Assembling the Radiosity Equation**

Let's write down the energy balance for a single surface `i`:

`B_i = E_i + ρ_i * H_i`

Where:
*   `B_i` is the radiosity of surface `i`.
*   `E_i` is the self-emitted energy of surface `i`.
*   `ρ_i` is the reflectivity of surface `i`.
*   `H_i` is the total radiation arriving at surface `i` per unit area (also called irradiation).

The incoming radiation `H_i` is the sum of radiosities from all other surfaces (`B_j`), adjusted by the view factor from `i` to `j`. Using a fundamental property called the reciprocity relation (`A_i * F_ij = A_j * F_ji`), we can express the total incident energy on `i` as the sum of energy leaving all other surfaces `j`, weighted by the view factor `F_ji`.

`H_i = Σ_j (F_ji * B_j)`

Substituting this back gives the full radiosity equation for surface `i`:

`B_i = E_i + ρ_i * Σ_j (F_ji * B_j)`

Since this equation exists for every surface in the system, we have a set of linear equations that can be solved simultaneously.

##### **4. Incorporating External Heating**

Surfaces don't just gain and lose energy by radiation. They might be heated or cooled by other means (e.g., conduction, convection, or an internal power source). We can define a `heatingPower_i` term, which is the net non-radiative energy added to the surface.

In a steady state, the energy for each surface must balance:
`Energy In = Energy Out`

`heatingPower_i + Energy Absorbed = Energy Emitted`
`heatingPower_i + (A_i * ε_i * H_i) = (A_i * E_i)`

By combining this energy balance with the radiosity definition, we can derive a powerful and elegant matrix equation that can be solved for the radiosities `B`.

#### **Design Document for Developers**

##### **1. System Overview**

The objective is to create a Python-based 2D radiosity solver. The system will take a set of 2D line segments as input, discretize them into smaller elements, compute the geometric view factors between them, solve the radiosity matrix equation, and finally calculate and visualize the temperature of each element.

##### **2. Data Structures**

*   **Input Surfaces:** A list of tuples, where each tuple defines a surface:
    `((x1, y1), (x2, y2), n_elements, epsilon, heating_power)`
*   **Element Array:** A single NumPy array will store the state of all discretized elements. This is crucial for efficient vectorized computation. The columns of the array will be:
    *   `pos`: (N, 2) array for the center coordinates `(x, y)` of each element.
    *   `p1`, `p2`: (N, 2) arrays for the start and end points of each element.
    *   `length`: (N,) array for the area (length in 2D) of each element.
    *   `normal`: (N, 2) array for the outward-facing normal vector.
    *   `epsilon`: (N,) array for emissivity.
    *   `rho`: (N,) array for reflectivity (`1 - epsilon`).
    *   `heating_power`: (N,) array for the non-radiative power injected into the element.

##### **3. Core Algorithm**

1.  **Discretization:**
    *   A function `discretize_surfaces` will take the input list and populate the master element array.
    *   For each input line, it will generate `n_elements` equally sized smaller segments.
    *   It will calculate and store the properties for each element (center, length, normal, etc.).
    *   **Note on Normals / Two-sided surfaces:** The "front" face normal is defined by a 90-degree counter-clockwise rotation of the vector from `(x1, y1)` to `(x2, y2)`. For two-sided surfaces we do NOT duplicate elements. Each geometric element stores two emissivities `eps_front`, `eps_back` and a single heat source `heat_in`. The system unknowns remain one per element (N unknowns), not `2N`.

##### **Two-sided surfaces with single-DOF (N elements, not 2N)**

We represent each discretized segment once (single geometric element), but allow different optical properties on the two sides:

*   `eps_front`, `eps_back` per element; reflectivities `rho_front = 1 - eps_front`, `rho_back = 1 - eps_back`.
*   Single non-radiative input `heat_in` per element; no split into front/back heats.
*   The geometry coupling from emitter `j` to receiver `i` is split into four directional combinations using signed cosines with the inter-center direction `r̂_ij`:
    - `K_ff = max(0, n_i·r̂_ij) max(0, n_j·r̂_ji) / (π r^2) * A_j`
    - `K_fb = max(0, n_i·r̂_ij) max(0, -n_j·r̂_ji) / (π r^2) * A_j`
    - `K_bf = max(0, -n_i·r̂_ij) max(0, n_j·r̂_ji) / (π r^2) * A_j`
    - `K_bb = max(0, -n_i·r̂_ij) max(0, -n_j·r̂_ji) / (π r^2) * A_j`

Receiver-side reflectivity is applied per receiving side, yielding the effective transport operator

`W = diag(rho_front) * (K_ff + K_fb) + diag(rho_back) * (K_bf + K_bb)`

The linear system is then

`(I - W) B = heat_in / area`

which has size `N × N` with one unknown radiosity `B_i` per element.

Temperature is recovered by a thin-sheet balance with a single heat source and two-side emissivities

`(eps_front + eps_back) T_i^4 = (heat_in/area)_i + eps_front * H_front,i + eps_back * H_back,i`

where `H_front = (K_ff + K_fb) B` and `H_back = (K_bf + K_bb) B`.

2.  **View Factor Calculation:**
    *   A function `compute_view_factors(elements)` returns an `(N, N)` geometry-only matrix `F`.
    *   Implementation note: we use a midpoint inverse-square kernel with Lambertian cosines and emitter area weighting
        `F_ij ≈ (|n_i·r̂_ij| |n_j·r̂_ji| / (π r^2)) A_j`, with zero diagonal and clipping to [0,1].
    *   **Vectorization Strategy:**
        *   Pre-calculate two large distance matrices: one for all `p1` to all `p1`/`p2`, and one for all `p2` to all `p1`/`p2`. This avoids repeated distance calculations in a loop.
        *   A directional visibility split is applied later (in the solver and temperature steps) via the four `K_*` combinations; no hard zeroing mask is used in the base geometry matrix.

3.  **Solving the System:**
    *   The governing matrix equation is `(I - W) B = P`, where `W` is the effective transport operator combining geometry and receiver-side reflectivity: `W = diag(rho_front) (K_ff + K_fb) + diag(rho_back) (K_bf + K_bb)`.
    *   `B` is the unknown vector (size N). The source vector is `P_i = heat_in_i / area_i` (single heat per element).
    *   The implementation uses `numpy.linalg.solve` to obtain `B`.

4.  **Temperature Calculation:**
    *   A function `calculate_temperatures(B, F, elements)` will determine the final temperatures.
    *   First, calculate the incident radiation vector: `H = F^T * B`.
    *   Then, calculate the `T^4` value for each element. The Stefan-Boltzmann constant `σ` is assumed to be 1 for simplicity, so the result is proportional to the true `T^4`.
        `T4 = (P / epsilon) + H`
    *   Handle potential division by zero if an element has zero emissivity.

5.  **Visualization:**
    *   A function `visualize_results` will use Matplotlib.
    *   It will draw the original line segments for context.
    *   It will create a scatter plot of the element centers, with the color of each point mapped to its calculated `T^4` value. A "hot" colormap will be used.

### Part 2: Python Implementation and Visualization

Here is the complete, documented Python code that implements the radiosity solver described above.

