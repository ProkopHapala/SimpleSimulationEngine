
## Efficient Representation and Interpolation of Spherical Functions using Octahedral Mapping

### Abstract

Representing and efficiently interpolating functions defined on the surface of a sphere is a fundamental problem in computer graphics, physics, and data visualization. Common methods like latitude-longitude maps suffer from polar singularities and non-uniform sampling, while Spherical Harmonics can be computationally expensive for high-frequency data. This document provides a complete guide to an alternative, elegant solution: **Octahedral Mapping**. We detail the mathematical projection from 3D space to a 2D texture, the construction of a continuous triangular grid for data storage, and a robust algorithm for performing fast, accurate barycentric interpolation. This method is particularly well-suited for applications requiring fast, random-access queries of spherical data, such as environment maps, light probes, and angular distribution functions in global illumination.

### 1. The Challenge: Data on a Sphere

Many real-world phenomena are best described as functions on a sphere. A function $f(\mathbf{d})$ takes a 3D direction vector $\mathbf{d}$ (where $\|\mathbf{d}\|=1$) and returns a value, such as color, intensity, or temperature. Examples include:

*   **Environment Maps:** The color of light arriving from every direction.
*   **Radiance Fields:** The outgoing light energy from a point in space for every direction (essential for radiosity and global illumination).
*   **Planetary Data:** Storing elevation or climate data on a globe.

The challenge is to store this continuous spherical function in a discrete, finite data structure that allows for efficient and accurate reconstruction (interpolation) of the function's value at any given direction.

#### 1.1. Common but Flawed Approaches

*   **Latitude-Longitude (Equirectangular) Maps:** This method maps spherical coordinates $(\theta, \phi)$ directly to a 2D grid. While simple, it suffers from severe distortion and oversampling near the poles, leading to wasted memory and interpolation artifacts.
*   **Cube Maps:** Projects the sphere onto the six faces of a cube. This is a significant improvement, but it introduces discontinuities at the cube edges and corners, which require special handling during sampling. The pixel area also varies by up to 40% from the center to the corner of a face.
*   **Spherical Harmonics (SH):** Represents the function as a sum of basis functions. SH is perfect for low-frequency, diffuse functions but becomes computationally expensive and prone to ringing artifacts when representing high-frequency (detailed) data. Evaluation requires a costly summation over all coefficients.

The octahedral map provides a superior balance of sampling uniformity, computational efficiency, and implementation simplicity.

### 2. The Octahedral Projection: From Sphere to Square

The core idea of octahedral mapping is to inscribe a regular octahedron within the sphere and "unfold" its eight triangular faces into a single 2D square texture.

* [Octahedron Environment Map,Thomas Engelhardt, Carsten Dachsbacher, 2008](https://www.readkong.com/page/fullscreen/octahedron-environment-maps-6054207)
* [Survey of Efficient Representations for Independent Unit Vectors, Zina H. Cigolle et al.](https://jcgt.org/published/0003/02/01/)
* [What is Octahedral Compression of Vertex Arrays?](https://stackoverflow.com/questions/74743644/what-is-octahedral-compression-of-vertex-arrays)
   * [Normals Compression - Octahedron](https://www.shadertoy.com/view/Mtfyzl)
   * [Normals Compression Comparison](https://www.shadertoy.com/view/4llcRl)
   * [RTXGI Octahedron Mapping ](https://www.shadertoy.com/view/Ddy3WG)


#### 2.1. Geometric Intuition

An octahedron has 6 vertices and 8 triangular faces. We can orient it such that its vertices align with the Cartesian axes: `(±1,0,0)`, `(0,±1,0)`, `(0,0,±1)`. The four faces with a positive Z component form the "upper hemisphere," and the four with a negative Z component form the "lower hemisphere." The mapping cleverly arranges these faces as follows:
*   The **upper hemisphere** faces form a central diamond in the UV square.
*   The **lower hemisphere** faces are "folded" outwards to fill the four corners of the UV square.

This arrangement ensures that vertices and edges that were adjacent on the 3D octahedron remain adjacent in the 2D map, which is the key to seamless interpolation.



#### 2.2. The Mathematical Mapping

The projection from a 3D direction vector $\mathbf{p} = (x, y, z)$ to a 2D UV coordinate $\mathbf{uv} = (u, v) \in [-1, 1]^2$ is a two-step process.

**Step 1: Project Sphere to Octahedron Surface**

First, normalize the input vector $\mathbf{p}$ so it lies on the unit sphere. The projection from the sphere to the surface of the octahedron is achieved by dividing the vector by its **L1 norm** (the sum of the absolute values of its components).

$$ \mathbf{p}' = \frac{\mathbf{p}}{\|\mathbf{p}\|_1} = \frac{(x, y, z)}{|x| + |y| + |z|} $$

This elegant operation maps any point on the sphere to a corresponding point on the octahedron.

**Step 2: Unfold the Octahedron**

The $(x, y)$ components of the projected point $\mathbf{p}'$ give us our initial UV coordinates. For the upper hemisphere ($z \ge 0$), this is the final mapping. For the lower hemisphere ($z < 0$), a "folding" transformation is required to map the point correctly into the corners of the UV square.

Let $(u', v') = (p'_x, p'_y)$.

$$ (u, v) = \begin{cases} (u', v') & \text{if } z \ge 0 \\ \left( (1 - |v'|) \cdot \text{sign}(u'), (1 - |u'|) \cdot \text{sign}(v') \right) & \text{if } z < 0 \end{cases} $$

This mapping is continuous, efficient (requiring only a few absolute values, additions, and multiplications), and forms the foundation of our system.

### 3. Discretization and Continuous Triangulation

To store data, we impose a uniform grid on the $[-1, 1]^2$ UV square. A grid of $(N+1) \times (N+1)$ vertices creates a grid of $N \times N$ quads (pixels). While bilinear interpolation within these quads is possible, a more robust method for ensuring continuity is to subdivide each quad into two triangles and use **barycentric interpolation**.

#### 3.1. The Diagonal Splitting Problem

A naive choice, such as splitting every quad with a diagonal from bottom-left to top-right, will create T-junctions and interpolation artifacts along the boundaries of the original octahedron faces (i.e., the U and V axes).

**The Solution:** The orientation of the diagonal must alternate in a specific pattern that ensures the triangulation is continuous across all edges. The correct pattern depends on the UV quadrant:

*   For quads in the **top-right** and **bottom-left** quadrants of the UV map, the diagonal runs from top-left to bottom-right (`\`).
*   For quads in the **top-left** and **bottom-right** quadrants, the diagonal runs from bottom-left to top-right (`/`).

This "starburst" pattern guarantees that any two quads sharing an edge will have a consistent triangulation, eliminating seams.



### 4. The Interpolation Algorithm: A Developer's Guide

This section provides a step-by-step algorithm to find the value of the spherical function for an arbitrary 3D direction vector $\mathbf{p}$.

**Input:**
*   A 3D direction vector `point_3d`.
*   A 2D array `grid_data` of size `(N+1, N+1)` storing the function values at each grid node.

**Output:**
*   The interpolated function value.

---

**Algorithm Steps:**

**1. Project to UV Space:**
   Apply the octahedral mapping from Section 2.2 to `point_3d` to obtain its corresponding `uv_point = (u, v)`.

**2. Localize the Containing Quad:**
   Convert the `uv_point` from the `[-1, 1]` range to grid index space `[0, N]`.
   
   ```
   u_grid = (u + 1.0) / 2.0 * N
   v_grid = (v + 1.0) / 2.0 * N
   
   col_idx = floor(u_grid)
   row_idx = floor(v_grid) 
   ```
   
   The point lies within the quad defined by grid nodes `(row_idx, col_idx)` and `(row_idx+1, col_idx+1)`.

**3. Identify the Correct Sub-Triangle:**
   This is the most critical step and must be consistent with the diagonal pattern.

   a.  Determine the quad's diagonal orientation based on its position relative to the grid's center (`N/2`).
   b.  The interpolation triangle's first two vertices are **always the two endpoints of that diagonal**.
   c.  Calculate the point's local coordinates `(u_local, v_local)` within the quad, ranging from `[0, 1]`.
      ```
      u_local = u_grid - col_idx
      v_local = v_grid - row_idx
      ```
   d.  Compare the point's position to the diagonal line equation to select the third, "off-diagonal" vertex.

      *   **If the diagonal is `/`:** The line is `v_local = u_local`. The triangle vertices are `(v_bl, v_tr, v_tl)` if `v_local > u_local`, or `(v_bl, v_tr, v_br)` otherwise.
      *   **If the diagonal is `\`:** The line is `v_local = 1 - u_local`. The triangle vertices are `(v_tl, v_br, v_tr)` if `v_local > 1 - u_local`, or `(v_tl, v_br, v_bl)` otherwise.

**4. Calculate Barycentric Coordinates:**
   With the three vertices of the triangle (`V0`, `V1`, `V2`) identified, solve for the barycentric weights `(w0, w1, w2)` such that:
   
   `uv_point = w0*V0 + w1*V1 + w2*V2`
   
   This can be solved as a small 2x2 linear system. The weights represent the "influence" of each vertex on the point.

**5. Perform Final Interpolation:**
   Retrieve the data values stored at the three grid node indices (`idx0`, `idx1`, `idx2`) corresponding to the triangle's vertices. The final interpolated value is the weighted average:

   `Value = w0 * grid_data[idx0] + w1 * grid_data[idx1] + w2 * grid_data[idx2]`

---

### 5. Conclusion and Applications

The octahedral mapping method provides a robust, efficient, and mathematically elegant framework for handling spherical data.

**Advantages:**
*   **Quasi-Uniform Sampling:** Offers a more uniform distribution of samples across the sphere than latitude-longitude maps.
*   **No Singularities:** Completely avoids the pole problem.
*   **Computational Efficiency:** The projection and interpolation algorithms use only simple arithmetic operations, making them extremely fast and suitable for real-time applications on GPUs.
*   **Continuity:** The specified triangulation scheme guarantees seamless C0 continuity across the entire spherical domain.

This technique is a powerful tool for any developer working on problems in **global illumination**, **radiosity**, **radiation scattering**, or any domain that requires the storage and rapid querying of angular distribution functions. By providing a direct map between 3D directions and a simple 2D array, it bridges the gap between the geometry of the sphere and the linear memory of a computer.