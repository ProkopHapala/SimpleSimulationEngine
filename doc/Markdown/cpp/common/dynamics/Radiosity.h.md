# Radiosity.h

This file defines the `Radiosity` class, which implements the radiosity algorithm for computing the light distribution in a scene. Radiosity is a global illumination algorithm that calculates the light arriving at each surface in the scene, considering the light emitted by other surfaces. The class inherits from `TriangleRayTracer` to handle ray tracing for occlusion calculations and from `LinSolver` to solve the linear system of equations that arise in the radiosity algorithm.

## Includes

- `Vec3.h`: Defines the `Vec3d` class for 3D vector operations.
- `Mat3.h`: Defines the `Mat3` class for 3D matrix operations.
- `geom3D.h`: Provides geometric primitives and functions, including the `Triangle3D` class.
- `raytrace.h`: Likely provides base classes or utilities related to ray tracing.
- `CMesh.h`: Defines the `CMesh` class for managing triangle meshes.
- `TriangleRayTracer.h`: Defines the `TriangleRayTracer` class, used for ray tracing and occlusion calculations in the radiosity algorithm.
- `Lingebra.h`: Provides linear algebra utilities, potentially used for vector and matrix operations.

---

## Types (classes and structs)

### class `Radiosity`

The `Radiosity` class implements the radiosity algorithm, computing the light distribution in a scene by solving a linear system of equations representing the energy exchange between surfaces. It inherits ray tracing capabilities from `TriangleRayTracer` for occlusion calculations and linear system solving functionality from `LinSolver`.

**Inheritance**

- `TriangleRayTracer`: Provides ray tracing functionality for occlusion calculations.
- `LinSolver`: Provides linear system solving functionality.

#### properties

- `couplingTrashold`:`double`:`public:` - A threshold value used to discard negligible coupling coefficients between surface elements in the radiosity matrix.
- `M`:`double*` - A pointer to the radiosity matrix, which represents the coupling coefficients between surface elements.
- `vals`:`double*` - A pointer to an array storing the radiosity values (unknowns) for each surface element. These are the values that the radiosity algorithm solves for. Corresponds to the `a` variable in the `LinSolver` base class.
- `sources`:`double*` - A pointer to an array storing the emission values (light sources) for each surface element. Corresponds to the `b` variable in the `LinSolver` base class.

#### methods

- `allocateWork` - Allocates memory for the `vals` and `sources` arrays, based on the number of surface elements in the scene.
- `prepare` - Prepares the radiosity algorithm for execution by allocating work arrays and setting up the linear problem in the `LinSolver` base class.
- `makeCouplingMatrix` - Computes the radiosity matrix, which represents the geometric relationships and visibility between surface elements. This matrix is a key component of the radiosity algorithm, and its elements determine the amount of light exchanged between surfaces.
- `processTriangles` - Processes a set of triangles by converting them into surface elements and constructing the radiosity matrix.
- `step_Direct` - Performs a single iteration of the radiosity algorithm using a direct method. This method updates the radiosity values for each surface element based on the radiosity matrix and the emission values. This is effectively a Gauss-Seidel iteration. The equation implemented is:

    $$
    B_i = \sum_{j=1}^{n} F_{ij} (B_j + E_j)
    $$

    where \(B_i\) is the radiosity of surface element \(i\), \(F_{ij}\) is the form factor (coupling coefficient) between surface elements \(i\) and \(j\), \(E_j\) is the emittance of surface element \(j\), and \(n\) is the total number of surface elements.  The code is missing a factor in source term, it should be: `vals[i] = Axi + sources[i];`