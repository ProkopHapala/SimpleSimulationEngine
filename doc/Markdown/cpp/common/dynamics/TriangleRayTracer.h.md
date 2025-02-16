# TriangleRayTracer.h

This file defines the `TriangleRayTracer` class, which is designed to manage a collection of triangles and perform ray tracing operations against them. It includes methods to convert triangles into smaller surface elements for more accurate calculations, add triangles to the scene, and determine if a ray is occluded by any of the triangles. It's used in applications such as rendering, collision detection, and visibility analysis.

## Includes

- `Vec3.h`: Defines the `Vec3d` class for 3D vector operations.
- `Mat3.h`: Defines the `Mat3` class for 3D matrix operations.
- `geom3D.h`: Provides geometric primitives and functions, including the `Triangle3D` class.
- `raytrace.h`: Likely provides base classes or utilities related to ray tracing.
- `CMesh.h`: Defines the `CMesh` class for managing triangle meshes.
- `Lingebra.h`: Provides linear algebra utilities, potentially used for vector and matrix operations.

---
## Types (classes and structs)
---

### class `SurfElement`

The `SurfElement` class represents a small surface element derived from a larger geometric primitive, such as a triangle. It stores the position, normal, area, and surface index of the element. These elements are used to approximate the surface of triangles for more accurate ray tracing or other calculations.

#### properties

- `pos`:`Vec3d`:`public:` - The 3D position of the surface element.
- `normal`:`Vec3d` - The 3D normal vector of the surface element, indicating its orientation.
- `area`:`double` - The area of the surface element.
- `isurf`:`int` - The index of the original surface (e.g., triangle) from which this element was derived.

#### methods

- `geomCoupling` - Calculates the geometric coupling between the surface element and another `SurfElement`, based on their positions, normals, and areas. The other `SurfElement` is passed as argument.
- `geomCoupling` - Calculates the geometric coupling between the surface element and another point and normal, based on their positions, normals, and areas. The formula implemented is:

    $ C = \frac{c_0 * c}{r^2 * (r^2 + A_0)} $

    where $ c_0 = \vec{d} \cdot \vec{n_0}$ , $ c = \vec{d} \cdot \vec{n}$, $ \vec{d} = \vec{p} - \vec{p_0} $, $r = ||\vec{d}||$, and $A_0$ is a combined area term.  This is used to model interactions between surface elements, and the coupling is inversely proportional to the squared distance and incorporates a correction for small sizes using the area term.

### class `TriangleRayTracer`

The `TriangleRayTracer` class manages a collection of triangles and provides methods for determining if a ray intersects with any of these triangles. It converts triangles into smaller `SurfElement` objects for increased accuracy in ray tracing and occlusion calculations.

#### properties

- `triangleObstacles`:`std::vector<Triangle3D>`:`public:` - A vector storing the 3D triangles that serve as obstacles for ray tracing.
- `elements`:`std::vector<SurfElement>` - A vector storing the surface elements derived from the triangles. These elements are used for more detailed calculations.

#### methods

- `clearTriangles` - Clears both the `triangleObstacles` and `elements` vectors, effectively removing all triangles and surface elements from the scene.
- `trinagleToElements` - Converts a single `Triangle3D` object into a set of `SurfElement` objects, distributing them evenly across the triangle's surface. It subdivides the triangle into smaller elements based on the input parameter `n`.
- `trapezElements` - Generates surface elements within a trapezoidal region defined by vectors and parameters. It's used as a helper function in `trinagleToElements2` to generate elements more efficiently.
- `trinagleToElements2` - Converts a `Triangle3D` object into a set of `SurfElement` objects, adapting the element size to the triangle's dimensions and using `trapezElements` for efficient generation. This method provides an alternative approach to `trinagleToElements` that may be more efficient for certain triangle shapes.
- `addTriangle` - Adds a `Triangle3D` object to the `triangleObstacles` vector and, if `active` is true, converts it into `SurfElement` objects using `trinagleToElements2`.
- `fromMesh` - Adds triangles from a raw mesh representation (vertices and triangle indices) to the `triangleObstacles` vector and converts them into `SurfElement` objects.
- `fromMesh` - Adds triangles from a `CMesh` object to the `triangleObstacles` vector and converts them into `SurfElement` objects.
- `getOcclusion` - Determines if a ray, defined by its origin (`ray0`) and direction (`hRay`), is occluded by any of the triangles in the scene within a specified maximum distance (`tmax`).  It ignores intersections with triangles specified by indices `ip1` and `ip2`. Returns 1.0 if occluded, 0.0 otherwise.