# raytrace.h

This file provides a collection of functions for performing ray tracing operations against various geometric primitives, including spheres, lines, planes, triangles, polygons, and bounding boxes. These functions are used to determine if a ray intersects with a given primitive and, if so, to calculate the point of intersection and the surface normal at that point. This is a fundamental component of ray tracing-based rendering engines and other applications that require simulating the interaction of light with objects in a 3D scene.

## Includes

- `<math.h>`: Provides standard mathematical functions, such as `sqrt` and `fabs`.
- `<cstdlib>`: Provides general utility functions, including memory allocation and random number generation.
- `<stdio.h>`: Provides standard input/output functions, such as `printf`.
- `Vec2.h`: Defines the `Vec2` class for 2D vector operations.
- `Vec3.h`: Defines the `Vec3` class for 3D vector operations.

---

## Free functions

- `setToRay` - Modifies a point to lie on a ray, by removing the component orthogonal to the ray direction.
- `rayPointDistance2` - Calculates the squared distance between a ray and a point, also computing the parameter `t` along the ray closest to the point.
- `linePointDistance2` - Calculates the squared distance between a line segment and a point.
- `raySphere` - Calculates the distance along a ray to the intersection point with a sphere. Returns `INFINITY` if no intersection occurs.
- `sphereNormal` - Calculates the normal vector on a sphere at a given point of intersection with a ray.
- `rayLineDist` - Calculates the distance between a ray and a line.
- `rayLine` - Calculates the closest points between a ray and a line, returning the parameters `t1` and `t2` along the ray and line, respectively, and the distance between the lines.
- `capsulaIntersect` - Calculates the intersection distance between a ray and a capsule (a cylinder with hemispherical caps).
- `rayPlane` - Calculates the distance along a ray to the intersection point with a plane. Returns `t_inf` if no intersection occurs (ray is parallel to the plane).
- `pointInTriangleEdges` - Determines if a point lies within a triangle, by checking if it is on the same side of each edge as the third vertex.
- `rayTriangle` - Calculates the distance along a ray to the intersection point with a triangle, also determining if the intersection point is inside the triangle and calculating the surface normal.
- `rayInTriangle` - Determines if a point (expressed as a vector from a triangle vertex) lies inside a triangle, by projecting the vectors onto a plane and checking the sign of the cross products.
- `rayTriangle2` - Calculates the distance along a ray to the intersection point with a triangle, using pre-computed orthogonal vectors to the ray direction for inside testing.
- `rayTriangle2` - Calculates the distance along a ray to the intersection point with a triangle, automatically generating orthogonal vectors to the ray direction.
- `rayPolygon` - Calculates the distance along a ray to the intersection point with a polygon, using pre-computed orthogonal vectors to the ray direction for inside testing.
- `rayTriangles` - Finds the closest intersection point between a ray and a set of triangles.
- `pickParticle` - Finds the closest particle (sphere) to a ray.
- `pickPoinMinDist` - Finds the closest point to a ray.
- `pickBondCenter` - Finds the closest bond center (midpoint between two points) to a ray.
- `rayPickBond` - Finds the closest bond (defined by a function providing endpoints) to a ray.
- `rayBox` - Calculates the distance along a ray to the intersection point with a bounding box.