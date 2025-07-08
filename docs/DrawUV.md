# DrawUV.h Function Documentation

## Function Templates (General Mesh Topology)

- **`UVFunc2smooth(...)`** - Creates a smooth mesh surface from a UV function by sampling points on a grid and generating triangles with normals.
- **`UVFunc2wire(...)`** - Generates a wireframe mesh from a UV function by sampling points and connecting them with edges in a grid pattern.
- **`UVFunc2wireExtruded(...)`** - Creates an extruded wireframe mesh from a UV function by generating two layers of points (base and offset by normal) and connecting them with edges and faces.

## Specialized Functions (Specific UV Surfaces)

- **`Cone2Mesh(...)`** - Creates a cone mesh with given radii and length.
- **`Sphere2Mesh(...)`** - Generates a sphere mesh with given radius.
- **`Torus2Mesh(...)`** - Creates a torus (donut) mesh with given minor and major radii.
- **`Teardrop2Mesh(...)`** - Generates a teardrop-shaped mesh with given radii and length.
- **`NACASegment2Mesh(...)`** - Creates a NACA airfoil segment mesh using given coefficients.
- **`HarmonicTube2Mesh(...)`** - Generates a tube with harmonic oscillations along its length.
- **`Parabola2Mesh(...)`** - Creates a parabolic surface mesh.
- **`Hyperbola2Mesh(...)`** - Generates a hyperbolic surface mesh.
- **`Cone_ExtrudedWire(...)`** - Creates an extruded wireframe version of a cone.
- **`Sphere_ExtrudedWire(...)`** - Generates an extruded wireframe version of a sphere.
- **`Torus_ExtrudedWire(...)`** - Creates an extruded wireframe version of a torus.
- **`Teardrop_ExtrudedWire(...)`** - Generates an extruded wireframe version of a teardrop shape.
- **`NACASegment_ExtrudedWire(...)`** - Creates an extruded wireframe version of a NACA airfoil.
- **`HarmonicTube_ExtrudedWire(...)`** - Generates an extruded wireframe version of a harmonic tube.
- **`Parabola_Wire(...)`** - Creates a basic wireframe version of a parabola.
- **`Parabola_ExtrudedWire(...)`** - Generates an extruded wireframe version of a parabola.
- **`Hyperbola_ExtrudedWire(...)`** - Creates an extruded wireframe version of a hyperbola.
- **`ConeFan(...)`** - Generates a cone-shaped fan of triangles.
- **`CylinderStrip(...)`** - Creates a cylindrical strip mesh between two points.
- **`CylinderStrip_wire(...)`** - Generates a wireframe cylindrical strip between two points.
- **`SphereTriangle_wire(...)`** - Creates a wireframe spherical triangle.
- **`SphereTriangle(...)`** - Generates a filled spherical triangle.
- **`Sphere_oct(...)`** - Creates a sphere by combining 8 spherical triangles (octahedral approximation).
- **`Capsula(...)`** - Generates a capsule mesh (cylinder with hemispherical caps) between two points.

## New Flexible Wireframe Functions

- **`UVFunc2wire_new(...)`** - Enhanced wireframe generator with:
  - `closeCenter`: Collapses center to single vertex (default false)
  - `periodicB`: Closes angular dimension (default false)
- **`Parabola_Wire_new(...)`** - Flexible parabola wireframe with:
  - `closeCenter`: Collapses vertex to single point (default true)
  - `periodicB`: Closes equatorial angle (default true)
