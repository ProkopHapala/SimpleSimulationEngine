# pyRay — SDF Ray Tracer (OpenCL)

Modular signed-distance-field (SDF) ray tracer for rendering. Not a physics simulation — uses implicit SDF primitives, not explicit triangle meshes.

## How It Works

Sphere tracing: `castRay()` iteratively steps through an SDF `sceneDistFunction(pos)` until hit or max distance. One OpenCL thread per pixel. Geometry is **implicit** (SDF primitives), not explicit triangles.

- `rayTrace_basic` kernel: per-pixel ray march, accumulates soft-shadow term `1/(1+80·d²)`, computes hit position + normal via finite differences
- Scene defined by composing SDF primitives (`primitives.cl`) with boolean operations (`operations.cl`)
- Normal computed by central-difference SDF gradient

## Files

- **`scene.py`** — Scene assembly: loads CL kernels, injects scene SDF and user functions via string substitution, parses object definitions
- **`ocl.py`** — OpenCL setup: context, queue, ray generation from camera, buffer management, kernel launch
- **`image.py`** — Image output utilities
- **`common.py`** — Shared defines, camera matrix, paths
- **`GUI.py`** — PyQt5 GUI for interactive scene editing
- **`cl/primitives.cl`** — SDF primitives (sphere, box, plane, etc.)
- **`cl/operations.cl`** — SDF boolean ops (union, subtract, intersect, smooth blend)
- **`cl/rayScene_basic.cl`** — Main ray marching kernel

## Relation to Radiosity / pyScatter

pyRay uses **SDF ray marching** — fundamentally different from the explicit triangle intersection used by Radiosity and pyScatter. It is not directly applicable to the radiosity/scattering occlusion problem. However, the camera/ray-generation utilities (`getRays`, `getCamMat`) could be reused for detector ray generation in pyScatter.

## Run

```bash
python GUI.py
```
