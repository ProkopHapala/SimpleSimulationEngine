---
type: TopicalAudit
title: OpenGL Rendering
tags: [topic, cpp, javascript, opengl, glsl, shader, draw3d, camera, framebuffer, vbo, texture]
---

## Summary

OpenGL rendering infrastructure across C++ (SDL2 + GLEW) and JavaScript (Three.js / WebGL). C++ side: `Draw3D` namespace for immediate-mode-style 3D primitives (points, lines, arrows, spheres, cones, boxes, text), `GLObject`/`GLBuff` for VBO/VAO management, `Shader` for GLSL program compilation, `GLFramebuffer` for render-to-texture, `Camera` for view/projection matrices. JavaScript side: Three.js-based rendering for MHD plasma, LandCraft WebGPU terrain, and N-body demos.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common_SDL/SDL2OGL3/Draw3D.h` / `.cpp` | active | `Draw3D` namespace: `drawPoint()`, `drawLine()`, `drawArrow()`, `drawTriangle()`, `drawSphere()`, `drawCone()`, `drawCylinderStrip()`, `drawBox()`, `drawBBox()`, `drawPolyLine()`, `drawTextBillboard()`, `drawScalarGrid()`, `drawAxis3D()`. Double-precision overloads (cast to float). `toGLMat()` / `toGLMatCam()` for rigid transforms. `drawPBC()` templated for periodic boundary visualization. |
| C++ | `cpp/common_SDL/SDL2OGL3/Michal/Draw3D_o3.cpp` | active | Alternative `Draw3D` using `GLMesh` + `MeshLibrary` (shader-based). Static `tmpMesh1/2/3` for dynamic geometry. Uses `setUniform3f()` for color. |
| C++ | `cpp/common_SDL/SDL2OGL3/Shader.h` | active | `Shader` class: vertex/fragment shader compilation, `setUniform*()` methods, program linking |
| C++ | `cpp/common_SDL/SDL2OGL3/GLObject.h` / `.cpp` | active | `GLObject`: VBO/VAO management with up to 4 attribute buffers (vertex, normal, color, UV). `GLBuff`: single buffer abstraction with `toGPU()` and `activate()`. Indexed drawing support. |
| C++ | `cpp/common_SDL/SDL2OGL3/Camera.h` | active | `Camera`: position, rotation (Mat3f/Quat4f), perspective/orthographic projection, `setGLProjection()`, `setGLModelView()`. |
| C++ | `cpp/common_SDL/SDL2OGL3/Michal/Framebuffer.h` | active | `GLFramebuffer`: render-to-texture with color + depth attachments. `begin()`/`end()`/`pause()`/`unpause()` for nested rendering. `GLTexture` for texture management. |
| C++ | `cpp/common_SDL/SDL2OGL3/MeshRenderOGL3.cpp` | active | Mesh rendering with materials and shaders |
| C++ | `cpp/common_SDL/SDL2OGL3/SceneOGL3.h` | active | Scene graph with camera, objects, lighting |
| C++ | `cpp/common_SDL/SDL2OGL3/Michal/GLMesh.h` | active | Templated `GLMesh<Attrs...>`: dynamic vertex buffers with `addVertex()`, `clear()`, `draw()`. Attribute types: `MPOS`, `MNORMAL`, `MCOLOR`. |
| C++ | `cpp/common_SDL/SDL2OGL3/Michal/MeshLibrary.h` | active | Library of pre-built meshes (point, line, pointCross) with associated shaders |
| C++ | `cpp/GUI/Draw3D.c` / `Draw3D.h` | active | C version of Draw3D for simpler contexts |
| C++ | `cpp/GUI/Camera.c` / `Camera.h` | active | C version of camera |
| JS | `js/mhd_demo/render.js` | active | Three.js renderer for MHD plasma: orthographic camera, coil rings, B-field vector grid, GLSL fragment shader for HSV-colored B-field. See `mhd-plasma.md`. |
| JS | `js/common_js/Draw3D.js` | active | JavaScript port of Draw3D for WebGL |
| JS | `js/LandCraft_web/LandCraft.html` | active | WebGPU/WGSL-based terrain rendering. See `landcraft.md`. |
| Python | `python/GLCL/NBody_glcl.py` | active | PyQt5 OpenGL widget: VBO/VAO, vertex/fragment shaders, point rendering. See `pyopencl-workflows.md`. |
| Python | `python/GLCL2/GLCLGUI.py` | active | `GLCLWidget`: QOpenGLWidget with CL↔GL sharing, texture/FBO display |
| Python | `python/GLCL2/OGLsystem.py` | active | `OGLSystem`: OpenGL context, shader compilation, texture/FBO management |
| Doc | `doc/Markdown/GUI.md` | doc | GUI documentation |
| Doc | `doc/Markdown/cpp/` | doc | C++ rendering documentation |

## Sub-topics

### Draw3D Primitives

Two implementations:
- **Legacy (`Draw3D.h/.cpp`)**: Uses `opengl1renderer` (fixed-function-style API with shader backend). Immediate-mode-style calls. `drawPBC()` template for periodic boundary visualization. Double overloads cast to float.
- **Michal (`Draw3D_o3.cpp`)**: Uses `GLMesh<MPOS,...>` with `MeshLibrary` shaders. Static `tmpMesh` buffers for dynamic geometry. `setUniform3f("uColor", ...)` for per-draw color.

### VBO/VAO Management

`GLObject`:
- Up to 4 `GLBuff` attributes (vertex, normal, color, UV)
- `setup(nVert)`: allocates vertex + normal arrays
- `setIndexes(nInd, cbuff)`: creates index buffer
- `draw()`: binds VAO, sets attributes, draws indexed
- `draw_instance()`: instanced rendering

`GLBuff`:
- `toGPU(n)`: `glGenBuffers` → `glBufferData`
- `activate()`: `glEnableVertexAttribArray` → `glVertexAttribPointer`

### Framebuffer (Render-to-Texture)

`GLFramebuffer`:
- Color + depth texture attachments
- `begin()`: push framebuffer, clear, set viewport
- `end()`: pop framebuffer
- `pause()`/`unpause()`: temporarily switch to default framebuffer
- Status check with error strings
- Auto-resize to screen size

### Camera

`Camera`:
- Position (`Vec3f`), rotation (`Mat3f` or `Quat4f`)
- Perspective: fov, aspect, near/far
- `setGLProjection()`: sets projection matrix
- `setGLModelView()`: sets modelview from camera pose
- `toGLMatCam()`: converts pose to GL matrix (transposed + negated for camera)

### GLSL Shaders

- `Shader.h`: compilation, linking, uniform setting
- `MeshLibrary`: pre-built shader programs for point, line, pointCross
- MHD plasma: fragment shader with elliptic integrals for B-field visualization
- LandCraft: WGSL compute shaders for terrain generation

## Parity Status

- **C++ `Draw3D` (legacy) ↔ C++ `Draw3D_o3` (Michal)**: Two parallel implementations with overlapping functionality. Legacy uses `opengl1renderer`, Michal uses `GLMesh`/`MeshLibrary`. No convergence.
- **C++ `Draw3D` ↔ JS `Draw3D.js`**: JavaScript port exists but coverage unknown.
- **C++ `GLObject` ↔ Python `GLCLWidget`**: Different VBO/VAO management approaches. No parity test.
- **C++ `GLFramebuffer` ↔ Python `OGLSystem`**: Both manage FBOs but different APIs.

## Open Issues

- Two parallel `Draw3D` implementations (legacy vs Michal) — DRY violation
- `GLObject::setup()` prints debug `printf` for every vertex in `GLBuff::setup(double*)` overload
- `GLFramebuffer::begin()` calls `exit(0)` on error — should throw
- `Draw3D` text rendering (`drawTextBillboard`) always delegates to billboard regardless of 3D text call
- `Camera` has no orthographic mode in all variants
- No modern OpenGL core profile in legacy `Draw3D` — uses `opengl1renderer` abstraction
- `GLMesh` static `tmpMesh` buffers — not thread-safe
- No depth peeling or order-independent transparency
