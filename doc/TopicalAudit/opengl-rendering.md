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
| C++ | `cpp/common_SDL/SDL2OGL3/MeshRenderOGL3.h/.cpp` | active | VAO-based mesh renderer: solid (lit) + line (const-color) shaders, `uploadMesh_d()` / `uploadLines_d()`, model transform, camera setup. Wraps `GLMesh` with VAOs. |
| C++ | `cpp/common_SDL/SDL2OGL3/TextRendererOGL3.h` | active | Instanced billboarded 3D text: VAO + static quad + dynamic instance VBO (`GL_STREAM_DRAW` with orphaning), `addLabel()`, `clear()`, `draw(cam)`. Inline vert/frag shaders. `makeDummyFontTex()` helper. Binary alpha discard — pixelated at zoom. |
| C++ | `cpp/common_SDL/SDL2OGL3/DebugView.h` | active | Global debug rendering: `DEBUG_mesh` (`GLMeshBuilder*`) + `DEBUG_shader`. `DEBUG_drawCamera()`, `DEBUG_draw()`. Macros `DEBUG_VIEW_INIT()` / `DEBUG_VIEW_DEFINE()`. **Inefficient: `new`/`delete` GLMesh every frame.** |
| C++ | `cpp/common_SDL/SDL2OGL3/GLobjects.h` | active | `GLMesh`: VBO-based mesh (pos/nor/col/UV + index), `GL_STATIC_DRAW` default, `init_d()` for double input, `init_wireframe()`, `init_hardface()`. `preDraw()`/`draw()` re-binds VBOs every call (no VAO). |
| C++ | `cpp/common_SDL/SDL2OGL3/DrawOGL3.h` | active | `GLMeshBuilder`: CPU-side builder with `addLine()`, `addPointCross()`, `addArrow()`, `addLines()`, transforms (`move`, `rotate`, `applyMatrix`), `makeLineMesh()` / `makeGLmesh()`. UV-mesh templates (`UVFunc2smooth`, `UVFunc2wire`, `drawExtrudedWireUVFunc`). Primitive generators (`Sphere2Mesh`, `Cone2Mesh`, `Torus2Mesh`, etc.). |
| C++ | `cpp/common_SDL/SDL2OGL3/GLInstances.h` | active | `GLInstances`: instanced meshes (pos/dir/up/scale per instance). `GLBillboards`: instanced billboards. `GLParticles`: particle system with streaming pos+color. **None use VAOs — re-bind VBOs every call.** |
| C++ | `cpp/common_SDL/SDL2OGL3/GLfunctions.h` | active | GL helpers: `newArrayBuffer()`, `uploadArrayBuffer()` (orphaning), `bindVertexAttribPointer()`, `newElementBuffer()`, `newTexture2D()`, `bindTexture()`, `checkOpenGLError()` / `GL_DEBUG` macro. |
| C++ | `cpp/common_SDL/SDL2OGL3/SceneOGL3.h` | active | `setCamera()` / `setCameraPersp()` / `setCameraOrtho()` / `useWithCamera()` helpers. `SceneNode3D` + `SceneOGL3` scene graph base. |
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
| C++ | `cpp/common_SDL/SDL2OGL3/SDFFont.h` | active | SDF font rendering: TTF atlas generation (SDL2_ttf + Felzenszwalb-Huttenlocher distance transform), instanced billboard renderer with `dFdx`/`dFdy` screen-space adaptive `smoothstep` edges. Crisp at any zoom from 64x64 textures. |

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

### Signed Distance Field (SDF) Rendering

**SDF Font Rendering (implemented):**
- `cpp/common_SDL/SDL2OGL3/SDFFont.h` — complete SDF font pipeline:
  - `makeSDFFontAtlas()` — renders TTF glyphs via SDL2_ttf at 8x resolution, computes two-pass Felzenszwalb-Huttenlocher distance transform (inside→outside + outside→inside), normalizes to [0,1] with 0.5 at edge, downsamples to 64x64 per glyph, uploads as GL_R32F horizontal strip atlas
  - `SDFFontRenderer` — instanced billboard renderer (same pattern as `TextRendererOGL3`), fragment shader uses `dFdx`/`dFdy` screen-space derivatives for adaptive `smoothstep` edge width (~1 texel wide at any zoom)
  - `makeDummySDFFontAtlas()` — analytical rectangle SDF for pipeline testing without a font file
- `cpp/sketches_SDL/OGL3/test_SDFFont.cpp` — interactive demo with real DejaVu Sans font
- Links against SDL2_ttf (added to `CMakeLists.txt`)

**SDF Shape Rendering (existing, shapes not fonts):**
- `python/GLCL2/shaders/sdf_gen.glslf` — generates SDF for a circle into RGBA32F texture (texel units)
- `python/GLCL2/shaders/sdf_view.glslf` — views SDF as grayscale (visualization, no alpha-test rendering)
- `python/GLCL2/scripts/sdf.py` / `sdf2.py` — demo configs for SDF generation pipeline
- `python/subPixelContour/subPixelContour.py` — bilinear fitting for sub-pixel contours (related concept)

**Key paper:** Chris Green (Valve), "Improved Alpha-Tested Magnification for Vector Textures and Special Effects", SIGGRAPH 2007. Extension: msdfgen (Viktor Chlumský) — 3-channel SDF for sharper corners.

**Comparison: SDF vs Binary Alpha Text Rendering**

| Feature | `TextRendererOGL3` (binary) | `SDFFontRenderer` (SDF) |
|---|---|---|
| Atlas resolution | High (e.g. 256x256/glyph) | Low (64x64/glyph) |
| Edge quality at zoom | Pixelated, blocky | Crisp, anti-aliased |
| Fragment shader | `if (alpha < 0.1) discard;` | `smoothstep(0.5-sm, 0.5+sm, dist)` |
| Adaptive smoothing | None | `dFdx`/`dFdy` screen-space derivatives |
| Memory | Large atlas | ~16x smaller atlas |
| Dependencies | None (dummy texture) | SDL2_ttf + TTF font file |

## Parity Status

- **C++ `Draw3D` (legacy) ↔ C++ `Draw3D_o3` (Michal)**: Two parallel implementations with overlapping functionality. Legacy uses `opengl1renderer`, Michal uses `GLMesh`/`MeshLibrary`. No convergence.
- **C++ `Draw3D` ↔ JS `Draw3D.js`**: JavaScript port exists but coverage unknown.
- **C++ `GLObject` ↔ Python `GLCLWidget`**: Different VBO/VAO management approaches. No parity test.
- **C++ `GLFramebuffer` ↔ Python `OGLSystem`**: Both manage FBOs but different APIs.

## Efficiency Evaluation & Improvement Suggestions

### What's done well

- **`MeshRenderOGL3`** uses VAOs properly — binds VAO once at upload, then `glDrawElements` with just `glBindVertexArray` at draw time. Correct modern pattern.
- **`TextRendererOGL3`** is well-designed: instanced rendering, VAO-based, buffer orphaning for streaming, `GL_STREAM_DRAW` for dynamic data, `GL_STATIC_DRAW` for the quad.
- **`GLMesh`** uses `GL_STATIC_DRAW` by default for VBOs — correct for static meshes.
- **`GLInstances` / `GLBillboards`** use `glVertexAttribDivisor` for instancing — good for many objects.

### Key inefficiencies

#### 1. `GLMesh::draw()` re-binds VBOs every call (no VAO)

`GLobjects.h:93-114` — `preDraw()` calls `bindVertexAttribPointer()` for each attribute every draw call. A VAO caches all this state — bind once, then just `glBindVertexArray(vao)` + `glDrawElements()`. `MeshRenderOGL3` does this correctly, but raw `GLMesh::draw()` does not. Every sketch using `GLMesh` directly pays this cost.

#### 2. `DebugView.h` — `new`/`delete` GLMesh every frame

`DebugView.h:52-61` — `DEBUG_draw()` calls `makeLineMesh()` (which does `new GLMesh()` + `glGenBuffers`), draws, then `delete msh` (which does `glDeleteBuffers`). This allocates/frees GPU buffers every frame. Should use a **persistent streaming VBO** — pre-allocate a large buffer, `glBufferSubData` or `glMapBuffer` each frame, never delete.

#### 3. `GLMeshBuilder` → `GLMesh` pipeline designed for build-once, not streaming

`DrawOGL3.h:105-113` — `makeGLmesh()` / `makeLineMesh()` copy CPU `std::vector` data to a new `GLMesh` with `glBufferData`. Fine for one-time static mesh creation. But `DebugView` does this every frame — CPU vector → GPU buffer → delete.

#### 4. No UBO (Uniform Buffer Objects) for camera

Every shader (`Shader.h`) has `set_camMat()`, `set_camPos()`, `set_modelPos()`, `set_modelMat()` called individually per draw call. With UBOs, you bind camera data once per frame to a shared buffer all shaders reference via `layout(std140, binding=0) uniform CameraBlock`. This reduces per-draw uniform calls.

#### 5. No shader caching / management

Each renderer creates its own `Shader` instances. `MeshRenderOGL3` has `sh_solid` + `sh_const`, `TextRendererOGL3` has its own `shader`, `DebugView` has `DEBUG_shader`. If two renderers need the same shader, it gets compiled twice. A simple `std::map<string, Shader*>` cache (`ShaderManager`) would fix this.

#### 6. `GLInstances` / `GLBillboards` / `GLParticles` don't use VAOs

`GLInstances.h:33-53` — `draw()` calls `bindVertexAttribPointer()` for 6 attributes every call. Same problem as raw `GLMesh`. Should wrap in VAO.

### Recommended improvements (priority order)

| Priority | Fix | Impact |
|---|---|---|
| **High** | Persistent streaming VBO for debug lines (no per-frame alloc) | Eliminates GPU buffer churn |
| **High** | Add VAO to `GLMesh` (optional method `setupVAO()`) | Removes per-draw VBO re-binding |
| **Medium** | UBO for camera uniforms | One bind per frame instead of per-shader-per-draw |
| **Medium** | Shader cache (`ShaderManager` with `get("color3D")`) | Avoid duplicate compilations |
| **Low** | VAO for `GLInstances` / `GLBillboards` | Same as GLMesh fix |
| **Low** | `glMultiDrawElements` for batched static meshes | Fewer draw calls for multi-sub meshes |

### Proposed unified renderer design

A single `RendererOGL3` class wrapping existing components:

```
RendererOGL3
├── ShaderManager        // cache by name, load from file or inline
├── CameraUBO            // one UBO, shared by all shaders (std140 layout)
├── DebugLineBatch       // persistent streaming VBO, addLine()/flush()
├── StaticMeshStore      // upload once, VAO-based draw
└── TextRendererOGL3     // already good, just plug in
```

This gives the 4 desired features (debug lines, efficient static meshes, shader management, text rendering) with minimal new code — mostly gluing existing components and fixing the debug line streaming.

## Open Issues

- Two parallel `Draw3D` implementations (legacy vs Michal) — DRY violation
- `GLObject::setup()` prints debug `printf` for every vertex in `GLBuff::setup(double*)` overload
- `GLFramebuffer::begin()` calls `exit(0)` on error — should throw
- `Draw3D` text rendering (`drawTextBillboard`) always delegates to billboard regardless of 3D text call
- `Camera` has no orthographic mode in all variants
- No modern OpenGL core profile in legacy `Draw3D` — uses `opengl1renderer` abstraction
- `GLMesh` static `tmpMesh` buffers — not thread-safe
- No depth peeling or order-independent transparency
- `GLMesh::draw()` re-binds all VBOs every call — no VAO (see Efficiency section)
- `DebugView.h` allocates/frees GPU buffers every frame — needs persistent streaming VBO
- No shader cache — same shader compiled multiple times across renderers
- No UBO for camera — per-draw uniform calls instead of per-frame bind
- `GLInstances` / `GLBillboards` / `GLParticles` — no VAOs, re-bind VBOs every draw
- ~~No SDF font rendering~~ — **Implemented**: `SDFFont.h` with `makeSDFFontAtlas()` + `SDFFontRenderer` (see SDF section above)
- `TextRendererOGL3` still uses binary alpha discard — consider migrating to SDF renderer
- SDF renderer uses single-channel SDF only — multi-channel (msdfgen) not yet implemented for sharper corners
- SDF atlas generation is CPU-only (distance transform) — could be GPU-accelerated with compute shader
