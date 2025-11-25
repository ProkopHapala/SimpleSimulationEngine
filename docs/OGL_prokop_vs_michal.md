# Comparison of OpenGL3 Rendering Systems

_Work in progress notes._

## Scope

This document compares the original OpenGL3 system (`cpp/common_SDL/SDL2OGL3`) with Michal's system (`cpp/common_SDL/SDL2OGL3/Michal`) in terms of:

- **Completeness / generality**
- **Performance efficiency**
- **Elegance / simplicity**

Findings will be refined iteratively as more of the code is inspected.

## High-level architecture

- **Prokop system (SDL2OGL3 root)**
  - **App / window / input**:
    - `AppSDL2OGL3`, `ScreenSDL2OGL3` manage SDL init, GL context, camera, and the main loop.
    - `ScreenSDL2OGL3` also sets up a deferred framebuffer and a full-screen quad for post‑processing.
  - **Scene graph / world structure**:
    - `SceneOGL3`, `SceneNode3D` provide a very lightweight scene abstraction: a VAO plus a virtual `draw(const Camera&)` that you override.
    - Higher‑level drawing for particular domains is split into headers such as `DrawOGL3.h`, `DrawUV.h`, `DrawField.h`, `DrawVegetation.h`, `TerrainOGL3.h`.
  - **Rendering primitives**:
    - `GLMesh`, `GLObject`, `GLObject.h`, `GLobjects.h`, `GL3Utils.h` expose low‑level VBO/VAO/texture/FBO helpers and a lot of mesh/terrain/grid generators.
    - `GLfunctions.h` contains raw GL utilities and error‑checking wrappers.
  - **Instancing, sprites, examples**:
    - `GLInstances.h` defines `GLInstances`, `GLBillboards`, `GLParticles` for various instanced layouts.
    - Sketches under `cpp/sketches_SDL/OGL3/` (e.g. `test_Instances.cpp`, `test_Atoms.cpp`, `test_Sprites.cpp`, `test_RenderToTexture.cpp`, `test_AntiAliasing.cpp`) demonstrate geometry instancing, billboards/sprites, and render‑to‑texture.
  - **Shader handling**:
    - A single `Shader` class with hard‑coded default GLSL and simple `setUniform*` helpers; more advanced shaders are typically managed directly by file path and manual uniform calls.

- **Michal system (`SDL2OGL3/Michal`)**
  - **Central GL state (`GLES`)**:
    - `GLES.cpp/.h` wrap GL calls and track current program, bound buffers, framebuffer stack, active camera, and screen size.
  - **Shader system**:
    - `Shader.h` wraps program creation, attribute lookups, and **uniform batching** via `GLuniform`/`GLuniformSet` with type‑safe enums and bitmasks.
    - Default shader source is generated based on requested attributes (position, normal, color, UV, etc.), so most simple meshes can use auto‑generated GLSL.
  - **Vertex formats, VBOs, meshes**:
    - `GLattribs.h`, `GLvbo.h` define compile‑time vertex layouts and map C++ types to GL enums and attribute sizes.
    - `GLMesh.h` owns a VAO + `GLvbo`, knows its `Shader`, and exposes `draw`, `draw2D`, and transform‑aware draws that combine with `GLES::active_camera`.
    - `GLInstancedMesh.h` generalizes instancing to arbitrary combinations of per‑vertex and per‑instance attributes.
  - **Textures & framebuffers**:
    - `GLTexture.h` manages textures with **lazy initialization** from file or size/format; sampler parameters can be queued before GL context exists.
    - `Framebuffer.h/.cpp` (`GLFramebuffer`) wraps an FBO with color/depth textures sized to `GLES::screen_size` and enforces correct `begin/end/pause/unpause` usage via a stack.
    - `FullscreenShader.h` provides a small helper for post‑processing passes.
  - **Debug drawing, primitives, and text**:
    - `Renderer.{h,cpp}` implements `OpenGL1Renderer opengl1renderer`, an OpenGL 1‑style layer implemented on top of `GLMesh` and internal matrix stacks.
    - `MeshLibrary.{h,cpp}` defines a set of reusable primitives (points, lines, crosses, wire cubes, cubes with normals, circles, rectangles, spheres, instanced spheres).
    - `Draw_o3.{h,cpp}` and `Draw3D_o3.{h,cpp}` provide many immediate‑mode‑like helpers using `MeshLibrary` and `OpenGL1Renderer` (arrows, axes, scalar fields, boxes, triclinic cells, etc.).
    - `TextRenderer.{h,cpp}` implement a small bitmap‑font subsystem on top of `GLMesh<MPOS,MUV>` and a font atlas texture, with `Draw_o3`/`Draw3D_o3` providing convenient 2D and 3D billboarded text helpers.

At a high level, both systems cover similar responsibilities (window/camera, shaders, meshes, textures, framebuffers, some form of scene/drawing helpers), but with different emphasis:

- Your side favors **direct C‑style utilities and specialized helpers** (especially for instancing, sprites, and deferred rendering sketches) with maximum low‑level control.
- Michal’s side focuses on **typed, RAII‑style wrappers** with a central `GLES` state, auto‑generated shaders, and a rich library of debug drawing and text utilities built on top.

In the following sections we will go over:

1. Surface API comparison (meshes, shaders, textures, framebuffers, scene/camera integration).
2. Completeness/generality: which features each side provides and where one covers cases the other doesn't.
3. Performance aspects: buffer management, uniform updates, draw call overhead, texture/FBO handling.
4. Elegance/simplicity: API clarity, separation of concerns, and how easy it is to use/extend.

## 1. Surface API comparison

### Mesh / geometry API

- **Prokop**
  - `GLMesh` (in `GLobjects.h`) is a thin POD-like holder of up to 4 VBOs (`vpos`, `vnor`, `vcol`, `vUVs`) plus optional index buffer.
  - Creation is usually: compute arrays externally → `init(...)` uploads data.
  - Drawing is: `preDraw()` to bind buffers + enable attribute arrays, then `drawRaw()` or `draw(draw_mode)` to call `glDrawArrays`/`glDrawElements`.
  - Higher-level geometry helpers live in `GL3Utils.h` and create `GLMesh` or `GLObject` instances for specific patterns (grids, patches, camera-frustum meshes, etc.).

- **Michal**
  - `GLvbo<attribs...>` owns CPU-side vectors of `GLvertex<...>` with compile-time attribute layout. It knows how to set up `glVertexAttribPointer` based on attribute list.
  - `GLMesh<attribs...>` owns a VAO + a `GLvbo<attribs...>*` and a `Shader*`.
  - Mesh users add data with strongly-typed calls like `addVertex(position, normal, color, uv)`; `GLvbo` handles buffer updates and attribute formats.
  - Draw calls (`draw`, `draw2D`, `draw(position, rotation, scale)`) combine geometry, transform, and current camera into an MVP and call `glDrawArrays` on the internal VBO.

**Implication:** your API is closer to raw OpenGL and gives full control per call; Michal’s API is more structured and type-safe, with camera and transforms integrated into the mesh object. For quick experiments and custom layouts, your approach is more flexible; for larger projects and reuse, Michal’s abstractions reduce boilerplate and mistakes.

### Shader API

- **Prokop**
  - Single `Shader` class that wraps creation of one program with hard-coded default GLSL for simple colored geometry.
  - Uniform locations for camera/model matrices and base color are stored in members (`uloc_camMat`, `uloc_modelMat`, etc.) and updated via helpers like `set_modelMat`, `set_camMat`, `set_baseColor`.
  - For additional uniforms, user calls `getUloc(name)` + `glUniform*` (or the thin wrappers `setUniform*`).

- **Michal**
  - `Shader` takes GLSL source strings, lazily compiles and links them on first use.
  - Uniforms are represented by `GLuniform`/`GLuniformSet` (strongly-typed union + map) and attached to the `Shader` by location.
  - A bitmask (`update_uniforms`) tracks which uniforms changed; on `use()` only those are re-sent to GL.
  - Texture uniforms are tracked separately in `tex_uniforms` and bound automatically (texture unit selection is handled in one place).
  - Default shader source is generated at compile time via templates, based on which attributes are present (`Position`, `Normal`, `Color`, `UV`, etc.), producing matching vertex + fragment shaders without hand-writing GLSL.

**Implication:** your shader API is minimal and explicit, good for learning and manual control. Michal’s shader/uniform system is more complex internally but much more scalable: less duplicated GLSL, safer uniform handling, and less chance of mismatched attribute layouts.

### Textures and framebuffers

- **Prokop**
  - `GLfunctions.h` supplies procedural helpers (`newTexture2D`, `bindTexture`, `makeRandomTexture`). Responsible code must manage lifetimes and sizes manually.
  - `FrameBuffer` in `GLobjects.h` builds an FBO from given or newly-created textures and provides `bind()` + simple initialization; resize and lifetime are left to the caller.
  - `ScreenSDL2OGL3` additionally sets up a deferred framebuffer (color+depth textures) and a canvas quad + shader for post-processing, but in a fairly ad-hoc way.

- **Michal**
  - `GLTexture` owns one texture, with **lazy** creation either from a BMP file path or a size/format; sampler parameters can be queued before context exists and applied later.
  - `GLFramebuffer` always has a color and depth texture sized to `GLES::screen_size`, and `begin/end/pause/unpause` enforce correct use via a framebuffer stack and error checks.
  - Fullscreen passes use `FullscreenShader` plus `GLMesh`/`GLTexture` abstractions instead of hard-coding quad VBOs.

**Implication:** both sides expose enough control for advanced effects, but Michal’s code has a more consistent RAII-style encapsulation of textures/FBOs, while your code gives more direct procedural control and specialized helpers (e.g. for deferred rendering in `ScreenSDL2OGL3`).


### Geometry instancing

- **Prokop**
  - Uses several specialized instancing helpers in `GLInstances.h`:
    - `GLInstances` for oriented/scaled solids with attributes `model_vpos`, `model_vnor`, `pose_pos`, `pose_dir`, `pose_Up`, `pose_sc`.
    - `GLBillboards` for screen-aligned quads with per-instance position/scale.
    - `GLParticles` for line-strip particles with per-instance position/color.
  - Layout is hard-coded on both sides:
    - C++ binds fixed attribute locations 0–5 and divisors manually.
    - Shaders like `Instance3D.glslv` assume that exact layout and compute rotation/scale there.
  - Very explicit, minimal abstraction, and tuned for your use cases in `test_Instances.cpp` and `test_Atoms.cpp`.

- **Michal**
  - Provides a generic instancing abstraction: `GLInstancedMeshBase<VertVboType, instAttribs...>` in `GLInstancedMesh.h`.
    - `VertVboType` is a `GLvbo<attribs...>` describing per-vertex attributes.
    - `instAttribs...` describe per-instance attributes (e.g. offsets, radii) at the type level.
  - Compile-time checks ensure:
    - Vert and instance attribute sets are disjoint and have no duplicates.
  - A default shader is auto-generated for the union of attributes, and `draw()` uses `glDrawArraysInstanced` with data coming from the typed VBOs.

**Implication:** you already support efficient geometry instancing with a very specific layout for atoms/particles, while Michal generalizes the same idea to arbitrary vertex/instance attribute combinations. Your approach gives maximum visibility and manual control for those demos; Michal’s approach makes it easier to add many different instanced systems (vegetation, debris, projectiles, etc.) without writing new low-level glue for each one.


### Scene graph / drawing API

- **Prokop**
  - Scene-level structure is provided by `SceneOGL3` / `SceneNode3D` and friends:
    - Each scene owns a VAO and defines a virtual `draw(const Camera&)` that you override in derived classes.
    - Higher-level domain-specific drawing is done in many separate headers: `DrawOGL3.h`, `DrawUV.h`, `DrawField.h`, `DrawVegetation.h`, `TerrainOGL3.h`.
  - The “scene graph” is simple but explicit: you typically store objects in vectors, iterate, and issue direct GL calls or use `GLMesh`/`GLObject` helpers.

- **Michal**
  - There is no dedicated scene-graph class in `Michal/`, but there is a rich **debug‑drawing and primitive layer**:
    - `Renderer.{h,cpp}` implements `OpenGL1Renderer opengl1renderer`, an OpenGL 1‑style wrapper built on top of `GLMesh<MPOS,MNORMAL,MCOLOR>` and matrix stacks.
    - `MeshLibrary.{h,cpp}` defines reusable meshes for points, lines, crosses, wire cubes, cubes with normals, circles, rectangles, spheres, and instanced spheres.
    - `Draw_o3.{h,cpp}` and `Draw3D_o3.{h,cpp}` provide many immediate‑mode‑style helpers (axes, arrows, scalar fields, boxes, triclinic cells, etc.) implemented using `MeshLibrary` and `OpenGL1Renderer`.
  - Application code using Michal’s system typically builds its “scene” using these mesh/draw helpers plus the central `GLES::active_camera`, instead of inheriting from a `Scene` base class.

**Implication:** your side has a slightly more explicit scene‑graph abstraction (`SceneOGL3`) integrated with your app framework, whereas Michal relies on composing `GLMesh`/`GLInstancedMesh`/`Draw3D_o3` calls around a global camera. For your style of simulation demos the difference is mostly about taste; both approaches are lightweight and keep you close to raw GL.


### Sprites / billboards

- **Prokop**
  - Dedicated instanced billboard class `GLBillboards` in `GLInstances.h`:
    - Per‑vertex attributes: `model_UVs` for a quad (6 vertices).
    - Per‑instance attributes: `pose_pos` (vec3), `pose_sc` (vec2) for position and scale.
    - Rendered with `glDrawArraysInstanced` and shaders like `Bilboard3D.glslv`.
  - `Bilboard3D.glslv` (JS version, conceptually same as C++ one):
    - Uses camera rotation (`camRot`) and instance position/scale to construct world‑space quad vertices that always face the camera.
  - Example usages: `test_Sprites.cpp`, `test_AntiAliasing.cpp`, `test_RenderToTexture.cpp` show sprite clouds and billboarded quads driven by this system.

- **Michal**
  - `Draw_o3` provides low‑level billboard matrix helpers:
    - `billboardCam()` zeros out the rotation part of the modelview matrix.
    - `billboardCamProj()` derives a screen‑space‑scaled billboard basis from the projection and modelview matrices.
  - `MeshLibrary` defines simple 2D primitives (`line2D`, `rect`, `cross`, `xmark`) and 3D primitives that can be drawn via `Draw3D_o3`.
  - There is no dedicated `GLBillboards`‑equivalent instanced class in the Michal tree; sprite‑like quads would typically be built as:
    - A `GLMesh<MPOS,MUV>` quad plus a texture (or atlas) in `GLTexture`.
    - Optional instancing via `GLInstancedMeshBase` with per‑instance offsets/scales.

**Implication:** you have a very concrete, tested billboard/sprite pipeline (instanced quads + shaders + C++ helpers). Michal’s side has the primitives and matrix utilities to implement billboards but no ready‑made sprite manager; his focus there is more on spheres and geometric debug shapes than sprites.


### Render‑to‑texture / post‑processing

- **Prokop**
  - `FrameBuffer` in `GLobjects.h` wraps an FBO with color/depth textures, but usage is manual.
  - `ScreenSDL2OGL3` sets up a deferred G‑buffer and uses a canvas quad plus a shader to do post‑processing in a fairly ad‑hoc but effective way.
  - Sketches like `test_RenderToTexture.cpp` demonstrate rendering a 3D scene into a texture and then drawing it back on a quad or billboard.

- **Michal**
  - `Framebuffer.{h,cpp}` (`GLFramebuffer`) encapsulates a color and depth texture that are **always kept in sync with `GLES::screen_size`**, with `begin()`/`end()`/`pause()`/`unpause()` managing a framebuffer stack and enforcing correct nesting.
  - `FullscreenShader.h` provides a convenient way to draw full‑screen passes using `GLMesh` and `GLTexture` rather than hand‑rolled quads.
  - `Renderer` + `MeshLibrary` can also be used to draw geometry into FBOs for debug visualizations.

**Implication:** both systems support render‑to‑texture and post‑processing. Your implementation is more custom to your deferred renderer and sketches; Michal’s is more systematic and stateful (RAII‑like FBO wrapper + fullscreen shader abstraction), which is attractive if you plan a more complex post‑processing chain, but gives up a bit of your free‑form control.


### Text rendering

- **Prokop**
  - Currently no dedicated text subsystem in the C++ `SDL2OGL3` tree.
  - You *do* have all the low‑level ingredients needed for text as sprites:
    - `GLBillboards` + instancing for many quads.
    - Texture and texture‑atlas sampling shaders (as used for sprites and minimaps).
  - To get text, you would still need a CPU‑side layer that:
    - Loads a font atlas (bitmap or signed‑distance field).
    - Maps characters to UV rectangles in that atlas.
    - Builds per‑glyph quads (position + UV) and uploads them via `GLBillboards` or a simple mesh.

- **Michal**
  - `TextRenderer` implements a complete bitmap‑font text pipeline:
    - Global `GLTexture text_font` loads `common_resources/dejvu_sans_mono_RGBA_pix-UpDown.bmp` and sets nearest‑neighbor filtering.
    - Internal `GLMesh<MPOS,MUV>` (`mesh`) uses `defaultTexColorShader<MPOS,MUV>`.
    - `set(txt, maxWidth, maxHeight)` builds quads per character using `MeshBuilder::addChar`, performing simple word‑wrap in a screen‑space grid.
    - `draw2D(pos, size, color)` renders text in 2D screen space; `draw3D(pos, color)` (in the original project) uses `GLES::active_camera` to place text in 3D.
  - `Draw_o3` wraps this via:
    - Global `TextRenderer text;`
    - `Draw::drawText(...)` helpers that call `text.set()` and `text.draw2D()`.
  - `Draw3D_o3` adds world‑space labels:
    - `drawTextBillboard(const char* str, Vec3f pos, float sz, bool ontop, int iend)` projects `pos` via `GLES::active_camera->world2Screen()` and then draws text in screen space at that position (optionally forcing it on top).

**Implication:** once these modules are present, Michal’s side clearly wins on *ready‑to‑use* text rendering: it already provides font loading, glyph UVs, 2D/3D layout, and billboarding. Your engine is closer to the metal (better low‑level billboard/sprite infrastructure) but still lacks the thin “string → glyph quads” layer, which you could either implement yourself on top of `GLBillboards` or adapt from Michal’s `TextRenderer` if you decide to integrate it.


## 2. Michal modules & integration notes

- **Current status in this repo**
  - All Michal files are under `cpp/common_SDL/SDL2OGL3/Michal` and, as of now, include:
    - Core wrappers: `GLES.h/.cpp`, `GLattribs.h`, `GLuniform.h`, `GLvbo.h`, `GLMesh.h`, `GLInstancedMesh.h`, `GLTexture.h`, `Framebuffer.h/.cpp`, `FullscreenShader.h`, `Shader.h`.
    - Shared math/camera types from your existing code: `Camera.h`, `Mat4.h`, `Vec3.h`, `quaternion.h`.
  - A pass over Michal’s folder shows **no remaining dangling `#include` headers** inside this project; the earlier missing `GLES.h` has been added.

- **Implications for future integration**
  - To actually *use* Michal’s system in the engine you still need to decide how it plugs into your existing structure:
    - How `GLES::active_camera` and `GLES::screen_size` relate to `ScreenSDL2OGL3::cam` and your window/screen management.
    - Whether post-processing should go through `FullscreenShader`/`GLFramebuffer` or your current deferred-quad code in `ScreenSDL2OGL3`.
  - Michal’s abstractions are built around a central `GLES` state; if you integrate them, you should:
    - Keep using your existing low-level helpers (`GLfunctions.h`, `GLobjects.h`, `GL3Utils.h`, `GLInstances.h`) where you need **maximum control and custom layouts** or already have optimized code.
    - Use Michal’s `GLMesh`/`GLvbo`/`GLInstancedMesh`/`Shader`/`GLTexture`/`GLFramebuffer` **selectively**, where the abstraction reduces boilerplate and mistakes without hiding performance-critical details.
  - From a performance/control perspective, nothing in Michal’s design *forces* extra overhead; the main cost is a bit of additional indirection and state bookkeeping in C++. You can still:
    - Keep tight control over when VBOs are updated (`GLvbo::sync()`),
    - Control draw modes and instance counts,
    - Decide per-subsystem whether to stay close to raw GL (your style) or to use the higher-level wrappers.

