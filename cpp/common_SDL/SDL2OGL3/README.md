SDL2OGL3 — OpenGL 3.x rendering infrastructure for SDL2 applications.

Core application framework, scene graph, camera, shader management, and
specialized renderers (mesh, text, SDF font, terrain, vegetation, debug).

- **AppSDL2OGL3.h** — Base application class: SDL2 init, event loop, keyboard/mouse camera controls, multi-screen support
- **ScreenSDL2OGL3.h** / **.cpp** — SDL window + OpenGL context manager, camera holder, scene list iteration
- **SceneOGL3.h** — Scene base class with `draw(Camera&)` virtual, `setCameraPersp`/`setCameraOrtho`/`setCamera` uniform helpers
- **SceneNode.cpp** — Scene graph node implementation
- **Shader.h** / **.cpp** — GLSL shader compilation, linking, uniform setting (`setUniformMat4f`, `setUniformVec3f`, etc.), `GLSL()` macro for inline source
- **GLfunctions.h** — Low-level GL helpers: buffer creation/upload, texture helpers, `checkOpenGLError()`
- **GLobjects.h** — `GLMesh`: VBO-based mesh (pos/nor/col/UV + index), `GL_STATIC_DRAW` default, no VAO caching
- **GLObject.h** / **.cpp** — `GLObject`: VBO/VAO management with up to 4 attribute buffers, `GLBuff` single-buffer abstraction
- **GLInstances.h** — `GLInstances` (instanced meshes), `GLBillboards` (billboarded quads), `GLParticles` (particle system) — all use `glVertexAttribDivisor`, no VAOs
- **DrawOGL3.h** — `GLMeshBuilder`: CPU-side mesh builder with `addLine`/`addPoint`/transforms, primitive generators (`Sphere2Mesh`, `Cone2Mesh`, `Torus2Mesh`)
- **MeshRenderOGL3.h** / **.cpp** — VAO-based mesh renderer with solid (lit) + line (const-color) shaders, wraps `GLMesh` with proper VAO caching
- **TextRendererOGL3.h** — Instanced billboarded 3D text: VAO + static quad + streaming instance VBO, binary alpha discard shader
- **SDFFont.h** — SDF (Signed Distance Field) font rendering: TTF atlas generation via SDL2_ttf + Felzenszwalb-Huttenlocher distance transform, instanced billboard renderer with screen-space adaptive `smoothstep` edges (crisp at any zoom)
- **DebugView.h** — Global debug rendering: `GLMeshBuilder` + `Shader`, `DEBUG_draw()` (inefficient: allocates/frees GLMesh every frame)
- **DrawField.h** — Scalar field visualization
- **DrawUV.h** — UV-mapped mesh generation from parametric functions
- **DrawVegetation.h** — Vegetation rendering
- **HorizontOGL3.h** — Horizon/sky rendering
- **TerrainOGL3.h** — Terrain mesh rendering
- **GL3Utils.h** — Miscellaneous GL3 utility functions
- **Michal/** — Alternative implementations by Michal: `Draw3D_o3.cpp`, `GLMesh.h` (templated), `MeshLibrary.h`, `Framebuffer.h`
