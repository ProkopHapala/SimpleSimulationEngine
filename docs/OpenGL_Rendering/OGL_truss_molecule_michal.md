# Truss & Molecule Rendering in Michal SDL2OGL3/Michal System

_Working notes on how to render mass–points / atoms, bonds, and labels using Michal’s abstractions (`GLES`, `GLMesh`, `GLInstancedMesh`, `MeshLibrary`, `TextRenderer`)._

## 1. Problem statement

We want the same capabilities as in the Prokop document, but implemented on top of Michal’s engine:

- **Atoms / nodes** as **sphere impostors** in the fragment shader (billboarded quads, `atom.glslf`-style shading).
- **Bonds / truss bars** as **hardware lines** (`GL_LINES`) between those atoms.
- **Labels** as **billboarded text** at each atom, using a font atlas and per-letter quads.

Additionally, we want the same architectural constraint:

- Atom positions live **once on the GPU** (single source of truth).
- Atoms, bonds, and labels all read from that same GPU buffer (texture or VBO/SSBO), so future GPU-based simulation can update positions without CPU round-trips.

Michal’s system provides higher-level building blocks for this:

- `GLES`: central GL state (active camera, framebuffer stack, screen size).
- `GLMesh` / `GLInstancedMeshBase`: typed meshes with built-in shaders.
- `MeshLibrary`: reusable primitives and a sphere impostor shader (`sphere`, `sphereInstanced`).
- `TextRenderer` + `Draw3D_o3`: text rendering and 3D billboards.

---

## 2. Shared GPU buffer for atom positions

We first decide how atom positions are stored on the GPU, since **everything else will read from that**.

### 2.1. Options for position storage

There are two main options, analogous to the Prokop document:

- **Option A – Position texture**
  - Create a `GLTexture` that stores positions in a float texture (e.g. RGBA32F or RGB32F):
    - You can follow the JS pattern (`uPosTex`, `uTexSize`) and update it from CPU with a `std::vector<Vec4f>`.
    - Shaders for atoms, bonds, and labels take an integer atom ID and fetch `pos = texelFetch(uPosTex, id, 0).xyz` (or equivalent 2D addressing).
  - Pros: matches JS `MeshRenderer` design closely.
  - Cons: slightly more plumbing in GLSL for texture coordinates; WebGL-style limitations are less relevant here.

- **Option B – Shared VBO / SSBO**
  - Use a single `GLvbo<MPOS>` (or an SSBO) that stores `vec3` positions for each atom.
  - Atoms:
    - Use this VBO as the per-instance attribute source for an instanced billboard/sphere mesh.
  - Bonds:
    - Use this VBO as the vertex attribute source for `GL_LINES`, with an index buffer of atom indices.
  - Labels:
    - Use this VBO as per-instance anchor positions for per-letter quads.
  - Pros: very natural in modern GL, integrates well with Michal’s `GLvbo`/`GLMesh`.

For Michal’s system, **Option B (shared VBO)** is usually more idiomatic, but the document keeps Option A in mind for easier comparison with JS.

### 2.2. Shared positions VBO

Assuming **Option B**:

- Create `GLvbo<MPOS>` (or some custom attribute structure) to hold all atom positions.
- When simulation updates positions on CPU:
  - Update the VBO once (e.g. `positionsVbo.verts[i] = newPos[i]`) and call `positionsVbo.sync()`.
- All rendering passes treat this VBO as the canonical source of `vec3` positions.

If you later move simulation to the GPU, the same VBO/SSBO can be updated by compute shaders or transform feedback without touching the CPU.

---

## 3. Atoms as sphere impostors

Instead of writing our own impostor shader from scratch, we can reuse or mirror Michal’s existing sphere impostor machinery:

- `MeshLibrary::sphere` and `MeshLibrary::sphereInstanced` are built on the generated shaders from `makeSphereVertexShader<instanced>()` and `makeSphereFragmentShader<instanced>()`.
- They already implement a billboard-like quad that is shaded as a sphere (view-dependent impostor) and compute `uMVPMatrix` / `uMVPinv` for correct depth and lighting.

### 3.1. Using `MeshLibrary::sphereInstanced`

`MeshLibrary` provides:

- `MeshLibrary::sphereInstanced`:
  - Type: `GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET, MRADIUS, MCOLOR>`.
  - Per-instance attributes:
    - `vPosOffset` (position offset / center) – `MPOSOFFSET`.
    - `vRadius` – `MRADIUS`.
    - `vColor` – `MCOLOR`.
  - Uniforms:
    - `uMVPMatrix`, `uMVPinv` (provided from `GLES::active_camera`).

To use it for atoms:

1. **Instantiate and configure**
   - Ensure `MeshLibrary::sphereInstanced` is created (it is in `MeshLibrary.cpp`).
   - Set global uniforms once per frame:
     - `sphereInstanced.uniforms.set4m("uMVPMatrix", GLES::active_camera->viewProjectionMatrix());`
     - `sphereInstanced.uniforms.set4m("uMVPinv",   GLES::active_camera->inverseViewProjectionMatrix());`

2. **Feed per-atom data**
   - Clear instances each frame: `MeshLibrary::sphereInstanced.instances->clear();`
   - For each atom `i`:
     - Add instance: `MeshLibrary::sphereInstanced.addInstance( center, radius, color );`
       - `center` is the atom position.
       - `radius` is atomic radius / visualization scale.
       - `color` is element/type color.

3. **Draw**
   - Simply call `MeshLibrary::sphereInstanced.draw();`
   - This issues a single instanced draw call for all atoms, using the impostor shader.

This gives you a **GPU-efficient atom impostor renderer** without having to write your own GLSL or attribute setup.

### 3.2. Sharing positions

If you want a **single shared positions buffer**:

- You can either:
  - Keep `sphereInstanced`’s own instance buffer as the canonical positions store (atoms, bonds, and labels all use the same instance indices), or
  - Keep a separate positions VBO/SSBO and copy positions into `sphereInstanced.instances` each frame (slightly redundant, but can still be driven from a GPU buffer later).

For strict "positions live only once" semantics, you’d adapt `GLInstancedMeshBase` to point at your shared positions buffer, but that is an implementation detail—the conceptual design remains: **instanced impostor spheres read from a shared positions source**.

---

## 4. Bonds / truss bars as GL_LINES

We want bonds drawn as simple hardware lines (`GL_LINES`), with endpoints derived from the **same positions buffer** used by atoms.

### 4.1. GPU positions + index buffer

Given a shared positions VBO (or texture/SSBO):

- Maintain an array of bond index pairs `(i, j)` on the CPU.
- Upload these to an index buffer bound to a `GLMesh<MPOS>` that uses the shared positions as its vertex attribute:
  - Positions come from `GLvbo<MPOS>`.
  - Indices `(i, j)` define which atoms are connected.

With Michal’s `GLMesh`:

1. **Create bond mesh**
   - `GLMesh<MPOS> bonds(GL_LINES);`
   - Instead of storing its own `verts`, you can:
     - Point its `verts` to the **same `GLvbo<MPOS>` as atoms** (shared positions), and
     - Give it an index buffer of bond pairs.
   - If the current `GLMesh` API does not directly expose shared VBOs + custom indices, an alternative is to create a dedicated `GLMesh<MPOS>` for bonds that is filled from the shared positions VBO (still no extra CPU copy if you update from GPU later).

2. **Vertex shader**
   - If using a shared VBO + index buffer, the default generated shader for `GLMesh<MPOS>` is enough:
     - `position` attribute comes from the shared positions VBO.
     - `glDrawElements(GL_LINES, ...)` with the bonds index buffer ensures each vertex uses the correct atom position.
   - If using a position texture (Option A):
     - Write a small custom shader that mirrors JS `bond.glslv`:
       - Take an `aAtomID` attribute (or index from the element buffer).
       - Fetch `pos` from `uPosTex` using `uTexSize` and `aAtomID`.

3. **Draw bonds**
   - Set a uniform color (`uColor`) or per-vertex colors if you need stress/strain visualization.
   - Call `bonds.draw();` (or equivalent `draw(GL_LINES)`).

The key property remains: **the positions are not re-uploaded or duplicated per bond**; bonds merely refer to atom indices.

---

## 5. Labels using TextRenderer & Draw3D_o3

Michal’s system already includes a small text subsystem:

- `TextRenderer`:
  - Uses a bitmap font atlas (`dejvu_sans_mono_RGBA_pix-UpDown.bmp`).
  - Builds a `GLMesh<MPOS,MUV>` of quads per character using `MeshBuilder::addChar`.
  - Can draw in 2D (`draw2D`) or 3D (`draw3D`).
- `Draw_o3` + `Draw3D_o3`:
  - Provide convenience wrappers for drawing text in world space but billboarded to the camera (`drawTextBillboard`).

### 5.1. Basic world-space labels

The simplest integration—without full per-letter instancing—is:

- For each atom, call `Draw3D::drawTextBillboard(label, pos, size, ontop)`.
- Internally this:
  - Projects `pos` via `GLES::active_camera->world2Screen(pos)`.
  - Uses `Draw::drawText` (which wraps `TextRenderer`) to draw a small label in screen space at that location.

This is already a **working, billboarded label system** that uses a font atlas, but does **not yet share the same position buffer** in a strongly coupled way (it takes positions as input arguments).

### 5.2. Per-letter instanced labels (analogous to JS `label.glslv`)

If you want the same level of control as in the Prokop + JS design (per-letter quads, labels driven from a shared positions buffer), you can:

1. **Use a shared positions buffer** (texture or VBO) as described in section 2.

2. **Design a label instanced mesh**
   - Use `GLInstancedMeshBase` with:
     - Per-vertex quad (for a unit glyph).
     - Per-instance attributes:
       - `atomID` (or direct `Vec3f` position if you don’t want to fetch from texture).
       - `charIndex` (glyph index into the font atlas).
       - `charPos` / `strLen` (position within the label and label length, for centering).
   - The vertex shader:
     - Uses `atomID` to fetch the anchor position from the shared positions buffer (or reads from a per-instance `Vec3f` if you choose that).
     - Computes per-character offset along the label baseline (`charPos`, `strLen`).
     - Turns `charIndex` into `vUv` via atlas grid (`uFontGrid`) as in JS `label.glslv`.

3. **Fragment shader**
   - Very similar to `label.glslf`:
     - Samples the font atlas texture.
     - Discards low alpha.
     - Outputs `vec4(uColor, alpha)`.

4. **Feeding label data**
   - On CPU, keep labels per atom (up to N characters) as `std::string`.
   - Build a per-letter instance buffer describing:
     - Which atom (`atomID`).
     - Which glyph (atlas index).
     - Character position within the string.
   - Upload as instanced attributes and draw once per frame.

This is more work than using `Draw3D::drawTextBillboard`, but:

- It matches the JS `label.glslv/label.glslf` design.
- It can be driven entirely by GPU-side positions if you choose Option A (position texture).

If you are okay with labels being slightly decoupled from the shared positions buffer (taking `Vec3f` positions per atom from CPU), `Draw3D::drawTextBillboard` may be sufficient and much simpler.

---

## 6. Mapping from JS `MeshRenderer` concepts to Michal system

JS `MeshRenderer` structure recap:

- **Position storage**: float texture `uPosTex` updated by `updatePositions`.
- **Atoms**: instanced plane geometry + `atom` shader pair, fetching position from `uPosTex` using `aAtomID`.
- **Bonds**: `LineSegments` whose vertex shader fetches positions from `uPosTex` (`bond.glslv`).
- **Labels**: instanced quads using `uPosTex` + font atlas (`label.glslv/label.glslf`).

In Michal’s C++ system, a consistent mapping is:

- **Shared positions buffer**:
  - Either:
    - A position texture owned by `GLTexture` (JS-like), or
    - A `GLvbo<MPOS>` / SSBO storing `Vec3f` positions.
- **Atoms**:
  - Reuse `MeshLibrary::sphereInstanced` as the atom impostor renderer.
  - Feed instance data (`pos`, `radius`, `color`) from the shared positions buffer.
- **Bonds**:
  - A `GLMesh<MPOS>` with draw mode `GL_LINES` that points to the shared positions buffer and uses an index buffer of atom indices, or
  - A custom shader (if using position texture) that mirrors `bond.glslv`.
- **Labels**:
  - Simple version: `Draw3D::drawTextBillboard` with `TextRenderer` and camera projection.
  - Advanced version: a `GLInstancedMeshBase` of per-letter quads that fetch anchor positions from the shared positions buffer and atlas UVs from a font grid (JS-like).

The core idea is preserved: **one GPU representation of atom positions, multiple high-level passes (spheres, lines, labels) built on top using Michal’s typed mesh and shader abstractions.**

---

## 7. Recommended Michal pipeline summary

Putting it together for Michal’s system:

- **Atom spheres**
  - Use `MeshLibrary::sphereInstanced`.
  - Maintain a shared positions buffer (texture or VBO); fill instances per atom with center, radius, color.

- **Bonds / truss bars**
  - Use a `GLMesh<MPOS>` with draw mode `GL_LINES`.
  - Positions sourced from the shared positions buffer.
  - Index buffer encodes which atoms are connected.

- **Labels**
  - Start with `Draw3D::drawTextBillboard` + `TextRenderer` for practical, debug-friendly labels.
  - Optionally move to per-letter instanced labels using `GLInstancedMeshBase` if you want closer parity with the JS system and full GPU-driven positioning.

This document focuses on how to **compose Michal’s higher-level building blocks** to solve the same truss/molecule visualization problem, while keeping the "single GPU source of atom coordinates" goal in mind for future GPU-side simulation.
