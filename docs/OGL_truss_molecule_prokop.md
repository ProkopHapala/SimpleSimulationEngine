# Truss & Molecule Rendering in Prokop SDL2OGL3 System

_Working notes on how to render mass–points / atoms and bonds with raycasted spheres and labels using the original SDL2OGL3 engine._

**Design note (MeshRenderOGL3 data binding)**

The planned `MeshRenderOGL3` helper in `cpp/common_SDL/SDL2OGL3` must **not** hard‑code any particular mesh, truss, or molecule geometry inside the renderer itself. Instead, it will:

- Accept **external raw C arrays** (pointers) and counts for:
  - Vertex / atom positions.
  - Bonds / truss bar connectivity (index pairs).
  - Optional per‑atom labels / glyph indices.
- Or accept data via **explicit setup functions** that upload these user‑provided arrays into VBOs / textures.

This way, the same renderer can be reused both:

- As a **standalone app** (e.g. `meshviewer.cpp` in `cpp/sketches_SDL/OGL3`).
- As a **library component** (e.g. from `cpp/libs_SDL/GLView` or `cpp/libs_SDL/Vis3D`).

All simulation / domain code remains responsible for owning and updating the arrays; `MeshRenderOGL3` is a thin visualization layer that only knows how to bind these arrays to GPU buffers and draw them.

## 1. Problem statement

We want a reusable pattern in the **Prokop** C++ engine to render:

- **Points / atoms / nodes** as **raycasted spheres** in the fragment shader (screen‑space impostors or analytic ray–sphere intersection).
- **Bonds / truss bars** as **line segments or cylinders** between those points.
- **Labels** (indices, element types, chemical symbols) as **billboarded text** at or near each point.

Conceptually this mirrors the JS `MeshRenderer` pattern:

- A storage of per‑atom positions (in JS: a float texture `uPosTex`).
- An instanced quad for atoms, shaded as spheres (`atom.glslv/atom.glslf` or `sphere_frag.c` style raycasting).
- A line segment mesh for bonds, also driven by the same positions (`bond.glslv`).
- An instanced label mesh that samples a font atlas and draws a short string per atom (`label.glslv/label.glslf`).

In C++ SDL2OGL3 we don’t have the exact same texture‑buffer setup, but we do have:

- Flexible VBOs/VAOs (`GLObject`, `GLMesh`, helpers in `GL3Utils.h`).
- Instancing helpers (`GLInstances`, `GLBillboards`, `GLParticles` in `GLInstances.h`).
- Working examples of sphere impostors (`test_SphereShader.cpp` + `shaders/sphere_frag.c`).
- Working sprite/billboard shaders (`hardSprite.glslf`, `Bilboard3D.glslv`, texture/atlas shaders).

## 2. Current OGL3 Text Rendering (Instanced Characters)

For Modern OpenGL 3 text / label rendering we use an **instanced-character** design (see `TextRendererOGL3`):

- **Geometry:** One rectangular quad in local `[0,1]×[0,1]` (two triangles).
- **Instances:** One instance per character (`GlyphInstance` with `pos`, `charIndex`, `ascii`).
- **Texture atlas:** A logically 1D font atlas in a 2D texture (one row of glyphs). UV mapping follows the old `Draw::drawText` logic, but implemented in a GLSL 330 shader (no fixed pipeline).
- **Shader:**
  - VS receives quad vertex (`aQuadPos`) and per-instance attributes (`aBasePos`, `aCharIndex`, `aAscii`).
  - Computes `worldPos` as a camera-aligned billboard using `uCamRight`/`uCamUp` and character index/size.
  - Computes atlas `u` as `(ascii - glyphOffset) * uAtlasStep + local_u * uAtlasStep`.
  - Outputs `vUV` and transforms by `uVP` (camera matrix built the same way as in `MeshRenderOGL3::draw`).
  - FS samples `uFontTex` and multiplies by `uTextColor`, discarding low-alpha pixels.

### 2.1. TextRendererOGL3 implementation details

- Lives in `cpp/common_SDL/SDL2OGL3/TextRendererOGL3.h`.
- Uses **VBO/EBO only** (no VAOs), mirroring `GLMesh` style:
  - `vboQuad` – static quad vertices.
  - `vboInst` – dynamic per-glyph instances.
  - `ebo` – indices `{0,1,2, 0,2,3}`.
- In `draw()` it explicitly sets vertex attributes every frame:
  - Bind `vboQuad`, enable attrib 0 (`aQuadPos`).
  - Bind `vboInst`, enable attribs 1–3 (`aBasePos`, `aCharIndex`, `aAscii`) and set divisors.
  - Bind `ebo` and call `glDrawElementsInstanced`.
  - Then disable attribs 0–3 and unbind buffers again.

This avoids interfering with whatever VAO / attribute state is implicitly used by `GLMesh` and other legacy OGL3 helpers.

### 2.2. Why VAOs are currently avoided (and future plan)

During initial experiments a separate VAO created in `TextRendererOGL3::init()`:

```cpp
glGenVertexArrays(1, &vao);
glBindVertexArray(vao);
// ... setup attribs ...
glBindVertexArray(0);
```

caused the main mesh to stop rendering or render incorrectly. The rest of the engine (especially `GLMesh`) does **not** use per-object VAOs; it binds buffers and calls `glVertexAttribPointer` directly on the default / global VAO.

Therefore, for now:

- `TextRendererOGL3` **does not own a VAO**; it only owns buffers.
- It configures attributes on whatever VAO is current, and then cleans them up.

**Future optimization note (VAO sharing):**

- For performance and cleaner state management it is desirable to reintroduce VAOs later.
- When that happens, the VAO must be **shared or at least restored** consistently:
  - Either: a global "OGL3 default VAO" that both `GLMesh` and `TextRendererOGL3` use and configure.
  - Or: each renderer binds its own VAO at the start of `draw()` and restores the previous VAO afterwards.
- The key lesson from this debugging session is: creating a new VAO and leaving it bound (or assuming others don’t rely on the default VAO) can silently break mesh rendering.

Until that refactor, the safe pattern is:

- Use VBO/EBO + explicit `glVertexAttribPointer` per draw (like `GLMesh::preDraw()` / `postDraw()`).
- Avoid `glGenVertexArrays` in new OGL3 helpers unless VAO ownership and sharing are carefully designed.

This document sketches how to build the same functionality in your C++ system.

---

## 2. Rendering atoms as raycasted spheres

There are two closely related approaches in your codebase:

- **Billboarded impostor sphere (preferred)**: instanced quads with a fragment shader like `atom.glslf` or `pointSprite.glslf`.
- **Screen‑space raycasted sphere (reference)**: `test_SphereShader.cpp` + `shaders/sphere_frag.c`.

### 2.1. Reference: `test_SphereShader.cpp` + `sphere_frag.c`

Files:

- `cpp/sketches_SDL/OGL3/test_SphereShader.cpp`
- `cpp/sketches_SDL/OGL3/shaders/sphere_frag.c`

Key ideas:

- The **geometry** is just a screen‑aligned quad (`GLObject` with 4 vertices, `GL_TRIANGLE_STRIP`).
- The **fragment shader** computes, per pixel, a ray from the camera through the pixel, intersects it with a sphere, and shades it:
  - Uniforms:
    - `resolution` – viewport size.
    - `sphere` – `vec4(center.xyz, radius)`.
    - `light_dir` – light direction.
  - Functions `rayPointDist`, `raySphere`, `sphereNormal` perform the analytic intersection and normal reconstruction.
  - `gl_FragDepth` is written so the sphere integrates with the depth buffer.

This example is useful mainly as a **mathematical reference** for ray–sphere intersection and analytic normals. For large atom counts, it’s not the right rendering pattern; we should instead use **instanced billboards** as in the JS `atom.glslv/atom.glslf` shaders.

### 2.2. Instanced billboard spheres (impostors, using `atom.glslf`)

This is directly analogous to the JS `atom.glslv/atom.glslf` pair and is the **recommended** way to draw many atoms:

- Per‑instance data:
  - Atom center position `pos[i] : vec3`.
  - Atom radius / scale `radius[i]`.
  - Color `color[i] : vec3`.
- Per‑vertex data (for the quad):
  - `position.xy` in `[-0.5, 0.5]^2`.

In your C++ engine, the natural building block is **`GLBillboards`** from `GLInstances.h`:

- `GLBillboards` already manages:
  - A quad mesh (UVs and positions for a unit billboard).
  - An instance buffer with per‑instance position and scale.
  - Instanced draw using `glDrawArraysInstanced`.
- The shader side (`Bilboard3D.glslv` family) knows how to:
  - Take camera rotation and per‑instance `pose_pos`, `pose_sc`.
  - Generate world‑space vertices that always face the camera.

To turn billboards into **spheres**:

- Reuse the billboard-style vertex shader (to get quads in front of the camera for each atom), patterned after `atom.glslv`:
  - Per‑instance `pos` and `radius` determine quad placement and size.
  - Per‑instance color is passed through to the fragment shader.
- Use a **sphere impostor** fragment shader modeled on `atom.glslf`:
  - Use interpolated `vUv` (or equivalent) to get local coordinates in the quad.
  - Discard pixels outside the circle.
  - Reconstruct the virtual normal on the sphere (`z = sqrt(r^2 - x^2 - y^2)`).
  - Do simple Phong/Blinn lighting, as in `atom.glslf`.

Implementation steps in Prokop system:

1. **Define atom instance data**
   - CPU arrays: `atom_pos[i]`, `atom_radius[i]`, `atom_color[i]`.
   - Choose layout consistent with `GLBillboards` (position+scale).

2. **Initialize `GLBillboards`**
   - Similar to `test_Sprites.cpp`:
     - Set `nInstances` to atom count.
     - Upload per‑instance positions and radii to the instance VBO.

3. **Create / load shaders**
   - Vertex shader: use your existing billboard vertex shader (`Bilboard3D.glslv` variant) or a simplified one that:
     - Takes per‑instance `pos` and `scale`.
     - Outputs `vUv` in `[0,1]^2` for the quad.
   - Fragment shader: port `atom.glslf` or `pointSprite.glslf` logic to GLSL 3.3 / your style.

4. **Draw atoms**
   - Bind the atom billboard shader.
   - Set camera uniforms (`camMat`, `camRot`, etc.).
   - Bind textures if needed (e.g. material LUT, environment).
   - Call `GLBillboards.draw(GL_TRIANGLES)`.

This yields an **efficient, instanced atom renderer** directly analogous to the JS `InstancedMesh` + `atom.glslv/atom.glslf` approach, without needing a fullscreen `sphere_frag.c` pass.

---

## 3. Rendering bonds / truss bars

JS reference: `bond.glslv` takes atom IDs (`aAtomID`), reads positions from a **single position texture** `uPosTex`, and emits line segment vertices. This means **atoms, bonds, and labels all share the same GPU position buffer**.

For the Prokop system we want the same property:

- Atom positions should live **once** on the GPU (VBO or position texture / TBO / SSBO).
- Atoms, bonds, and labels should all read from that same buffer.
- Bonds are rendered as `GL_LINES` using **indices of atoms**, not a separate CPU‑side copy of positions.

### 3.1. Single GPU source of truth for positions

On the CPU you still maintain a logical array `pos[i] : Vec3d` for simulation, but for rendering you upload it to **one GPU buffer**:

- Option A (more like JS):
  - Pack positions into a float texture (e.g. `GL_TEXTURE_2D` or `GL_TEXTURE_BUFFER`) and index it in shaders with an integer atom ID, just like `uPosTex` and `aAtomID`.
- Option B (more classic GL):
  - Use a VBO with `vec3 position` per atom and use that VBO as:
    - A per‑instance attribute source for atoms.
    - A vertex attribute source for bonds, accessed via an **index buffer of atom indices** (`GL_ELEMENT_ARRAY_BUFFER`).

Both options satisfy the requirement: **one GPU buffer/texture is the source of truth for atom coordinates**, shared by atom quads, bond lines, and label billboards.

### 3.2. Bonds as indexed `GL_LINES`

Given per‑atom positions on the GPU, bonds/truss bars are defined by a list of **pairs of atom indices** `(i, j)`:

- You create an index buffer `bondsIBO` storing `[i0, j0, i1, j1, ...]` as unsigned ints.
- You render with draw mode `GL_LINES`:
  - The vertex shader receives the atom position for index `i`/`j` from the **shared positions buffer**.
  - No per‑bond position duplication is needed; only indices are stored.

In a JS‑style **position texture** approach (Option A):

- The C++ vertex shader for bonds is almost identical to `bond.glslv`:
  - An attribute or element index provides `aAtomID`.
  - The shader computes texture coordinates from `aAtomID` and fetches the atom position from the shared `uPosTex`.
- The index buffer still holds `(i, j)` pairs, but the shader turns them into positions by sampling the texture.

In a **shared VBO + index buffer** approach (Option B):

- Positions are in a `GL_ARRAY_BUFFER` bound as vertex attribute `position`.
- Bonds are drawn with `glDrawElements(GL_LINES, ...)` using the bonds index buffer.
- Atoms and labels can still use the same VBO as a per‑instance attribute source.

Either way, the important points are:

- Bonds are drawn as **simple hardware lines** (`GL_LINES`).
- All three passes (atoms, bonds, labels) read positions from the **same GPU buffer**.
- The CPU only uploads positions once per frame (or once per simulation step) and never round‑trips atom coordinates between CPU and GPU for different passes.

---

## 4. Rendering labels at atom positions

We want labels like in the JS `label.glsl*` pair:

- Each atom has up to 8 characters.
- Labels rendered as small quads, either in screen space or in world space but billboarded.

In your C++ system, you already have all necessary primitives:

- **Billboarded quads** (`GLBillboards` again).
- **Texture atlas shaders** (`Bilboard3D.glslv`, texture shaders, atlas variants).
- **Text rendering in JS** (and now, in Michal’s C++ `TextRenderer`), which you could conceptually mirror.

For a pure Prokop implementation, a straightforward approach is:

1. **Choose a font atlas**
   - Reuse the existing bitmap font texture `js/common_resources/dejvu_sans_mono_RGBA_inv.bmp`, where characters are arranged in a regular grid (compatible with a `uFontGrid` uniform as in the JS shaders).
   - Load this texture using your existing texture helpers (`newTexture2D` etc.).

2. **Per‑character billboards (no pre‑rendering)**
   - We do **not** pre‑render whole labels to textures; instead each **letter is one sprite (one billboard)**.
   - For each atom and each character in its label:
     - Create an instance with:
       - World position (same as atom position, plus a small offset).
       - Per‑instance data encoding which glyph to show (character code / atlas index) and where along the label baseline it belongs.
   - Use a billboard vertex shader similar to `label.glslv` / your own `Bilboard3D.glslv`:
     - Start from the atom position as anchor.
     - Compute per‑character offset along a baseline (like `aCharPos` and `aStrLen` in `label.glslv`).
     - Compute `vUv` from `charIndex` and `uFontGrid` so that the correct cell of `dejvu_sans_mono_RGBA_inv.bmp` is sampled.
     - Optionally support `uScreenSpace` vs world‑space scaling, as in `label.glslv`.

3. **Fragment shader**
   - Sample the font atlas texture.
   - Use alpha thresholding like `label.glslf` and the old `SDL2OGL/Draw.cpp` code:
     - Discard if alpha is low.
     - Output `vec4(uColor, alpha)`.

4. **Data feeding**
   - On CPU, store labels per atom (`std::string` or fixed `[8]` buffer).
   - Build instance buffers for characters (similar to JS `aAtomID`, `aLabel1`, `aLabel2`, `aCharPos`, `aStrLen`), but mapping to your own attribute layout in `GLInstances`/`GLBillboards`.
   - Compute the glyph indices exactly as you already do in `SDL2OGL/Draw.cpp` (where characters are mapped to atlas columns by subtracting a base ASCII code and multiplying by `persprite`).
   - Upload these attributes to your instancing VBO and draw with `glDrawArraysInstanced`.

Given your preference for performance and control, you might initially implement a **simpler variant** (e.g. fixed‑length labels, 1–2 characters per atom) and extend the per‑character instancing scheme later once the atlas mapping is in place.

---

## 5. Mapping from JS `MeshRenderer` concepts

JS `MeshRenderer` structure:

- **Position storage**: float texture `uPosTex` filled by `updatePositions`.
- **Atoms**: instanced quad mesh with per‑instance `instanceColor` and `instanceScale`, vertex shader fetches atom position from `uPosTex` using `aAtomID`.
- **Bonds**: line segments with per‑vertex `aAtomID` and vertex shader fetching positions.
- **Labels**: instanced mesh using `uPosTex` + font atlas, attributes describing up to 8 characters.

In your C++ engine:

- You can still keep **positions in CPU arrays** for simulation, but for rendering they should be uploaded to **one GPU positions buffer** and then reused everywhere.

One sensible mapping mirroring the JS design is:

- Keep the **single GPU source of positions**:
  - Either as a position texture (JS style) or as a VBO that plays the same role.
- Use integer atom IDs / indices consistently:
  - Atoms: per‑instance ID or direct indexing into the positions buffer in the vertex shader.
  - Bonds: pairs of atom indices in an index buffer, with the bond vertex shader fetching positions by ID (texture) or letting the index select from the shared VBO.
  - Labels: per‑letter instances that carry an atom ID to locate the anchor position, fetching from the same positions buffer.
- Keep the **shader ideas** similar to JS:
  - Atom impostors: `atom.glslv/atom.glslf` logic.
  - Bonds: `bond.glslv` logic (adapted to your positions buffer representation).
  - Labels: `label.glslv/label.glslf` glyph indexing and positioning, plus your atlas UV mapping.

This preserves your preference for control while ensuring the critical property you want: **GPU‑only atom positions shared across atoms, bonds, and labels**, enabling fully GPU‑side simulation of positions with no redundant CPU↔GPU transfers.

---

## 6. Recommended Prokop pipeline summary

A practical implementation in the Prokop SDL2OGL3 system could look like:

- **Atoms**
  - Use `GLBillboards` with a custom `AtomBillboard` shader (vertex similar to `Bilboard3D.glslv`, fragment similar to `atom.glslf`).
  - Per instance: `pos`, `radius`, `color`.

- **Bonds / truss bars**
  - Use a `GLMesh` built as lines between `pos[i]` and `pos[j]`.
  - Simple colored line shader, optionally with per‑bond color.

* **Labels**
  - Use the existing font atlas `dejvu_sans_mono_RGBA_inv.bmp` and a GLSL pair modeled on `label.glslv/label.glslf`.
  - Implement **per‑letter instanced billboards** (each glyph is one sprite) anchored at atom positions, with UVs computed from character codes exactly as in the OpenGL1 `Draw::drawText` implementation, but using your modern instancing infrastructure. No pre‑rendering of whole labels and no external text system is required.

This document focuses on **what** to do and how your existing modules fit together; the exact GLSL and C++ code can then be derived from the referenced examples (`test_SphereShader.cpp`, `sphere_frag.c`, `GLInstances.h`, sprite/atlas shaders, and the JS shaders).
