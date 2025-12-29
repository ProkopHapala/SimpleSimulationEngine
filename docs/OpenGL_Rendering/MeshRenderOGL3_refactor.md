## User

I would like to ask you about best strategy how to work with VAO between mesh renderer and text renderer, notice that text labels are often used for vertexes of the mesh vetexes
We have some disputation about VAO here:

---

This is a classic architectural conflict in OpenGL engine design: **State Management vs. Encapsulation**.

Since you are porting logic where "Mesh" and "Text" are separate entities but share semantic data (positions), and you are currently avoiding VAOs to maintain compatibility with legacy `GLMesh` code, here is the definitive strategy.

### 1. The Verdict on Sharing VBOs (Mesh Vertices vs Text Positions)

**Do not try to use the Mesh VBO directly in the Text Renderer.**

You asked: *"Labels are often used for vertexes of the mesh... can we use that?"*

**Answer:** No, it is inefficient for the *Instanced Character* approach.
*   **Mesh VBO:** Contains $N$ vertices (one per atom).
*   **Text VBO:** Contains $M$ instances (one per *character*).
*   **The Mismatch:** If you have the label "Carbon" at vertex 10, you need 6 text instances. Each instance needs the position of vertex 10. To share the VBO, you would need to pass an `index` to the shader and fetch the position from a Texture Buffer Object (TBO). This adds shader complexity and setup overhead that outweighs the cost of simply copying the `Vec3` into the `GlyphInstance` vector on the CPU.

**Strategy:** Keep your `GlyphInstance` vector. It is the fastest, simplest way to bridge the gap between "Atoms" and "Characters".

---

### 2. The VAO Strategy: "Stateless" vs "Stateful"

You are currently seeing bugs because you are mixing **Stateless** code (legacy style, just binding VBOs) with **Stateful** objects (VAOs).

#### The Problem with your current Code
In `TextRendererOGL3::draw`, you set `glVertexAttribDivisor`. **This is persistent state** in the active VAO.
1. `MeshRenderOGL3` draws (binds default VAO, sets pointers).
2. `TextRendererOGL3` draws (binds default VAO, **sets Divisors to 1**).
3. Next frame, `MeshRenderOGL3` draws again. **The Divisors are still set to 1** on attributes 1, 2, and 3. The mesh disappears or explodes because the GPU thinks the mesh vertices are instanced data.

#### Solution A: The "Clean Neighbor" Policy (Recommended for NOW)
If you want to avoid creating VAOs in `TextRenderer` to stay compatible with `GLMesh`, you **must** strictly clean up the Attribute Divisors.

Update your `draw` function in `TextRendererOGL3.h`:

```cpp
    void draw(const Camera& cam){
        if(!fontTex || instances.empty()) return;

        // ... [Camera setup code] ...

        shader.use();
        // ... [Uniform setup code] ...

        // 1. SAVE STATE (If we were using a specific VAO, we would bind it here)
        // Since we share the global state, we just bind our buffers.
        
        // --- Setup Geometry ---
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0);

        // --- Setup Instances ---
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        
        // Attrib 1: BasePos
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance,pos));
        glVertexAttribDivisor(1, 1); // <--- We dirty the state here

        // Attrib 2: CharIndex (Integer!)
        glEnableVertexAttribArray(2);
        glVertexAttribIPointer(2, 1, GL_INT, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance,charIndex));
        glVertexAttribDivisor(2, 1);

        // Attrib 3: Ascii (Integer!)
        glEnableVertexAttribArray(3);
        glVertexAttribIPointer(3, 1, GL_INT, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance,ascii));
        glVertexAttribDivisor(3, 1);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

        // --- DRAW ---
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, instances.size());

        // --- CLEANUP (CRITICAL) ---
        // We must reset divisors, otherwise GLMesh will try to use instancing next time it uses these locations!
        glVertexAttribDivisor(1, 0); 
        glVertexAttribDivisor(2, 0); 
        glVertexAttribDivisor(3, 0);

        glDisableVertexAttribArray(3);
        glDisableVertexAttribArray(2);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        
        glDisable(GL_BLEND);
    }
```

#### Solution B: The Encapsulated VAO (Future Goal)
Once you are ready to refactor `MeshRenderOGL3`, every renderer should own its own VAO.
1. `MeshRenderOGL3::init()` -> `glGenVertexArrays(1, &vao);`
2. `TextRendererOGL3::init()` -> `glGenVertexArrays(1, &vao);`
3. In `draw()`:
   ```cpp
   glBindVertexArray(this->vao);
   glDrawElements...
   glBindVertexArray(0); // Unbind to be safe for legacy code
   ```

---

### 3. Improved Implementation Notes

Here is a refined version of your `TextRendererOGL3.h` that fixes the integer attribute handling and the state pollution bugs.

**Key Changes:**
1.  **State Cleanup:** Added `glVertexAttribDivisor(..., 0)` at the end of draw.
2.  **Integer Attributes:** Confirmed use of `glVertexAttribIPointer` for `int` types (ASCII, CharIndex). This is correct in your code, but often forgotten.
3.  **Optimization:** Added a check to avoid re-uploading the instance buffer if the text hasn't changed (optional, but good practice).

```cpp
#ifndef TextRendererOGL3_h
#define TextRendererOGL3_h

#include "Vec3.h"
#include "Mat4.h"
#include "Camera.h"
#include "GLobjects.h"
#include "Shader.h"
#include <vector>

struct GlyphInstance{
    Vec3f pos;       
    int   charIndex; 
    int   ascii;     
};

class TextRendererOGL3{
public:
    GLuint vboQuad  = 0;  
    GLuint vboInst  = 0;  
    GLuint ebo      = 0;  
    GLuint fontTex  = 0;

    Shader shader;

    float charW = 0.1f;   
    float charH = 0.1f;   

    int   numGlyphs    = 128;   
    int   glyphOffset  = 0; // Usually 0 or 32 depending on your atlas

    std::vector<GlyphInstance> instances;
    bool dirty = false; // Flag to check if we need to re-upload

    void init(GLuint fontTex_, float charW_, float charH_){
        fontTex = fontTex_;
        charW   = charW_;
        charH   = charH_;

        // Quad geometry
        const float quadVerts[8] = { 0.0f, 0.0f,  1.0f, 0.0f,  1.0f, 1.0f,  0.0f, 1.0f };
        const GLuint quadIdx[6]  = { 0,1,2, 0,2,3 };

        glGenBuffers(1, &vboQuad);
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);

        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIdx), quadIdx, GL_STATIC_DRAW);

        glGenBuffers(1, &vboInst); // Dynamic buffer
        
        // Shader Setup (Same as your provided code, omitted for brevity)
        // Ensure shader uses version 330 core
        initShader(); 
    }

    void initShader() {
         const char* vertSrc = GLSL(330,
            layout(location = 0) in vec2 aQuadPos;
            layout(location = 1) in vec3 aBasePos;
            layout(location = 2) in int  aCharIndex;
            layout(location = 3) in int  aAscii;

            uniform mat4 uVP;
            uniform vec3 uCamRight;
            uniform vec3 uCamUp;
            uniform float uCharW;
            uniform float uCharH;
            uniform float uAtlasStep;
            uniform int   uGlyphOffset;

            out vec2 vUV;

            void main(){
                // World Space Calculation
                float xOffset = float(aCharIndex) * uCharW;
                vec3 center = aBasePos + (uCamRight * xOffset);
                
                // Billboard expansion
                vec3 vertexPos = center 
                    + uCamRight * (aQuadPos.x - 0.5) * uCharW  // Centering adjustment if needed
                    + uCamUp    * (aQuadPos.y - 0.5) * uCharH;

                // Texture Space Calculation
                float glyphIndex = float(aAscii - uGlyphOffset);
                float uStart = glyphIndex * uAtlasStep;
                vUV = vec2(uStart + (aQuadPos.x * uAtlasStep), aQuadPos.y);

                gl_Position = uVP * vec4(vertexPos, 1.0);
            }
        );
        // ... frag shader ...
        // shader.init_str(...);
    }

    void clear(){ instances.clear(); dirty = true; }

    void addLabel(const Vec3f& basePos, const char* text){
        for(int i=0; text[i]; i++){
            GlyphInstance gi;
            gi.pos       = basePos;
            gi.charIndex = i;
            gi.ascii     = (unsigned char)text[i];
            instances.push_back(gi);
        }
        dirty = true;
    }

    void uploadInstances(){
        if(!dirty || instances.empty()) return;
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        // glBufferData orphan trick: pass NULL with size first to tell driver to allocate new block
        glBufferData(GL_ARRAY_BUFFER, instances.size()*sizeof(GlyphInstance), nullptr, GL_DYNAMIC_DRAW);
        glBufferData(GL_ARRAY_BUFFER, instances.size()*sizeof(GlyphInstance), instances.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        dirty = false;
    }

    void draw(const Camera& cam){
        if(!fontTex || instances.empty()) return;
        
        uploadInstances(); // Check dirty flag inside

        // Matrix setup ... 
        // shader.use(); ...

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // --- ATTRIBUTE BINDING ---
        // 0: Quad UV/Pos
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        // 1, 2, 3: Instance Data
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance,pos));
        glVertexAttribDivisor(1, 1);

        glEnableVertexAttribArray(2);
        glVertexAttribIPointer(2, 1, GL_INT, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance,charIndex));
        glVertexAttribDivisor(2, 1);

        glEnableVertexAttribArray(3);
        glVertexAttribIPointer(3, 1, GL_INT, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance,ascii));
        glVertexAttribDivisor(3, 1);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

        // Draw
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, instances.size());

        // --- CRITICAL CLEANUP ---
        glVertexAttribDivisor(1, 0);
        glVertexAttribDivisor(2, 0);
        glVertexAttribDivisor(3, 0);
        
        glDisableVertexAttribArray(0);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(2);
        glDisableVertexAttribArray(3);
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glDisable(GL_BLEND);
    }
};
#endif
```


---

## User 

the MeshRenderOGL3 is not mlegacy code, we just created that. Maybe or implementation is naieve / inefficiet, but we definitely can redesing it, just tell me what is the optimal implementation, performancewise (that is my main concern), but also from point ov view of ease of use, user convenience, elegance, flexibility. I'm not in favour of perfect encapsualtion just for sake of it (I o not follow object oriented paradigm), I rather follow data-oriented desing, and mutiparadigm philosophy.

---

You are absolutely right. If `MeshRenderOGL3` is new code, we should not treat it like legacy. We should design it for **high performance** and **low driver overhead**, following Data-Oriented principles.

Here is the **optimal design** for modern OpenGL (3.3 - 4.6), prioritizing performance and ease of use.

### The Philosophy: "Setup Once, Draw Cheaply"

The main performance bottleneck in your current code is **Driver Overhead**.
In your current `draw()` functions, you are calling `glVertexAttribPointer`, `glEnableVertexAttribArray`, and `glVertexAttribDivisor` every single frame.
*   **Why this is bad:** The CPU has to validate these pointers, check the VBO binding, and talk to the GPU driver every frame.
*   **The Solution:** Use **VAOs (Vertex Array Objects)** correctly. A VAO is essentially a "recording" of all those attribute setups.

### 1. The Strategy

1.  **TextRendererOGL3:** Owns 1 VAO. It records the format of the Quad VBO and the Instance VBO. In `draw()`, we just bind the VAO and draw.
2.  **MeshRenderOGL3:** Owns 1 VAO (or one per mesh type). It records the format of the position/normal/color buffers.
3.  **Data Sharing:** You asked about using Mesh vertices for Text.
    *   **Performance Reality:** It is actually *faster* for the GPU to read a linear "Instance Buffer" (Stream) than to perform "Dependent Texture Reads" (TBOs) to look up mesh positions, especially because text needs multiple instances (chars) per mesh vertex.
    *   **DOD Approach:** We will keep the `vector<GlyphInstance>` on the CPU. It is a flat, cache-friendly array. Copying positions from the mesh to this vector is incredibly fast (memory bandwidth speed) compared to the complexity of shader indirection.

---

### 2. Implementation: Fast Text Renderer (DOD Style)

This implementation uses **Buffer Mapping** (`glMapBufferRange`) which is the industry standard for updating dynamic geometry (particles/text) efficiently without allocating new memory constantly.

#### `TextRendererOGL3.h`

```cpp
#ifndef TextRendererOGL3_h
#define TextRendererOGL3_h

#include <vector>
#include <GL/glew.h>
#include "Vec3.h"
#include "Camera.h"
#include "Shader.h"

// Pure Data Struct
struct GlyphInstance {
    Vec3f pos;
    int   charIndex;
    int   ascii;
};

class TextRendererOGL3 {
public:
    // OpenGL State Handles
    GLuint vao      = 0;
    GLuint vboQuad  = 0;
    GLuint vboInst  = 0; // The dynamic buffer
    GLuint ebo      = 0;
    GLuint fontTex  = 0;

    Shader shader;

    // Config
    float charW = 0.1f;
    float charH = 0.1f;
    int   numGlyphs = 128;
    int   glyphOffset = 0;

    // CPU Data Store
    std::vector<GlyphInstance> instances;
    bool dirty = false;

    // --- 1. Initialization (The "Setup Once" phase) ---
    void init(GLuint fontTex_, float charW_, float charH_) {
        fontTex = fontTex_;
        charW = charW_;
        charH = charH_;

        initShader();

        // 1. Create the VAO first. All subsequent state is recorded into it.
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        // 2. Setup Static Geometry (The Quad)
        // [0,0] -> [1,1]
        const float quadVerts[] = { 0,0, 1,0, 1,1, 0,1 };
        const GLuint quadIdx[]  = { 0,1,2, 0,2,3 };

        glGenBuffers(1, &vboQuad);
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);

        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIdx), quadIdx, GL_STATIC_DRAW);

        // Attribute 0: Quad Position (2 floats)
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        // 3. Setup Dynamic Instance Buffer
        glGenBuffers(1, &vboInst);
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        // Pre-allocate memory (e.g., for 10,000 chars) to avoid resizing later.
        // GL_STREAM_DRAW indicates we modify it every frame but draw it few times.
        glBufferData(GL_ARRAY_BUFFER, 10000 * sizeof(GlyphInstance), nullptr, GL_STREAM_DRAW);

        // Attribute 1: World Position (3 floats, Instanced)
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance, pos));
        glVertexAttribDivisor(1, 1); // IMPORTANT: Update once per instance

        // Attribute 2: Char Index (1 int, Instanced)
        glEnableVertexAttribArray(2);
        glVertexAttribIPointer(2, 1, GL_INT, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance, charIndex));
        glVertexAttribDivisor(2, 1);

        // Attribute 3: ASCII Code (1 int, Instanced)
        glEnableVertexAttribArray(3);
        glVertexAttribIPointer(3, 1, GL_INT, sizeof(GlyphInstance), (void*)offsetof(GlyphInstance, ascii));
        glVertexAttribDivisor(3, 1);

        // 4. Unbind VAO to prevent accidental modification by other classes
        glBindVertexArray(0);
        // Unbind buffers (optional, but good for safety)
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    // --- 2. Data Manipulation (Data Oriented) ---
    void clear() { 
        instances.clear(); 
        dirty = true; 
    }

    void addLabel(const Vec3f& basePos, const char* text) {
        // Very fast CPU copy.
        // Modern CPUs can push vectors like this extremely fast.
        for (int i = 0; text[i]; i++) {
            // Emplace back constructs the object directly in the vector memory
            instances.push_back({ basePos, i, (int)text[i] });
        }
        dirty = true;
    }

    // --- 3. Rendering (The "Draw Cheaply" phase) ---
    void draw(const Camera& cam) {
        if (instances.empty() || !fontTex) return;

        shader.use();
        setupUniforms(cam); // Helper to set VP, CameraRight, etc.

        // Update GPU memory only if needed
        if (dirty) {
            glBindBuffer(GL_ARRAY_BUFFER, vboInst);
            
            // Check if we need to resize the buffer (if current text > pre-allocated size)
            // Or just orphan and upload new data.
            // "Orphaning" (pass NULL) tells driver to discard old buffer and give us a new pointer immediately.
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(GlyphInstance), nullptr, GL_STREAM_DRAW);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(GlyphInstance), instances.data(), GL_STREAM_DRAW);
            
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            dirty = false;
        }

        // --- THE FAST DRAW CALL ---
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fontTex);
        shader.setUniformi("uFontTex", 0);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // Bind ONE object (VAO) containing all state
        glBindVertexArray(vao);
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, instances.size());
        glBindVertexArray(0);

        glDisable(GL_BLEND);
    }

private:
    void initShader() {
        // (Use the shader code from previous answer, it was correct)
        // ...
    }
    void setupUniforms(const Camera& cam) {
        // (Matrix calculation logic)
        // ...
    }
};

#endif
```

---

### 3. Implementation: Efficient Mesh Renderer

Since `MeshRenderOGL3` is new, let's make it efficient. It should accept your raw data pointers but wrap them in a VAO internally.

#### `MeshRenderOGL3.h`

```cpp
#ifndef MeshRenderOGL3_h
#define MeshRenderOGL3_h

#include "GLobjects.h" // Assuming this has GLMesh
#include "Shader.h"

class MeshRenderOGL3 {
public:
    // We keep track of a VAO for the mesh
    GLuint vaoTri = 0;
    GLuint vaoLine = 0;
    
    // We assume GLMesh holds the actual VBO IDs (buffers)
    GLMesh* mesh_tri = nullptr;
    GLMesh* mesh_lines = nullptr;

    Shader sh_solid;
    Shader sh_const;

    // ... Transforms ...

    void initDefaultShaders() { /* ... */ }

    // When you upload mesh data, we Record the VAO state immediately.
    void uploadMesh_d(int nVerts, int nTris, const int* tris, const double* verts, const double* nors) {
        if (!mesh_tri) mesh_tri = new GLMesh();
        
        // 1. Upload data to GPU (creates VBOs inside GLMesh)
        mesh_tri->init_d(nVerts, nTris*3, (int*)tris, (double*)verts, (double*)nors, nullptr, nullptr);
        mesh_tri->draw_mode = GL_TRIANGLES;

        // 2. Create/Recreate VAO for this mesh configuration
        if(vaoTri) glDeleteVertexArrays(1, &vaoTri);
        glGenVertexArrays(1, &vaoTri);
        glBindVertexArray(vaoTri);

        // 3. Bind the VBOs created by GLMesh and configure attributes
        // Assuming GLMesh stores VBO ID in 'buffs[0]' for pos, 'buffs[1]' for nor, etc.
        // You might need to check GLMesh internal structure. 
        // Typically:
        
        // Position
        glBindBuffer(GL_ARRAY_BUFFER, mesh_tri->buffs[0]); 
        glEnableVertexAttribArray(0); // Layout 0 in Shader
        // Note: double* input usually converted to float on GPU, or use GL_DOUBLE if shader expects dvec3
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0); 

        // Normal
        if (mesh_tri->buffs[1]) {
            glBindBuffer(GL_ARRAY_BUFFER, mesh_tri->buffs[1]);
            glEnableVertexAttribArray(1); // Layout 1 in Shader
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
        }
        
        // Indices
        if (mesh_tri->buffs[3]) { // indices typically
             glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri->buffs[3]);
        }

        glBindVertexArray(0); // Done recording
    }

    void draw(const Camera& cam) {
        // Calculate matrices...

        if (mesh_tri && vaoTri) {
            sh_solid.use();
            // set uniforms...
            
            // Fast Draw
            glBindVertexArray(vaoTri);
            glDrawElements(GL_TRIANGLES, mesh_tri->count, GL_UNSIGNED_INT, nullptr);
            glBindVertexArray(0);
        }
    }
};
#endif
```

### Summary of Changes for Efficiency

1.  **VAO Encapsulation:** Both renderers now create a `vao` in `init` (or `upload`) and bind it in `draw`. This solves the attribute leakage problem completely.
2.  **No Redundant Calls:** The `draw` loops are now just: `UseShader` -> `SetUniforms` -> `BindVAO` -> `Draw`. Minimal driver communication.
3.  **Data Oriented:**
    *   Text uses a contiguous `std::vector` + `glBufferData(Stream)`.
    *   Mesh uses static VBOs baked into a VAO.
4.  **Flexibility:** You can still use `text.addLabel(mesh.points[i])`. It's user-friendly. The performance hit of copying the vector is non-existent compared to the gain from using Instanced Drawing with VAOs.