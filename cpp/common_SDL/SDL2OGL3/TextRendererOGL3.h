#ifndef TextRendererOGL3_h
#define TextRendererOGL3_h

#include "Camera.h"
#include "Mat4.h"
#include "Shader.h"
#include "Vec3.h"
#include <GL/glew.h>
#include <vector>

// Pure Data Struct
struct GlyphInstance {
  Vec3f pos;
  int charIndex;
  int ascii;
};

class TextRendererOGL3 {
public:
  // OpenGL State Handles
  GLuint vao = 0;
  GLuint vboQuad = 0;
  GLuint vboInst = 0; // The dynamic buffer
  GLuint ebo = 0;
  GLuint fontTex = 0;

  Shader shader;

  // Config
  float charW = 0.1f;
  float charH = 0.1f;
  int numGlyphs = 128;
  int glyphOffset = 0;

  // CPU Data Store
  std::vector<GlyphInstance> instances;
  bool dirty = false;

  Vec3f modelPos;
  Mat3f modelMat;

  void setModelPos(const Vec3f &pos) { modelPos = pos; }
  void setModelMat(const Mat3f &mat) { modelMat = mat; }

  // --- 1. Initialization (The "Setup Once" phase) ---
  void init(GLuint fontTex_, float charW_, float charH_) {
    fontTex = fontTex_;
    charW = charW_;
    charH = charH_;

    modelPos.set(0.0f, 0.0f, 0.0f);
    modelMat.setOne();

    initShader();

    // 1. Create the VAO first. All subsequent state is recorded into it.
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // 2. Setup Static Geometry (The Quad)
    // [0,0] -> [1,1]
    const float quadVerts[] = {0, 0, 1, 0, 1, 1, 0, 1};
    const GLuint quadIdx[] = {0, 1, 2, 0, 2, 3};

    glGenBuffers(1, &vboQuad);
    glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);

    glGenBuffers(1, &ebo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIdx), quadIdx,
                 GL_STATIC_DRAW);

    // Attribute 0: Quad Position (2 floats)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void *)0);

    // 3. Setup Dynamic Instance Buffer
    glGenBuffers(1, &vboInst);
    glBindBuffer(GL_ARRAY_BUFFER, vboInst);
    // Pre-allocate memory (e.g., for 10,000 chars) to avoid resizing later.
    // GL_STREAM_DRAW indicates we modify it every frame but draw it few times.
    glBufferData(GL_ARRAY_BUFFER, 10000 * sizeof(GlyphInstance), nullptr,
                 GL_STREAM_DRAW);

    // Attribute 1: World Position (3 floats, Instanced)
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(GlyphInstance),
                          (void *)offsetof(GlyphInstance, pos));
    glVertexAttribDivisor(1, 1); // IMPORTANT: Update once per instance

    // Attribute 2: Char Index (1 int, Instanced)
    glEnableVertexAttribArray(2);
    glVertexAttribIPointer(2, 1, GL_INT, sizeof(GlyphInstance),
                           (void *)offsetof(GlyphInstance, charIndex));
    glVertexAttribDivisor(2, 1);

    // Attribute 3: ASCII Code (1 int, Instanced)
    glEnableVertexAttribArray(3);
    glVertexAttribIPointer(3, 1, GL_INT, sizeof(GlyphInstance),
                           (void *)offsetof(GlyphInstance, ascii));
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

  void addLabel(const Vec3f &basePos, const char *text) {
    // Very fast CPU copy.
    // Modern CPUs can push vectors like this extremely fast.
    for (int i = 0; text[i]; i++) {
      // Emplace back constructs the object directly in the vector memory
      instances.push_back({basePos, i, (int)text[i]});
    }
    dirty = true;
  }

  // --- 3. Rendering (The "Draw Cheaply" phase) ---
  void draw(const Camera &cam) {
    if (instances.empty() || !fontTex)
      return;

    shader.use();
    setupUniforms(cam); // Helper to set VP, CameraRight, etc.

    // Update GPU memory only if needed
    if (dirty) {
      glBindBuffer(GL_ARRAY_BUFFER, vboInst);

      // Check if we need to resize the buffer (if current text > pre-allocated
      // size) Or just orphan and upload new data. "Orphaning" (pass NULL) tells
      // driver to discard old buffer and give us a new pointer immediately.
      glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(GlyphInstance),
                   nullptr, GL_STREAM_DRAW);
      glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(GlyphInstance),
                   instances.data(), GL_STREAM_DRAW);

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
    glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0,
                            instances.size());
    glBindVertexArray(0);

    glDisable(GL_BLEND);
  }

private:
  void initShader() {
    const char *vertSrc = GLSL(
        330, layout(location = 0) in vec2 aQuadPos;
        layout(location = 1) in vec3 aBasePos;
        layout(location = 2) in int aCharIndex;
        layout(location = 3) in int aAscii;

        uniform mat4 uVP; uniform vec3 uCamRight; uniform vec3 uCamUp;
        uniform vec3 uCamPos; uniform vec3 uModelPos; uniform mat3 uModelMat;
        uniform float uCharW; uniform float uCharH; uniform float uAtlasStep;
        uniform int uGlyphOffset;

        out vec2 vUV;

        void main() {
          // World Space Calculation
          vec3 worldAnchor = uModelMat * aBasePos + uModelPos;

          float xOffset = float(aCharIndex) * uCharW;
          vec3 center = worldAnchor + (uCamRight * xOffset);

          // Billboard expansion
          vec3 vertexPos = center +
                           uCamRight * (aQuadPos.x - 0.5) *
                               uCharW // Centering adjustment if needed
                           + uCamUp * (aQuadPos.y - 0.5) * uCharH;

          // Texture Space Calculation
          float glyphIndex = float(aAscii - uGlyphOffset);
          float uStart = glyphIndex * uAtlasStep;
          vUV = vec2(uStart + (aQuadPos.x * uAtlasStep), aQuadPos.y);

          gl_Position = uVP * vec4(vertexPos - uCamPos, 1.0);
        });

    const char *fragSrc = GLSL(
        330, in vec2 vUV; out vec4 FragColor;

        uniform sampler2D uFontTex; uniform vec4 uTextColor;

        void main() {
          vec4 texCol = texture(uFontTex, vUV);
          if (texCol.a < 0.1)
            discard;
          FragColor = uTextColor * texCol;
        });

    shader.init_str(vertSrc, fragSrc, nullptr);
  }

  void setupUniforms(const Camera &cam) {
    // Build camera matrix and right/up vectors from Camera (same as
    // MeshRenderOGL3)
    Mat4f camMat, mRot, mPersp;
    if (cam.persp) {
      mPersp.setPerspective(cam.aspect * cam.zoom, cam.zoom, cam.zmin,
                            cam.zmax);
      mRot.setOne();
      mRot.setRot(cam.rot);
      camMat.set_mmul_TN(mRot, mPersp);
    } else {
      camMat.setOrthographic(cam.zoom, cam.zoom * cam.aspect, cam.zmin,
                             cam.zmax);
    }

    // Extract camera basis from rotation
    // Assuming Mat3T stores columns a,b,c as Right, Up, Forward/Back
    // However, Mat3.h aliases 'a' to 'lf' (Left?), so we might need to invert
    // it to get Right. User reported opposite rotation, suggesting Right vector
    // is inverted.
    Vec3f camRight = cam.rot.a * -1.0f;
    Vec3f camUp = cam.rot.b;

    shader.setUniformMat4f("uVP", camMat);
    shader.setUniformVec3f("uCamRight", camRight);
    shader.setUniformVec3f("uCamUp", camUp);
    shader.setUniformVec3f("uCamPos", cam.pos);

    // Default Model Transform (Identity) if not set externally,
    // but better to add setModelPos/Mat methods.
    // For now, we use the member variables if we add them, or just Identity.
    // Let's add member variables to be safe and consistent.
    shader.setUniformVec3f("uModelPos", modelPos);
    shader.setUniformMat3f("uModelMat", modelMat);

    shader.setUniformf("uCharW", charW);
    shader.setUniformf("uCharH", charH);
    shader.setUniformf("uAtlasStep", 1.0f / (float)numGlyphs);
    shader.setUniformi("uGlyphOffset", glyphOffset);
    shader.setUniformVec4f("uTextColor", (Quat4f){1.0f, 1.0f, 1.0f, 1.0f});
  }
};

// Create a simple dummy font texture (white boxes) so that geometry can be
// tested
static GLuint makeDummyFontTex(int numGlyphs) {
  const int nx = numGlyphs;
  const int ny = 1;
  std::vector<unsigned char> data(nx * ny * 4);
  for (int i = 0; i < nx * ny; i++) {
    data[i * 4 + 0] = 255;
    data[i * 4 + 1] = 255;
    data[i * 4 + 2] = 255;
    data[i * 4 + 3] = 255;
  }
  GLuint tex = 0;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, nx, ny, 0, GL_RGBA, GL_UNSIGNED_BYTE,
               data.data());
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  return tex;
}

#endif
