#ifndef TextRendererOGL3_h
#define TextRendererOGL3_h

#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"

#include "Camera.h"

#include "GLobjects.h"
#include "Shader.h"


// =====================================================
//   Minimal OGL3 instanced-character text renderer
// =====================================================

struct GlyphInstance{
    Vec3f pos;       // base position of the label (world space)
    int   charIndex; // index of character in the string
    int   ascii;     // ASCII code of the character
};

class TextRendererOGL3{
public:
    // We do NOT use VAOs here on purpose, to match GLMesh style and avoid
    // interfering with whatever VAO/attribute state the rest of the engine uses.
    GLuint vboQuad  = 0;  // static quad vertices [0,1]x[0,1]
    GLuint vboInst  = 0;  // per-glyph instances
    GLuint ebo      = 0;  // quad indices
    GLuint fontTex  = 0;

    Shader shader;

    float charW = 0.1f;   // world-space width
    float charH = 0.1f;   // world-space height

    int   numGlyphs    = 96;   // e.g. 95 printable ASCII + margin
    int   glyphOffset  = 32;   // first ASCII in atlas (space)

    std::vector<GlyphInstance> instances;

    void init(GLuint fontTex_, float charW_, float charH_){
        fontTex = fontTex_;
        charW   = charW_;
        charH   = charH_;

        // Quad vertices in local [0,1] space
        const float quadVerts[8] = {
            0.0f, 0.0f,
            1.0f, 0.0f,
            1.0f, 1.0f,
            0.0f, 1.0f
        };
        const GLuint quadIdx[6] = {0,1,2, 0,2,3};

        glGenBuffers(1, &vboQuad);
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);

        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIdx), quadIdx, GL_STATIC_DRAW);

        glGenBuffers(1, &vboInst);
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);

        // Simple text shader (vertex + fragment) using 1D-style atlas in a 2D texture
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
                float xOffset = float(aCharIndex) * uCharW;
                vec3 charOrigin = aBasePos + uCamRight * xOffset;

                vec3 worldPos = charOrigin
                    + uCamRight * (aQuadPos.x * uCharW)
                    + uCamUp    * (aQuadPos.y * uCharH);

                float glyphIndex = float(aAscii - uGlyphOffset);
                float u0 = glyphIndex * uAtlasStep;
                float u1 = u0 + uAtlasStep;
                float u  = mix(u0, u1, aQuadPos.x);
                float v  = aQuadPos.y;
                vUV = vec2(u, v);

                gl_Position = uVP * vec4(worldPos, 1.0);
            }
        );

        const char* fragSrc = GLSL(330,
            in vec2 vUV;
            out vec4 FragColor;

            uniform sampler2D uFontTex;
            uniform vec4 uTextColor;

            void main(){
                vec4 texCol = texture(uFontTex, vUV);
                if(texCol.a < 0.1) discard;
                FragColor = uTextColor * texCol;
            }
        );

        shader.init_str(vertSrc, fragSrc, nullptr);
    }

    void clear(){ instances.clear(); }

    void addLabel(const Vec3f& basePos, const char* text){
        for(int i=0; text[i]; i++){
            GlyphInstance gi;
            gi.pos       = basePos;
            gi.charIndex = i;
            gi.ascii     = (unsigned char)text[i];
            instances.push_back(gi);
        }
    }

    void uploadInstances(){
        if(instances.empty()) return;
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        glBufferData(GL_ARRAY_BUFFER, instances.size()*sizeof(GlyphInstance), instances.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void draw(const Camera& cam){
        if(!fontTex) return;
        if(instances.empty()) return;

        // Build camera matrix and right/up vectors from Camera (same as MeshRenderOGL3)
        Mat4f camMat, mRot, mPersp;
        if(cam.persp){
            mPersp.setPerspective(cam.aspect*cam.zoom, cam.zoom, cam.zmin, cam.zmax);
            mRot.setOne();
            mRot.setRot(cam.rot);
            camMat.set_mmul_TN(mRot, mPersp);
        }else{
            camMat.setOrthographic(cam.zoom, cam.zoom*cam.aspect, cam.zmin, cam.zmax);
        }

        // Extract camera basis from rotation
        Vec3f camRight = cam.rot.b; // depending on convention in Camera/Quat; adjust if needed
        Vec3f camUp    = cam.rot.c;

        shader.use();
        shader.setUniformMat4f("uVP", camMat);
        shader.setUniformVec3f("uCamRight", camRight);
        shader.setUniformVec3f("uCamUp",    camUp);
        shader.setUniformf("uCharW", charW);
        shader.setUniformf("uCharH", charH);
        shader.setUniformf("uAtlasStep", 1.0f/(float)numGlyphs);
        shader.setUniformi("uGlyphOffset", glyphOffset);
        shader.setUniformVec4f("uTextColor", (Quat4f){1.0f,1.0f,1.0f,1.0f});

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fontTex);
        shader.setUniformi("uFontTex", 0);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // Set up vertex attributes explicitly, like GLMesh::preDraw
        // Quad vertices (location 0)
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0);

        // Instance data (locations 1,2,3)
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

        // Index buffer
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);

        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, instances.size());

        // Clean up attrib state
        glDisableVertexAttribArray(3);
        glDisableVertexAttribArray(2);
        glDisableVertexAttribArray(1);
        glDisableVertexAttribArray(0);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        glDisable(GL_BLEND);
    }
};

// Create a simple dummy font texture (white boxes) so that geometry can be tested
static GLuint makeDummyFontTex(int numGlyphs){
    const int nx = numGlyphs;
    const int ny = 1;
    std::vector<unsigned char> data(nx*ny*4);
    for(int i=0;i<nx*ny;i++){
        data[i*4+0] = 255;
        data[i*4+1] = 255;
        data[i*4+2] = 255;
        data[i*4+3] = 255;
    }
    GLuint tex=0;
    glGenTextures(1,&tex);
    glBindTexture(GL_TEXTURE_2D,tex);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,nx,ny,0,GL_RGBA,GL_UNSIGNED_BYTE,data.data());
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    return tex;
}

#endif
