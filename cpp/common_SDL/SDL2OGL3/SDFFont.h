#ifndef SDFFont_h
#define SDFFont_h

/// === AUTO-DOC BEGIN ===
/// @file SDFFont.h
/// @brief SDF (Signed Distance Field) font rendering — crisp text at any zoom from low-res textures.
///
/// Stores signed distance-to-edge in a tiny texture (e.g. 64x64 per glyph)
/// instead of binary alpha in a high-res atlas. The fragment shader uses
/// `smoothstep(0.5 - sm, 0.5 + sm, dist)` with screen-space adaptive smoothing
/// via `dFdx`/`dFdy` derivatives, giving anti-aliased edges at any zoom level.
///
/// **Pipeline:**
/// - `makeSDFFontAtlas()` — renders TTF glyphs via SDL2_ttf at 8x resolution,
///   computes two-pass Felzenszwalb-Huttenlocher distance transform, downsamples
///   to low-res SDF texture, uploads as GL_R32F strip atlas
/// - `SDFFontRenderer` — instanced billboard renderer (same pattern as
///   **TextRendererOGL3** but with SDF shader instead of binary alpha discard)
///
/// **Key insight:** The SDF encoding is resolution-independent — a 64x64 glyph
/// texture produces the same crisp edges whether the glyph is 16px or 1600px
/// on screen. The `spread` parameter controls how far outside the glyph the
/// distance field extends (larger spread = softer outline glow, but less
/// precision inside the glyph).
///
/// **Based on:** Chris Green (Valve), "Improved Alpha-Tested Magnification for
/// Vector Textures and Special Effects", SIGGRAPH 2007.
/// Extension: msdfgen (Viktor Chlumský) — 3-channel SDF for sharper corners.
///
/// **Used by:** `cpp/sketches_SDL/OGL3/test_SDFFont.cpp` — interactive demo
/// === AUTO-DOC END ===

#include "Camera.h"
#include "Mat4.h"
#include "Shader.h"
#include "Vec3.h"
#include <GL/glew.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_ttf.h>
#include <vector>
#include <cstring>
#include <cmath>
#include <algorithm>

// ============================================================
// CPU-side SDF Atlas Generation
// ============================================================

/// @brief 1D squared Euclidean distance transform (Felzenszwalb & Huttenlocher).
///
/// Lower envelope of parabolas — O(n) using intersection-of-parabolas algorithm.
/// Called twice (rows then columns) by `distanceTransform2D` for exact 2D EDT.
/// `f` = per-point cost (0 for object, 1e30 for non-object), `d` = output squared distance.
inline void distanceTransform1D(const float* f, int n, float* d, int* v, float* z) {
    int k = 0;
    v[0] = 0;
    z[0] = -1e30f;
    z[1] = 1e30f;
    for (int q = 1; q < n; q++) {
        float s;
        while (k >= 0) {
            int vp = v[k];
            float denom = 2.0f * (float)(q - vp);
            if (denom == 0.0f) denom = 1e-10f;
            s = ((float)q*(float)q - (float)vp*(float)vp + f[q] - f[vp]) / denom;
            if (s <= z[k]) k--;
            else break;
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = 1e30f;
    }
    k = 0;
    for (int q = 0; q < n; q++) {
        while (z[k + 1] < (float)q) k++;
        float dist = (float)q - (float)v[k];
        d[q] = dist * dist + f[v[k]];
    }
}

/// @brief 2D squared Euclidean distance transform of a binary image.
///
/// Separable: applies `distanceTransform1D` along rows then columns.
/// Non-zero pixels get distance 0; zero pixels get distance to nearest non-zero.
inline void distanceTransform2D(const unsigned char* img, int W, int H, float* out) {
    std::vector<float> f(std::max(W, H));
    std::vector<float> d(std::max(W, H));
    std::vector<int>   v(std::max(W, H));
    std::vector<float> z(std::max(W, H) + 1);

    // Pass 1: along rows (x)
    std::vector<float> rowSqrDist(W * H);
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            f[x] = img[y * W + x] ? 0.0f : 1e30f;
        }
        distanceTransform1D(f.data(), W, d.data(), v.data(), z.data());
        for (int x = 0; x < W; x++) {
            rowSqrDist[y * W + x] = d[x];
        }
    }
    // Pass 2: along columns (y)
    for (int x = 0; x < W; x++) {
        for (int y = 0; y < H; y++) {
            f[y] = rowSqrDist[y * W + x];
        }
        distanceTransform1D(f.data(), H, d.data(), v.data(), z.data());
        for (int y = 0; y < H; y++) {
            out[y * W + x] = d[y]; // squared distance
        }
    }
}

/// @brief Generate single-channel SDF atlas from pre-rendered binary glyph images.
///
/// Computes signed distance (positive inside, negative outside) via two distance
/// transforms (inside→outside, outside→inside), normalizes to [0,1] with 0.5 at edge,
/// downsamples via bilinear interpolation.
/// Output: `atlas` is loW*nGlyphs × loH, values [0,1]. 0.5=edge, >0.5=inside, <0.5=outside.
inline void generateSDFFontAtlas(
    const unsigned char* glyphs, int nGlyphs, int hiW, int hiH,
    int loW, int loH, float spread,
    float* atlas  // output: loW*nGlyphs * loH, values [0,1]
) {
    std::vector<float> dist(hiW * hiH);
    std::vector<float> outDist(hiW * hiH);
    float maxDist = (float)spread;

    for (int g = 0; g < nGlyphs; g++) {
        const unsigned char* glyph = glyphs + g * hiW * hiH;

        // Distance from inside pixels to nearest outside (distanceTransform of inverted)
        std::vector<unsigned char> inv(hiW * hiH);
        for (int i = 0; i < hiW * hiH; i++) inv[i] = glyph[i] ? 0 : 1;
        std::vector<float> distOut(hiW * hiH);
        distanceTransform2D(inv.data(), hiW, hiH, distOut.data()); // dist from outside to nearest inside

        // Distance from outside pixels to nearest inside
        std::vector<float> distIn(hiW * hiH);
        distanceTransform2D(glyph, hiW, hiH, distIn.data()); // dist from inside to nearest outside

        // Signed distance: positive inside, negative outside
        for (int i = 0; i < hiW * hiH; i++) {
            float dOut = std::sqrt(distOut[i]); // distance from outside pixel to nearest inside
            float dIn  = std::sqrt(distIn[i]);  // distance from inside pixel to nearest outside
            float sd = glyph[i] ? dIn : -dOut;  // positive inside, negative outside
            // Normalize to [0,1]: 0.5 at edge, 1.0 at spread inside, 0.0 at spread outside
            float normalized = 0.5f + 0.5f * (sd / maxDist);
            outDist[i] = std::clamp(normalized, 0.0f, 1.0f);
        }

        // Downsample: sample the hi-res SDF at lo-res grid positions
        // Simple box filter averaging
        for (int ly = 0; ly < loH; ly++) {
            for (int lx = 0; lx < loW; lx++) {
                float hx = (float)lx / (float)(loW - 1) * (float)(hiW - 1);
                float hy = (float)ly / (float)(loH - 1) * (float)(hiH - 1);
                int ix = (int)hx, iy = (int)hy;
                int ix1 = std::min(ix + 1, hiW - 1);
                int iy1 = std::min(iy + 1, hiH - 1);
                float fx = hx - (float)ix, fy = hy - (float)iy;
                // Bilinear interpolation
                float v00 = outDist[iy  * hiW + ix ];
                float v01 = outDist[iy  * hiW + ix1];
                float v10 = outDist[iy1 * hiW + ix ];
                float v11 = outDist[iy1 * hiW + ix1];
                float v0 = v00 + (v01 - v00) * fx;
                float v1 = v10 + (v11 - v10) * fx;
                float val = v0 + (v1 - v0) * fy;
                atlas[ly * (loW * nGlyphs) + g * loW + lx] = val;
            }
        }
    }
}

/// @brief Upload SDF atlas as a horizontal texture strip (nGlyphs*loW × loH).
///
/// Uses GL_R32F (single-channel float) with linear filtering and clamp-to-edge.
/// Linear filtering is essential — it interpolates the distance field smoothly.
inline GLuint uploadSDFAtlas(const float* atlas, int nGlyphs, int loW, int loH) {
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, loW * nGlyphs, loH, 0, GL_RED, GL_FLOAT, atlas);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, 0);
    return tex;
}

// ============================================================
// SDF Font Renderer (instanced billboard, same pattern as TextRendererOGL3)
// ============================================================

/// @brief Per-instance data for one glyph quad.
/// Each text character becomes one instance: world-space anchor position,
/// character index within its string (for x-offset), and ASCII code.
struct SDFFontInstance {
    Vec3f pos;        ///< World-space anchor position of the label
    int   charIndex;  ///< Character index within string (for horizontal offset)
    int   ascii;      ///< ASCII code of the character
};

/// @brief Instanced billboard SDF font renderer.
///
/// Same architecture as **TextRendererOGL3** (VAO + static quad + streaming
/// instance VBO with `GL_STREAM_DRAW` orphaning), but uses an SDF fragment
/// shader with screen-space adaptive `smoothstep` instead of binary alpha discard.
///
/// **Usage pattern:**
/// 1. `init()` with atlas texture from `makeSDFFontAtlas()`
/// 2. Each frame: `clear()` → `addLabel()` × N → `draw(cam)`
///
/// **Shader:** Vertex shader positions billboarded quads in world space using
/// camera right/up vectors (always faces camera). Fragment shader samples the
/// SDF texture and applies `smoothstep(0.5-sm, 0.5+sm, dist)` where `sm` is
/// derived from `dFdx`/`dFdy` screen-space derivatives — ~1 texel wide edge
/// at any zoom level.
class SDFFontRenderer {
public:
    GLuint vao = 0, vboQuad = 0, vboInst = 0, ebo = 0;
    GLuint fontTex = 0;
    Shader shader;

    float charW = 0.1f;   // world-space char width
    float charH = 0.1f;   // world-space char height
    int   nGlyphs = 95;   // printable ASCII range
    int   glyphOffset = 32; // ASCII offset (space = 32)
    int   atlasW = 0;     // atlas width  (nGlyphs * loW)
    int   atlasH = 0;     // atlas height (loH)
    int   perGlyphW = 32; // SDF texture width per glyph

    std::vector<SDFFontInstance> instances;
    bool dirty = false;

    Vec3f modelPos;
    Mat3f modelMat;

    /// @brief Set model-space offset for all labels (default: origin).
    void setModelPos(const Vec3f &pos) { modelPos = pos; }
    /// @brief Set model rotation matrix for all labels (default: identity).
    void setModelMat(const Mat3f &mat) { modelMat = mat; }

    /// @brief Initialize renderer: compile shader, create VAO/VBOs, set glyph metrics.
    /// `fontTex_` from `makeSDFFontAtlas()`, `charW_`/`charH_` = world-space glyph size.
    void init(GLuint fontTex_, int nGlyphs_, int glyphOffset_, int perGlyphW_, int atlasH_, float charW_, float charH_) {
        fontTex = fontTex_;
        nGlyphs = nGlyphs_;
        glyphOffset = glyphOffset_;
        perGlyphW = perGlyphW_;
        atlasW = nGlyphs * perGlyphW;
        atlasH = atlasH_;
        charW = charW_;
        charH = charH_;

        modelPos.set(0.0f, 0.0f, 0.0f);
        modelMat.setOne();

        initShader();

        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        // Static quad
        const float quadVerts[] = {0, 0, 1, 0, 1, 1, 0, 1};
        const GLuint quadIdx[] = {0, 1, 2, 0, 2, 3};
        glGenBuffers(1, &vboQuad);
        glBindBuffer(GL_ARRAY_BUFFER, vboQuad);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);
        glGenBuffers(1, &ebo);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIdx), quadIdx, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);

        // Dynamic instance buffer
        glGenBuffers(1, &vboInst);
        glBindBuffer(GL_ARRAY_BUFFER, vboInst);
        glBufferData(GL_ARRAY_BUFFER, 10000 * sizeof(SDFFontInstance), nullptr, GL_STREAM_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(SDFFontInstance), (void*)offsetof(SDFFontInstance, pos));
        glVertexAttribDivisor(1, 1);
        glEnableVertexAttribArray(2);
        glVertexAttribIPointer(2, 1, GL_INT, sizeof(SDFFontInstance), (void*)offsetof(SDFFontInstance, charIndex));
        glVertexAttribDivisor(2, 1);
        glEnableVertexAttribArray(3);
        glVertexAttribIPointer(3, 1, GL_INT, sizeof(SDFFontInstance), (void*)offsetof(SDFFontInstance, ascii));
        glVertexAttribDivisor(3, 1);

        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

    /// @brief Clear all labels (call at start of each frame).
    void clear() { instances.clear(); dirty = true; }

    /// @brief Add a text label at a 3D world position. Each character becomes one instance.
    void addLabel(const Vec3f &basePos, const char *text) {
        for (int i = 0; text[i]; i++) {
            instances.push_back({basePos, i, (int)text[i]});
        }
        dirty = true;
    }

    /// @brief Upload instances if dirty, bind shader/texture/VAO, draw instanced quads.
    /// Uses buffer orphaning for streaming (glBufferData with nullptr then real data).
    void draw(const Camera &cam) {
        if (instances.empty() || !fontTex) return;
        shader.use();
        setupUniforms(cam);

        if (dirty) {
            glBindBuffer(GL_ARRAY_BUFFER, vboInst);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(SDFFontInstance), nullptr, GL_STREAM_DRAW);
            glBufferData(GL_ARRAY_BUFFER, instances.size() * sizeof(SDFFontInstance), instances.data(), GL_STREAM_DRAW);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            dirty = false;
        }

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, fontTex);
        shader.setUniformi("uFontTex", 0);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glBindVertexArray(vao);
        glDrawElementsInstanced(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0, instances.size());
        glBindVertexArray(0);

        glDisable(GL_BLEND);
    }

private:
    void initShader() {
        const char *vertSrc = GLSL(330,
            layout(location = 0) in vec2 aQuadPos;
            layout(location = 1) in vec3 aBasePos;
            layout(location = 2) in int aCharIndex;
            layout(location = 3) in int aAscii;

            uniform mat4 uVP;
            uniform vec3 uCamRight;
            uniform vec3 uCamUp;
            uniform vec3 uCamPos;
            uniform vec3 uModelPos;
            uniform mat3 uModelMat;
            uniform float uCharW;
            uniform float uCharH;
            uniform float uAtlasStep;   // 1.0 / nGlyphs (width of one glyph in atlas UV)
            uniform int   uGlyphOffset;

            out vec2 vUV;
            out float vPixelScale;  // for adaptive smoothing

            void main() {
                vec3 worldAnchor = uModelMat * aBasePos + uModelPos;
                float xOffset = float(aCharIndex) * uCharW;
                vec3 center = worldAnchor + uCamRight * xOffset;
                vec3 vertexPos = center
                    + uCamRight * (aQuadPos.x - 0.5) * uCharW
                    + uCamUp    * (aQuadPos.y - 0.5) * uCharH;

                // Atlas UV: map ascii to horizontal strip in atlas
                float glyphIndex = float(aAscii - uGlyphOffset);
                float uStart = glyphIndex * uAtlasStep;
                vUV = vec2(uStart + aQuadPos.x * uAtlasStep, 1.0 - aQuadPos.y);

                // Estimate on-screen pixel size for adaptive edge smoothing
                // Project two points to get pixel-per-char ratio
                vec4 clipL = uVP * vec4(center - uCamPos, 1.0);
                vec4 clipR = uVP * vec4(center + uCamRight * uCharW - uCamPos, 1.0);
                float screenW = uCharW / max(abs(clipL.w), 0.001);
                vPixelScale = screenW;

                gl_Position = uVP * vec4(vertexPos - uCamPos, 1.0);
            }
        );

        const char *fragSrc = GLSL(330,
            in vec2 vUV;
            in float vPixelScale;
            out vec4 FragColor;

            uniform sampler2D uFontTex;
            uniform vec4  uTextColor;
            uniform float uAtlasStep;  // 1.0 / nGlyphs
            uniform vec2  uAtlasSize; // total atlas texture size (texels)

            void main() {
                // Sample SDF texture (R channel = signed distance, [0,1])
                float dist = texture(uFontTex, vUV).r;

                // Screen-space adaptive smoothing using derivatives.
                // fwidth gives us how much the UV changes per screen pixel.
                // Multiply by atlas size to get texels-per-screen-pixel.
                // The edge is at dist=0.5; smoothing band = ~1 texel wide on screen.
                vec2  dUV   = vec2(dFdx(vUV.x), dFdy(vUV.y));
                float pxSize = max(length(dUV * uAtlasSize), 1.0);
                float smoothing = 0.5 / pxSize * 2.0;  // ~1 texel wide edge
                float alpha = smoothstep(0.5 - smoothing, 0.5 + smoothing, dist);

                if (alpha < 0.01) discard;

                FragColor = vec4(uTextColor.rgb, uTextColor.a * alpha);
            }
        );

        shader.init_str(vertSrc, fragSrc, nullptr);
    }

    void setupUniforms(const Camera &cam) {
        Mat4f camMat, mRot, mPersp;
        if (cam.persp) {
            mPersp.setPerspective(cam.aspect * cam.zoom, cam.zoom, cam.zmin, cam.zmax);
            mRot.setOne();
            mRot.setRot(cam.rot);
            camMat.set_mmul_TN(mRot, mPersp);
        } else {
            camMat.setOrthographic(cam.zoom, cam.zoom * cam.aspect, cam.zmin, cam.zmax);
        }

        Vec3f camRight = cam.rot.a;
        Vec3f camUp = cam.rot.b;

        shader.setUniformMat4f("uVP", camMat);
        shader.setUniformVec3f("uCamRight", camRight);
        shader.setUniformVec3f("uCamUp", camUp);
        shader.setUniformVec3f("uCamPos", cam.pos);
        shader.setUniformVec3f("uModelPos", modelPos);
        shader.setUniformMat3f("uModelMat", modelMat);
        shader.setUniformf("uCharW", charW);
        shader.setUniformf("uCharH", charH);
        shader.setUniformf("uAtlasStep", 1.0f / (float)nGlyphs);
        shader.setUniformi("uGlyphOffset", glyphOffset);
        shader.setUniformVec4f("uTextColor", (Quat4f){1.0f, 1.0f, 1.0f, 1.0f});
        shader.setUniformVec2f("uAtlasSize", (Vec2f){(float)atlasW, (float)atlasH});
    }

};

// ============================================================
// Helper: generate a dummy SDF atlas for testing (no font needed)
// ============================================================

/// @brief Generate SDF font atlas from a real TTF font using SDL2_ttf.
///
/// Renders each glyph at 8x resolution via `TTF_RenderGlyph_Blended`, extracts
/// alpha channel via `SDL_GetRGBA` (handles ARGB/BGRA pixel formats), computes
/// signed distance field via two-pass distance transform, downsamples to
/// low-res SDF texture with bilinear interpolation.
///
/// **Critical detail:** Two distance transforms are needed — one from inside
/// pixels to nearest outside (gives positive distance for inside), one from
/// outside to nearest inside (gives negative distance for outside). The signed
/// distance is `mask ? dOut : -dIn` (positive inside, negative outside).
///
/// `spread` controls max distance encoded — larger = wider edge band but
/// less precision inside glyph. Typical: 8 texels at output resolution.
///
/// @return OpenGL texture ID (GL_R32F horizontal strip atlas), or 0 on error.
inline GLuint makeSDFFontAtlas(const char* fontPath, int nGlyphs, int glyphOffset,
                               int perGlyphW, int perGlyphH, float spread = 8.0f) {
    if (TTF_Init() < 0) {
        printf("SDFFont: ERROR - TTF_Init failed: %s\n", TTF_GetError());
        return 0;
    }

    // High-res rendering size: render large, then downsample to SDF
    int hiScale = 8;
    int hiW = perGlyphW * hiScale;
    int hiH = perGlyphH * hiScale;

    TTF_Font* font = TTF_OpenFont(fontPath, hiH * 3 / 4);  // font size ~75% of cell height
    if (!font) {
        printf("SDFFont: ERROR - cannot open font %s: %s\n", fontPath, TTF_GetError());
        return 0;
    }

    printf("SDFFont: rendering %d glyphs from %s at hi-res %dx%d -> SDF %dx%d\n",
           nGlyphs, fontPath, hiW, hiH, perGlyphW, perGlyphH);

    std::vector<float> atlas(perGlyphW * nGlyphs * perGlyphH, 0.0f);

    // Temporary buffers for distance transform
    std::vector<float> f(std::max(hiW, hiH));
    std::vector<float> d(std::max(hiW, hiH));
    std::vector<int>   v(std::max(hiW, hiH));
    std::vector<float> z(std::max(hiW, hiH) + 1);
    std::vector<float> rowSqrDist(hiW * hiH);
    std::vector<float> distIn(hiW * hiH);
    std::vector<float> distOut(hiW * hiH);
    std::vector<unsigned char> inv(hiW * hiH);

    float maxDist = spread * (float)hiScale;  // spread in hi-res pixels

    for (int g = 0; g < nGlyphs; g++) {
        char ch = (char)(glyphOffset + g);

        // Render glyph to SDL surface (white on black)
        SDL_Color white = {255, 255, 255, 255};
        SDL_Surface* surf = TTF_RenderGlyph_Blended(font, (Uint16)ch, white);
        if (!surf) {
            // Empty glyph (space) — fill with all-outside
            for (int ly = 0; ly < perGlyphH; ly++)
                for (int lx = 0; lx < perGlyphW; lx++)
                    atlas[ly * (perGlyphW * nGlyphs) + g * perGlyphW + lx] = 0.0f;
            continue;
        }

        // Convert to binary mask (alpha > 128 = inside)
        std::vector<unsigned char> mask(hiW * hiH, 0);
        SDL_LockSurface(surf);
        int srcW = surf->w;
        int srcH = surf->h;
        int srcPitch = surf->pitch;
        Uint8* srcPixels = (Uint8*)surf->pixels;
        SDL_PixelFormat* fmt = surf->format;
        // Center the glyph in the hi-res cell
        int ox = (hiW - srcW) / 2;
        int oy = (hiH - srcH) / 2;
        for (int y = 0; y < srcH && y + oy < hiH; y++) {
            for (int x = 0; x < srcW && x + ox < hiW; x++) {
                Uint32* px = (Uint32*)(srcPixels + y * srcPitch + x * fmt->BytesPerPixel);
                Uint8 r, g, b, a;
                SDL_GetRGBA(*px, fmt, &r, &g, &b, &a);
                if (a > 128) mask[(y + oy) * hiW + (x + ox)] = 255;
            }
        }
        SDL_UnlockSurface(surf);
        SDL_FreeSurface(surf);

        // Distance transform: inside -> nearest outside
        for (int i = 0; i < hiW * hiH; i++) inv[i] = mask[i] ? 0 : 1;
        // Pass 1: rows
        for (int y = 0; y < hiH; y++) {
            for (int x = 0; x < hiW; x++) f[x] = inv[y * hiW + x] ? 0.0f : 1e30f;
            distanceTransform1D(f.data(), hiW, d.data(), v.data(), z.data());
            for (int x = 0; x < hiW; x++) rowSqrDist[y * hiW + x] = d[x];
        }
        // Pass 2: cols
        for (int x = 0; x < hiW; x++) {
            for (int y = 0; y < hiH; y++) f[y] = rowSqrDist[y * hiW + x];
            distanceTransform1D(f.data(), hiH, d.data(), v.data(), z.data());
            for (int y = 0; y < hiH; y++) distOut[y * hiW + x] = d[y];
        }

        // Distance transform: outside -> nearest inside
        for (int y = 0; y < hiH; y++) {
            for (int x = 0; x < hiW; x++) f[x] = mask[y * hiW + x] ? 0.0f : 1e30f;
            distanceTransform1D(f.data(), hiW, d.data(), v.data(), z.data());
            for (int x = 0; x < hiW; x++) rowSqrDist[y * hiW + x] = d[x];
        }
        for (int x = 0; x < hiW; x++) {
            for (int y = 0; y < hiH; y++) f[y] = rowSqrDist[y * hiW + x];
            distanceTransform1D(f.data(), hiH, d.data(), v.data(), z.data());
            for (int y = 0; y < hiH; y++) distIn[y * hiW + x] = d[y];
        }

        // Compute signed distance and downsample to SDF
        for (int ly = 0; ly < perGlyphH; ly++) {
            for (int lx = 0; lx < perGlyphW; lx++) {
                // Sample hi-res SDF at bilinear interpolated position
                float hx = (float)lx / (float)(perGlyphW - 1) * (float)(hiW - 1);
                float hy = (float)ly / (float)(perGlyphH - 1) * (float)(hiH - 1);
                int ix = (int)hx, iy = (int)hy;
                int ix1 = std::min(ix + 1, hiW - 1);
                int iy1 = std::min(iy + 1, hiH - 1);
                float fx = hx - (float)ix, fy = hy - (float)iy;

                // Signed distance at 4 neighbors
                auto getSDF = [&](int x, int y) -> float {
                    float dOut = std::sqrt((float)distOut[y * hiW + x]); // dist to nearest outside
                    float dIn  = std::sqrt((float)distIn [y * hiW + x]); // dist to nearest inside
                    return mask[y * hiW + x] ? dOut : -dIn;  // +inside, -outside
                };
                float s00 = getSDF(ix,   iy);
                float s01 = getSDF(ix1,  iy);
                float s10 = getSDF(ix,   iy1);
                float s11 = getSDF(ix1,  iy1);
                float s0 = s00 + (s01 - s00) * fx;
                float s1 = s10 + (s11 - s10) * fx;
                float sd = s0 + (s1 - s0) * fy;

                // Normalize to [0,1]: 0.5 at edge
                float normalized = 0.5f + 0.5f * (sd / maxDist);
                atlas[ly * (perGlyphW * nGlyphs) + g * perGlyphW + lx] = std::clamp(normalized, 0.0f, 1.0f);
            }
        }
    }

    TTF_CloseFont(font);
    // Don't TTF_Quit() here — might be used elsewhere

    printf("SDFFont: atlas generated, %d glyphs, %dx%d each\n", nGlyphs, perGlyphW, perGlyphH);

    return uploadSDFAtlas(atlas.data(), nGlyphs, perGlyphW, perGlyphH);
}

/// @brief Generate dummy SDF atlas with rectangle glyphs (no font file needed).
///
/// Uses analytical signed distance to a centered rectangle — no distance transform.
/// Useful for testing the rendering pipeline without SDL2_ttf or a TTF file.
inline GLuint makeDummySDFFontAtlas(int nGlyphs, int perGlyphW, int perGlyphH, float spread = 0.25f) {
    float margin = 0.15f;
    float maxDist = spread;
    std::vector<float> atlas(perGlyphW * nGlyphs * perGlyphH, 0.0f);
    for (int g = 0; g < nGlyphs; g++) {
        for (int ly = 0; ly < perGlyphH; ly++) {
            for (int lx = 0; lx < perGlyphW; lx++) {
                float u = (float)lx / (float)(perGlyphW - 1);
                float v = (float)ly / (float)(perGlyphH - 1);
                float dx = std::min(u - margin, (1.0f - margin) - u);
                float dy = std::min(v - margin, (1.0f - margin) - v);
                float sd;
                if (dx > 0 && dy > 0) sd = std::min(dx, dy);
                else {
                    float ox = std::max(margin - u, u - (1.0f - margin));
                    float oy = std::max(margin - v, v - (1.0f - margin));
                    sd = -std::sqrt(std::max(ox,0.0f)*std::max(ox,0.0f) + std::max(oy,0.0f)*std::max(oy,0.0f));
                }
                atlas[ly * (perGlyphW * nGlyphs) + g * perGlyphW + lx] = std::clamp(0.5f + 0.5f * (sd / maxDist), 0.0f, 1.0f);
            }
        }
    }
    return uploadSDFAtlas(atlas.data(), nGlyphs, perGlyphW, perGlyphH);
}

#endif // SDFFont_h
