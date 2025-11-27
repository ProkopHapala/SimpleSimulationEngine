// Per-instance attributes (instanced on JS side via createLabelInstancedMesh):
// - aAtomID : which atom this label is anchored to
// - aLabel1 : packed ASCII codes for characters 0-3
// - aLabel2 : packed ASCII codes for characters 4-7
// - aStrLen : actual string length (<= 8), used for centering
attribute float aAtomID;
attribute vec4 aLabel1; // Chars 0-3 (per-instance)
attribute vec4 aLabel2; // Chars 4-7 (per-instance)
attribute float aStrLen; // Length of the string for this instance (per-instance)

// Per-vertex attribute (varies across the 8 character quads of an instance):
// - aCharPos: 0..7, selects which of the packed characters this quad should render
attribute float aCharPos; // 0..7, per-vertex

varying vec2 vUv;
varying float vAlpha;

uniform sampler2D uPosTex;
uniform vec2 uTexSize; // (width, height) of pos texture
uniform vec2 uFontGrid; // e.g. (16, 8)
uniform float uScale;
uniform bool uScreenSpace;
uniform float uAspect;

// Helper: select component from vec4 by small float index in [0,3].
// This encodes the mapping 0->x, 1->y, 2->z, 3->w in a readable way.
float indexToVec4(vec4 v, float idx) {
    if   (idx < 1.5){ return (idx<0.5) ? v.x: v.y; }
    else            { return (idx<2.5) ? v.z: v.w; }
}

void main() {
    // 1. Determine Character Index (ASCII code) for this quad.
    // aCharPos is 0..7 across the 8 character quads of an instance.
    float charIndex = 0.0;
    float idx = floor(aCharPos + 0.5); // round to nearest int in [0,7]
    if (idx < 4.0) { charIndex = indexToVec4(aLabel1, idx      ); } // Characters 0..3 packed in aLabel1
    else           { charIndex = indexToVec4(aLabel2, idx - 4.0); } // Characters 4..7 packed in aLabel2

    // If charIndex is 0, hide
    if (charIndex < 1.0) {
        vAlpha = 0.0;
        gl_Position = vec4(2.0, 2.0, 2.0, 1.0); // Clip
        return;
    }
    vAlpha = 1.0;

    // 2. Fetch Atom Position
    float tx = (mod(aAtomID, uTexSize.x) + 0.5) / uTexSize.x;
    float ty = (floor(aAtomID / uTexSize.x) + 0.5) / uTexSize.y;
    vec3 atomPos = texture2D(uPosTex, vec2(tx, ty)).xyz;

    // 3. UV Calculation
    float cols = uFontGrid.x;
    float rows = uFontGrid.y;
    float col = mod(charIndex, cols);
    float row = floor(charIndex / cols);
    float cellW = 1.0 / cols;
    float cellH = 1.0 / rows;

    // Shrink glyph horizontally inside the atlas cell to eliminate side padding.
    // uv.x is originally in [0,1]; remap to [0.25, 0.75] (centered), so the visible
    // glyph uses only the middle part of the cell.
    float glyphFill = 0.7;           // fraction of cell width used for glyph
    float localX = 0.5 + (uv.x - 0.5) * glyphFill;

    float u = (col + localX) * cellW;
    float v = (rows - 1.0 - row + uv.y) * cellH;
    vUv = vec2(u, v);

    // 4. View-Space Billboarding
    vec4 mvPosition = viewMatrix * vec4(atomPos, 1.0);
    vec4 projectedPos = projectionMatrix * mvPosition;

    // Character layout:
    // - Quad size stays 1.0 in X (position.x in [-0.5,0.5]) to avoid distorting glyphs.
    // - charAdvance controls spacing between character centers (< 1.0 to move them closer).
    float charAdvance = 0.6;   // horizontal spacing between character centers

    float centerOffset = (aStrLen - 1.0) * 0.5;
    float charOffset   = (aCharPos - centerOffset) * charAdvance;

    vec2 baseOffset = vec2(charOffset + position.x, position.y);

    if (uScreenSpace) {
        // Screen Space Scaling (Constant size on screen)
        // uScale is interpreted as size relative to screen height
        
        vec2 screenOffset = baseOffset * uScale;
        screenOffset.x /= uAspect; // Correct for aspect ratio
        
        // Apply to projected position (Clip Space)
        // We multiply by w because gl_Position will be divided by w
        projectedPos.xy += screenOffset * projectedPos.w;

    } else {
        // World Space Scaling (Standard 3D)
        vec2 viewOffset = baseOffset * uScale;
        mvPosition.xy += viewOffset;
        projectedPos = projectionMatrix * mvPosition;
    }

    gl_Position = projectedPos;
}
