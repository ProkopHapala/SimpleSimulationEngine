attribute float aAtomID;
attribute vec4 aLabel1; // Chars 0-3
attribute vec4 aLabel2; // Chars 4-7
attribute float aCharPos; // 0..7, per-vertex
attribute float aStrLen; // Length of the string for this instance

varying vec2 vUv;
varying float vAlpha;

uniform sampler2D uPosTex;
uniform vec2 uTexSize; // (width, height) of pos texture
uniform vec2 uFontGrid; // e.g. (16, 8)
uniform float uScale;
uniform bool uScreenSpace;
uniform float uAspect;

void main() {
    // 1. Determine Character Index
    float charIndex = 0.0;
    if (aCharPos < 3.5) {
        if (aCharPos < 0.5) charIndex = aLabel1.x;
        else if (aCharPos < 1.5) charIndex = aLabel1.y;
        else if (aCharPos < 2.5) charIndex = aLabel1.z;
        else charIndex = aLabel1.w;
    } else {
        if (aCharPos < 4.5) charIndex = aLabel2.x;
        else if (aCharPos < 5.5) charIndex = aLabel2.y;
        else if (aCharPos < 6.5) charIndex = aLabel2.z;
        else charIndex = aLabel2.w;
    }

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
    
    float u = (col + uv.x) * cellW;
    float v = (rows - 1.0 - row + uv.y) * cellH;
    vUv = vec2(u, v);

    // 4. View-Space Billboarding
    vec4 mvPosition = viewMatrix * vec4(atomPos, 1.0);
    vec4 projectedPos = projectionMatrix * mvPosition;

    // Calculate offset
    float charWidth = 0.6; 
    float centerOffset = (aStrLen - 1.0) * 0.5;
    float charOffset = (aCharPos - centerOffset) * charWidth; 
    
    vec2 baseOffset = vec2(charOffset + position.x * charWidth, position.y);

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
