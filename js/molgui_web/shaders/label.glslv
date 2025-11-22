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
    
    // Calculate offset in View Space
    // Center the string based on its length
    // aStrLen is the number of characters
    // We want to center the block of characters around 0
    
    float charWidth = 0.6; // Width of one character relative to height (aspect ratio)
    float centerOffset = (aStrLen - 1.0) * 0.5;
    float charOffset = (aCharPos - centerOffset) * charWidth; 
    
    // Apply scale
    vec2 offset = vec2(charOffset + (position.x - 0.5) * charWidth, position.y - 0.5) * uScale;
    // Note: position.x is 0..1 or -0.5..0.5? PlaneBufferGeometry(1,1) is centered at 0, so -0.5 to 0.5.
    // So (position.x) is enough if we want char center at charOffset.
    // But let's be precise:
    // We want the quad for char `i` to be at `charOffset`.
    // The quad itself spans -0.5 to 0.5.
    // So we just add position.x (which is the local quad coord).
    
    // Re-eval:
    // charOffset is the center of the character slot `i`.
    // position.x is the offset within that slot.
    // So `charOffset + position.x` is correct if charWidth is 1.0.
    // If charWidth is 0.6, we should scale position.x too? 
    // No, the texture is mapped to the quad. If we squish the quad, we squish the texture.
    // Font texture letters are usually not square. They are taller than wide.
    // Aspect ratio ~0.6.
    
    vec2 finalOffset = vec2(charOffset + position.x * charWidth, position.y) * uScale;
    
    mvPosition.xy += finalOffset;

    gl_Position = projectionMatrix * mvPosition;
}
