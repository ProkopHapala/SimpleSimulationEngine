attribute vec3 instanceColor;
attribute float instanceScale;

varying vec3 vColor;
varying vec2 vUv;
varying vec3 vViewPosition;
varying float vRadius;

uniform sampler2D uPosTex;
uniform vec2 uTexSize; // (width, height)
attribute float aAtomID;

void main() {
    vColor = instanceColor;
    vUv = uv;
    vRadius = instanceScale * 0.5;

    // Fetch position from texture
    // Calculate UV coordinate for the atom ID
    // Assuming texture is N x 1 or N x M
    // Let's assume 1D texture for simplicity first, or simple row-major mapping
    // But DataTexture is usually 2D. Let's assume width is set by renderer.
    
    float tx = (mod(aAtomID, uTexSize.x) + 0.5) / uTexSize.x;
    float ty = (floor(aAtomID / uTexSize.x) + 0.5) / uTexSize.y;
    
    vec3 atomPos = texture2D(uPosTex, vec2(tx, ty)).xyz;

    // Billboard technique:
    vec4 mvPosition = modelViewMatrix * vec4(atomPos, 1.0);
    
    // Offset by vertex position (scaled)
    mvPosition.xy += position.xy * instanceScale;
    
    vViewPosition = -mvPosition.xyz;
    
    gl_Position = projectionMatrix * mvPosition;
}
