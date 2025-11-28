attribute vec3 instanceColor;
attribute float instanceScale;

varying vec3 vColor;
varying vec2 vUv;
varying vec3 vViewPosition;
varying float vRadius;

uniform sampler2D uPosTex;
uniform vec2 uTexSize; // (width, height)
attribute float aAtomID;

// Global point size scale (set separately for normal and selected vertices
// via material uniforms, without touching textures or instance attributes).
uniform float uPointScale;

void main() {
    vColor = instanceColor;
    vUv = uv;

    // Fetch position (xyz) and base radius (w) from texture
    float tx = (mod(aAtomID, uTexSize.x) + 0.5) / uTexSize.x;
    float ty = (floor(aAtomID / uTexSize.x) + 0.5) / uTexSize.y;
    vec4 posTexel = texture2D(uPosTex, vec2(tx, ty));
    vec3 atomPos = posTexel.xyz;
    float baseRadius = posTexel.w; // default set from JS, typically 1.0

    // Effective scale combines:
    //   - per-instance scale (instanceScale)
    //   - per-vertex base radius from texture (baseRadius)
    //   - global uniform scale (uPointScale)
    float effScale = instanceScale * baseRadius * uPointScale;
    vRadius = effScale * 0.5;

    // Billboard technique:
    vec4 mvPosition = modelViewMatrix * vec4(atomPos, 1.0);

    // Offset by vertex position (scaled)
    mvPosition.xy += position.xy * effScale;

    vViewPosition = -mvPosition.xyz;

    gl_Position = projectionMatrix * mvPosition;
}
