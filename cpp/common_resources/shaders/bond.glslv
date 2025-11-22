uniform sampler2D uPosTex;
uniform vec2 uTexSize;
attribute float aAtomID;

void main() {
    float tx = (mod(aAtomID, uTexSize.x) + 0.5) / uTexSize.x;
    float ty = (floor(aAtomID / uTexSize.x) + 0.5) / uTexSize.y;
    vec3 pos = texture2D(uPosTex, vec2(tx, ty)).xyz;
    
    vec4 mvPosition = modelViewMatrix * vec4(pos, 1.0);
    gl_Position = projectionMatrix * mvPosition;
}
