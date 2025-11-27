uniform sampler2D uPosTex;
uniform vec2 uTexSize;

// Palette of possible bond/edge colors; small fixed-size array is fine.
// The JS side should set uMatColors[matID] for each material type.
uniform vec4 uMatColors[8];

attribute float aAtomID;
attribute float aMatID;

varying vec4 vColor;

void main() {
    float tx = (mod(aAtomID, uTexSize.x) + 0.5) / uTexSize.x;
    float ty = (floor(aAtomID / uTexSize.x) + 0.5) / uTexSize.y;
    vec3 pos = texture2D(uPosTex, vec2(tx, ty)).xyz;

    int idx = int(aMatID + 0.5);
    //idx = clamp(idx, 0, 7);
    vColor = uMatColors[idx];

    vec4 mvPosition = modelViewMatrix * vec4(pos, 1.0);
    gl_Position = projectionMatrix * mvPosition;
}
