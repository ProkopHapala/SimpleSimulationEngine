varying vec2 vUv;
varying float vAlpha;

uniform sampler2D uFontTex;
uniform vec3 uColor;

void main() {
    if (vAlpha < 0.5) discard;

    vec4 texColor = texture2D(uFontTex, vUv);
    
    // Assuming font texture is white text on transparent background
    // or grayscale.
    // If grayscale (alpha only), use texColor.a or texColor.r
    
    float alpha = texColor.a;
    if (alpha < 0.1) discard;
    
    gl_FragColor = vec4(uColor, alpha);
}
