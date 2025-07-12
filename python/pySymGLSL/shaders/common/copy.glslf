uniform vec3      iResolution;
uniform sampler2D iChannel0;

void main(){
    gl_FragColor = textureLod(iChannel0, gl_FragCoord.xy/iResolution.xy, 0.);
}
