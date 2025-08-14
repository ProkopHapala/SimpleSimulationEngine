#version 330 core
in vec2 uv;
out vec4 fragColor;
uniform sampler2D iChannel0;
void main(){
    vec4 c = texture(iChannel0, uv);
    fragColor = c;
}
