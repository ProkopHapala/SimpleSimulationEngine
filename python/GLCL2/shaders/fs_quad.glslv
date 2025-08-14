#version 330 core
layout(location=0) in vec2 pos;
out vec2 uv;
void main(){
    uv = pos*0.5 + 0.5; // map from [-1,1] to [0,1]
    gl_Position = vec4(pos, 0.0, 1.0);
}
