#version 330 core
layout (location = 0) in vec2 pos;      // NDC position
layout (location = 1) in vec4 inColor;  // RGBA per-vertex color
out vec4 vColor;
void main(){
    vColor = inColor;
    gl_Position = vec4(pos, 0.0, 1.0);
}
