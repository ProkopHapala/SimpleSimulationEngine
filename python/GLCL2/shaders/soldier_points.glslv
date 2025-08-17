#version 330 core
layout (location = 0) in vec4 pos_dir; // xy = position (NDC), zw = orientation unit vector
void main() {
    gl_Position = vec4(pos_dir.xy, 0.0, 1.0);
    gl_PointSize = 5.0; // default size
}
