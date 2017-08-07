#version 330 core

layout(location = 0) in vec4  Position;
layout(location = 1) in vec3  Normal;
layout(location = 2) in float PathCoord;

out float vPathCoord;

uniform vec3 modelPos;

void main(){
    gl_Position = vec4(Position.xy + modelPos.xy, 0.0, 1.0);
    vPathCoord  = PathCoord;
}
