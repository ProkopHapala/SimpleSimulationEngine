#version 330 core

layout(location = 0) in vec3  Position;
layout(location = 1) in vec3  Normal;
layout(location = 2) in float PathCoord;

out float vPathCoord;

uniform vec3 modelPos;
uniform mat3 modelMat;

void main(){
    gl_Position = vec4( modelPos + modelMat * Position, 1.0 );
    vPathCoord  = PathCoord;
}
