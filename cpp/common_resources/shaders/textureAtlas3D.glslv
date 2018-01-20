#version 330 core
layout(location = 0) in vec3 vpos;
layout(location = 1) in vec2 uv;

uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;

uniform vec2  uv0;
uniform vec2  du;
uniform vec2  dv;

smooth out vec2 fUV;

void main(){
    vec3 vpos_world = modelPos + modelMat * vpos;
    vec3 vdir    = modelMat[0];

    fUV = uv0 + du*uv.x + dv*uv.y;
    gl_Position  = camMat * vec4( vpos_world-camPos, 1 );
}
