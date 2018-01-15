#version 330 core

layout(location = 0) in vec2 vUV;

noperspective out vec3 world_pos;
//noperspective out vec3 world_nor;
smooth out vec3 fragColor;
//out float logz;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform vec3 camPos;
uniform mat4 camMat;

uniform vec2 uv0;
uniform vec2 pa;
uniform vec2 pb;
uniform vec3 mapScale;
uniform sampler2D txHeight;

void main(){
    vec2 p      = pa*(1.0-vUV.x) + pb*vUV.x;
    //p        += normalize(p)*vUV.y + modelPos.xz;
    p           = p*vUV.y + modelPos.xz;
    float h     = textureLod( txHeight, p*mapScale.xy + uv0, 0 ).r * mapScale.z;
    world_pos   = vec3( p.x ,h+modelPos.y, p.y );
    //world_pos   = vec3( p.x ,0.0, p.y );
    //world_pos = vec3( vUV.x, 0.0, vUV.y );
    gl_Position = camMat * vec4( world_pos-camPos, 1 );
}
