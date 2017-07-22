#version 330 core

// http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter02.html

// IN --- Input vertex data, different for all executions of this shader.
layout(location = 0) in vec2 vUV;

// OUT --- Output data ; will be interpolated for each fragment.
noperspective out vec3 world_pos;
//noperspective out vec3 world_nor;
smooth out vec3 fragColor;
//out float logz;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform vec3 camPos;
uniform mat4 camMat;

uniform vec2 uv_0;
uniform vec2 uv_da;
uniform vec2 uv_db;
uniform vec3 mapScale;
uniform sampler2D txHeight;

void main(){
    vec2 puv    = uv_da*vUV.x + uv_db*vUV.y + modelPos.xz;
    float h     = textureLod( txHeight, puv*mapScale.xy + uv_0, 0 ).r * mapScale.z;
    
    world_pos   = vec3( puv.x ,h+modelPos.y, puv.y );
    gl_Position = camMat * vec4( world_pos-camPos, 1 );
    
    //fragColor = sin( vec3( h*0.1,h, 10.0*h) );
    
    /*
    // http://outerra.blogspot.cz/2012/11/maximizing-depth-buffer-range-and.html
    const float far = 1000000.0;
    const float C   = 0.001;
    const float FC  = 1.0/log(far*C + 1);
    //logz = gl_Position.w*C + 1;  //version with fragment code
    logz = log(gl_Position.w*C + 1)*FC;
    gl_Position.z = (2*logz-1)*gl_Position.w;
    */

}
