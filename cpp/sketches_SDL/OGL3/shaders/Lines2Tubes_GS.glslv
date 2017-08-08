#version 330 core

// pass vert_data this way:
// https://gamedev.stackexchange.com/questions/49463/simple-pass-through-geometry-shader-with-normal-and-color


layout(lines) in;
layout(triangle_strip, max_vertices = 6) out;

//in  vData{ vec3 normal; vec4 color; }vertices[];
//out fData{ vec3 normal; vec4 color; }frag;

in  float vPathCoord[];

flat   out vec3  p_start;
flat   out vec3  p_end;
smooth out vec3  fpos_world;
smooth out float fPathCoord;

uniform vec3 camPos;
uniform mat4 camMat;

void main(){
    vec3 p0 = gl_in[0].gl_Position.xyz;
    vec3 p1 = gl_in[1].gl_Position.xyz;
    
    vec3 camDir = camMat[3].xyz;
    vec3 cdir   = (p1-p0).xyz;
    vec3 perp   = normalize( cross( cdir, camDir ) ) * 0.12;
    //vec3 perp = vec3(1.0,1.0,0.0)*0.5; // we should figure this out later from camMat
    
    vec3 p00 = (p0 - perp) - cdir*0.05;
    vec3 p01 = (p0 + perp) - cdir*0.05;
    vec3 p10 = (p1 - perp) + cdir*0.05;
    vec3 p11 = (p1 + perp) + cdir*0.05;
    
    p_start = p0;
    p_end   = p1;
    
    fPathCoord  = vPathCoord[0];
    fpos_world  = p00; gl_Position = camMat*vec4(p00-camPos,1);  EmitVertex();
    fpos_world  = p01; gl_Position = camMat*vec4(p01-camPos,1);  EmitVertex();
    
    fPathCoord  = vPathCoord[1];
    fpos_world  = p10; gl_Position = camMat*vec4(p10-camPos,1);  EmitVertex();
    fpos_world  = p11; gl_Position = camMat*vec4(p11-camPos,1);  EmitVertex();
    
    /*
    gl_Position = camMat*vec4(p00-camPos,1);  EmitVertex();
    gl_Position = camMat*vec4(p01-camPos,1);  EmitVertex();
    gl_Position = camMat*vec4(p10-camPos,1);  EmitVertex();
    gl_Position = camMat*vec4(p11-camPos,1);  EmitVertex();
    */
    
    /*
    gl_Position = vec4(p00,1);  EmitVertex();
    gl_Position = vec4(p01,1);  EmitVertex();
    gl_Position = vec4(p10,1);  EmitVertex();
    gl_Position = vec4(p11,1);  EmitVertex();
    */

}
