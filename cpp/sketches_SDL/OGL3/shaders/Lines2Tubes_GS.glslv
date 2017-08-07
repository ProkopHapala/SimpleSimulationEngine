#version 330 core

// pass vert_data this way:
// https://gamedev.stackexchange.com/questions/49463/simple-pass-through-geometry-shader-with-normal-and-color


layout(lines) in;
layout(triangle_strip, max_vertices = 6) out;

//in  vData{ vec3 normal; vec4 color; }vertices[];
//out fData{ vec3 normal; vec4 color; }frag;

in  float vPathCoord[];
out float fPathCoord;

void main(){
    vec4 p0 = gl_in[0].gl_Position;
    vec4 p1 = gl_in[1].gl_Position;
 
    vec4 perp = vec4( vec3(1.0,1.0,0.0)*0.1  ,0.0 ); // we should figure this out later from camMat
    
    vec4 p00 = p0 - perp*0;
    vec4 p01 = p0 + perp;
    vec4 p10 = p1 - perp*0;
    vec4 p11 = p1 + perp;
    
    fPathCoord = vPathCoord[0];
    gl_Position = p00;  EmitVertex();
    gl_Position = p01;  EmitVertex();
    fPathCoord = vPathCoord[1];
    gl_Position = p10;  EmitVertex();
    gl_Position = p11;  EmitVertex();
    
}
