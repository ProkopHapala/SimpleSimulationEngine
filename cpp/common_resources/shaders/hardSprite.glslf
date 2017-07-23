#version 330 core
// from : http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/

in      vec2      fUV;       // Interpolated values from the vertex shaders
uniform sampler2D texture_1;   // Values that stay constant for the whole mesh.

out vec4 gl_FragColor;

uniform vec4 keyColor;

void main(){
    vec3   clr   = textureLod( texture_1, fUV, 0 ).rgb;
    vec3   dclr  = clr - keyColor.rgb;
    float  err2  = dot(dclr,dclr)*keyColor.w; 
    if( err2<0.5 ){
        discard;
    }else{
        gl_FragColor = vec4( clr, clamp((err2-0.5)*2.0,0.0,1.0) ); 
        //gl_FragColor = vec4( err2, 0.0, 0.0, 1.0 ); 
        //gl_FragColor = vec4( dclr, 1.0 ); 
    }
}


