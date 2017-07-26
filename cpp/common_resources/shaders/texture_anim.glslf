#version 330 core

in      vec2      fUV1;
in      vec2      fUV2;
in      vec2      Cmix;
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

void main(){
	gl_FragColor   = Cmix.x*textureLod( texture_1, fUV1, 0 ) + Cmix.y*textureLod( texture_1, fUV2, 0 );
}


