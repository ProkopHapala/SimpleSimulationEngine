#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

void main(){
	gl_FragColor   = textureLod( texture_1, fUV, 0 );
	//gl_FragColor   = vec4( fUV, sin(fUV.x)*sin(fUV.y), 1.0 );
}


