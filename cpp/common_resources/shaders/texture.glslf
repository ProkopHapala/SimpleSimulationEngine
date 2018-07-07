#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

void main(){
	gl_FragColor   = textureLod( texture_1, fUV, 0 );
	//gl_FragColor   = vec4( fUV, sin(fUV.x*10.0)*sin(fUV.y*10.0), 1.0 );
	
	//gl_FragColor = vec4(fUV.x,fUV.y,0.0,1.0);
}


