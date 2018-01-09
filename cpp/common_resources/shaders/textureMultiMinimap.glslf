#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

void main(){
    vec4 c  = textureLod( texture_1, fUV, 0 );
    c      += textureLod( texture_1, fUV, 1 );
    c      += textureLod( texture_1, fUV, 2 );
    c      += textureLod( texture_1, fUV, 3 );
	gl_FragColor   = c*0.25;
	
	//gl_FragColor   = textureLod( texture_1, fUV, 32 );
	//gl_FragColor   = vec4( fUV, sin(fUV.x)*sin(fUV.y), 1.0 );
}


