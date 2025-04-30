#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

out vec4 fragColor; 

void main(){
    vec4 c  = textureLod( texture_1, fUV, 0 );
    c      += textureLod( texture_1, fUV, 1 );
    c      += textureLod( texture_1, fUV, 2 );
    c      += textureLod( texture_1, fUV, 3 );
	fragColor   = c*0.25;
	
	//gl_FragColor   = textureLod( texture_1, fUV, 32 );
	//gl_FragColor   = vec4( fUV, sin(fUV.x*10.0)*sin(fUV.y*10.0), 1.0 );
}


