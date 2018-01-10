#version 330 core

smooth in vec3 fragColor;
uniform sampler3D texture_1; 

uniform vec3 txOffset;

out vec4 gl_FragColor;

void main(){
	//gl_FragColor   = textureLod( texture_1, fragColor, 0 );
	gl_FragColor   = texture( texture_1, fragColor+txOffset );
	//vec4 texelFetch( 	gsampler2DArray sampler, ivec3 P, int lod);
	//gl_FragColor   = vec4( fragColor, 1.0 );
	//gl_FragColor   = vec4( fUV, sin(fUV.x)*sin(fUV.y), 1.0 );
}





