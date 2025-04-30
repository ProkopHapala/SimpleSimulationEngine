#version 330 core

smooth in vec3 color1;
uniform sampler3D texture_1; 

uniform vec3 txOffset;

out vec4 fragColor;

void main(){
	//fragColor   = textureLod( texture_1, fragColor, 0 );
	fragColor   = texture( texture_1, color1+txOffset );
	//vec4 texelFetch( 	gsampler2DArray sampler, ivec3 P, int lod);
	//fragColor   = vec4( fragColor, 1.0 );
	//fragColor   = vec4( fUV, sin(fUV.x)*sin(fUV.y), 1.0 );
}





