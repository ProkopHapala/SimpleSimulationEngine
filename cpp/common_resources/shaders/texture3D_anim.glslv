#version 330 core
layout(location = 0) in vec3 vpos;
layout(location = 1) in vec2 uv;

smooth out vec2 fUV1;
smooth out vec2 fUV2;
smooth out vec2 Cmix;

uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;

uniform float angle0;
const   float nPhases = 8;

void main(){
	vec3 vpos_world = modelPos + modelMat * vpos;
	vec3 vdir    = modelMat[0];
	float phase  = ( angle0 + atan ( vdir.x, vdir.y) ) * ( nPhases/6.28318530718 );
	float iphase = floor( phase );
	fUV1         = uv     + vec2( (iphase  )/nPhases,  1 );
	fUV2         = uv     + vec2( (iphase+1)/nPhases,  1 );
	Cmix.y       = phase - iphase;
	Cmix.x       = 1-Cmix.y;
	gl_Position  = camMat * vec4( vpos_world-camPos, 1 );
}



