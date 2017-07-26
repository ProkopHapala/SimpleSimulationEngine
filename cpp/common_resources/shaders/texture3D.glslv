#version 330 core
layout(location = 0) in vec3 vpos;
layout(location = 1) in vec2 uv;

smooth out vec2 fUV;

uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;

void main(){
	vec3 vpos_world = modelPos + modelMat * vpos;
	fUV = uv;
	gl_Position     = camMat   * vec4( vpos_world-camPos, 1 );
}



