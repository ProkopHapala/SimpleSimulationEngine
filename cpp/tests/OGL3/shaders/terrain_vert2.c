#version 330 core

// http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter02.html

// IN --- Input vertex data, different for all executions of this shader.
layout(location = 0) in vec2 vertPos_model;

// OUT --- Output data ; will be interpolated for each fragment.
noperspective out vec3 fragPos_world;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform mat4 camMat;
uniform float scale;
uniform vec2 terrain_tex0;

uniform sampler2D texRGB_vert;

void main(){

	vec2 xy_world = modelPos.xz + (vertPos_model*scale);

	vec2 size = vec2(10.0,10.0);
	float h = textureLod( texRGB_vert, (terrain_tex0+xy_world)/size, 0 ).r;

	vec3 position_world =  vec3( xy_world.x, modelPos.y+h,  -1 + xy_world.y*-1 );
	//vec3 position_world = modelPos + vec3( xy_world,  -5.0 );

	gl_Position         = camMat   * vec4( position_world, 1 );
	fragPos_world       = position_world;

	gl_Position.z   = -2.0*log(-position_world.z)-1.0;
}



