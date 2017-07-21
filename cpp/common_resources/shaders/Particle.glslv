#version 330 core

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 model_vpos;
layout(location = 1) in vec4 instance_pos;    // Position of the center of the particule and size of the square
layout(location = 2) in vec4 instance_color;  // Position of the center of the particule and size of the square

// Output data ; will be interpolated for each fragment.
//out vec2 UV;
out vec4 vcolor;

//uniform vec3 modelPos;
//uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;

void main(){
    //printf("(%f,%f,%f)\n",camPos.x,camPos.y,camPos.z);
    //const vec3 camPos = vec3(1.0,1.0,1.0);
    //const vec3 camPos = vec3(0.0,0.0,0.0);
	vec3 world_vpos = instance_pos.xyz + model_vpos * instance_pos.w;
	gl_Position     = camMat      * vec4( world_vpos-camPos, 1.0 );
	vcolor          = instance_color;
}

