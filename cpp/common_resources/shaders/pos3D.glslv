#version 330 core

//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertPos_model;

// OUT --- Output data ; will be interpolated for each fragment.
noperspective out vec3 world_pos;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;

void main(){
	world_pos    = modelPos + modelMat * vertPos_model;
	gl_Position  = camMat   * vec4( world_pos-camPos, 1 );
}



