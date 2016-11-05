#version 330 core

//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertPos_model;
layout(location = 1) in vec3 vertNormal_model;

// OUT --- Output data ; will be interpolated for each fragment.
noperspective out vec3 fragNormal_world;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform mat4 camMat;

void main(){
	//gl_Position =  MVC * vec4(vertexPosition_modelspace,1);
	vec3 position_world = modelPos + modelMat * vertPos_model;
	gl_Position         = camMat   * vec4( position_world, 1 );
	fragNormal_world    = modelMat * vertNormal_model;
}



