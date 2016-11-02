#version 330 core

//http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertPos_model;
layout(location = 1) in vec3 vertColor;
//layout(location = 1) in vec3 vertNormal_model;

// Output data ; will be interpolated for each fragment.
out vec3 fragColor;
//out vec3 position_cam;
//out vec3 normal_world;

// Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform mat4 camMat;

void main(){
	//gl_Position =  MVC * vec4(vertexPosition_modelspace,1);
	vec3 position_world = modelPos + modelMat * vertPos_model;
	gl_Position  = camMat   * vec4( position_world, 1 );
	fragColor = vertColor;
	//position_cam        =  position_cam_.xyz;

	//vec3 n = normalize( Normal_cameraspace );
	//vec3 l = normalize( LightDirection_cameraspace );
}



