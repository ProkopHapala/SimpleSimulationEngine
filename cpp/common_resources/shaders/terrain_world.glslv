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
//uniform vec3 terrainScale;

void main(){
    vec3 terrainScale = vec3( 0.02,0.01,40.0 );

    vec3 vert = vertPos_model;
    vert.y = sin( vertPos_model.x*terrainScale.x )*cos( vertPos_model.z*terrainScale.y ) * terrainScale.z; 
    
	world_pos       = modelPos + modelMat * vert;
	gl_Position     = camMat   * vec4( world_pos-camPos, 1 );
}



