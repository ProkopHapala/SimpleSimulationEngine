#version 330 core

// IN --- Input vertex data, different for all executions of this shader.
layout(location = 0) in vec3 vertPos_model;
layout(location = 1) in vec3 vertColor;

// OUT --- will be interpolated for each fragment.
smooth out vec3 fragColor;
noperspective out vec3 world_pos;
out float logz;

// UNI --- Values that stay constant for the whole mesh.
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;

//uniform float far = 1000000.0;
//uniform float C   = 0.001;

void main(){
	world_pos    = modelPos + modelMat * vertPos_model;
	gl_Position  = camMat   * vec4( world_pos-camPos, 1 );
	fragColor    = vertColor;
	
	// http://outerra.blogspot.cz/2012/11/maximizing-depth-buffer-range-and.html
    const float far = 1000000.0;
    const float C   = 0.001;
    const float FC  = 1.0/log(far*C + 1);
    //logz = gl_Position.w*C + 1;  //version with fragment code
    logz = log(gl_Position.w*C + 1)*FC;
    gl_Position.z = (2*logz-1)*gl_Position.w;
}



