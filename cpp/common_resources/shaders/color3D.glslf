#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in vec4 gl_FragCoord;
smooth in vec3 fragColor;

// OUT --- Ouput data
out vec3 color;

// UNI --- Values that stay constant for the whole mesh
//uniform vec3 light_dir;

void main(){
	color = fragColor;
	gl_FragDepth = gl_FragCoord.z;
}
