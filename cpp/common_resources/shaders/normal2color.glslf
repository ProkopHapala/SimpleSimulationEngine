#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in vec4 gl_FragCoord;
smooth in vec3 fragColor;

// OUT --- Ouput data
out vec3 color;

void main(){
	color = fragColor*0.5 + vec3(0.5);
	gl_FragDepth = gl_FragCoord.z;
}
