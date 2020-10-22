#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in vec4 gl_FragCoord;

// OUT --- Ouput data
//out vec3 color;
out vec4 color;

// UNI --- Values that stay constant for the whole mesh.
uniform vec4 baseColor;

void main(){
	//color = baseColor.xyz;
	color = baseColor;
	gl_FragDepth = gl_FragCoord.z;
}
