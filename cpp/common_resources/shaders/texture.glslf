#version 330 core
// from : http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/

in      vec2      fUV;       // Interpolated values from the vertex shaders
uniform sampler2D texture_1;   // Values that stay constant for the whole mesh.

out vec4 gl_FragColor;

void main(){
	gl_FragColor   = textureLod( texture_1, fUV, 0 );
}


