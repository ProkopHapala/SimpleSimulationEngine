// from : http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/

#version 330 core

layout(location = 0) out vec3 color;

uniform vec2  resolution;
in vec2 UV;                 // Interpolated values from the vertex shaders
uniform sampler2D texRGB;   // Values that stay constant for the whole mesh.

void main(){
	
	vec2 uv = (gl_FragCoord.xy / resolution).xy;
	color   = textureLod( texRGB, uv, 0 ).rgb;

}


