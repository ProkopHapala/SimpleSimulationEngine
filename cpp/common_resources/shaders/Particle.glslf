#version 330 core

//in vec2 UV;
in vec4 vcolor;

out vec4 color;

//uniform sampler2D myTextureSampler;

void main(){
	// Output color = color of the texture at the specified UV
	// color = texture( myTextureSampler, UV ) * particlecolor;
	color = vcolor;
}
