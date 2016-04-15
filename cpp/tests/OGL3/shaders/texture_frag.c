// from : http://www.opengl-tutorial.org/beginners-tutorials/tutorial-5-a-textured-cube/

#version 330 core

layout(location = 0) out vec3 color;

uniform vec2  resolution;
in vec2 UV;                   // Interpolated values from the vertex shaders
//out vec3 color;             // Ouput data
uniform sampler2D texture1;   // Values that stay constant for the whole mesh.

void main(){
   // color = texture( myTextureSampler, UV ).rgb;    // Output color = color of the texture at the specified UV
	//color = texture( myTextureSampler, (gl_FragCoord.xy/resolution.xy).xy ).rgb;
	
	vec2 uv = (gl_FragCoord.xy / resolution).xy;
	//color = texture( texture1, uv ).rgb;
	//float z = sin( 40000.0*texture( texture1, uv ).r )+1.0;
	//color = vec3(z,z,z);
	float z = texture( texture1, uv ).r;
	color  = sin( 10000.0* vec3(z,2*z,4*z)  )+1.0;
	
}


