#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// Interpolated values from the vertex shaders
in vec4 gl_FragCoord;
//in vec3 Normal_world;
in vec3 fragColor;

// Ouput data
out vec3 color;
//out float gl_FragDepth;

// Values that stay constant for the whole mesh.
uniform mat4 camMat_frag;
uniform vec3 light_dir;

//vec3 LightColor = vec3(1,1,1);
//float LightPower = 50.0f;
//uniform	vec3 MaterialDiffuseColor = texture( myTextureSampler, UV ).rgb;
//uniform	vec3 MaterialAmbientColor = vec3(0.1,0.1,0.1) * MaterialDiffuseColor;
//uniform	vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);

void main(){

	//co0lor = vec3(1.0,1.0,1.0);

	color = fragColor;
	//float c = gl_FragCoord.z;
	//color = vec3(c*fragColor,c*fragColor,c*fragColor);

	gl_FragDepth = gl_FragCoord.z;
	
	//color = 
	//	MaterialAmbientColor +
	//	MaterialDiffuseColor * LightColor * LightPower * cosTheta / (distance*distance) +
	//	MaterialSpecularColor * LightColor * LightPower * pow(cosAlpha,5) / (distance*distance);

}
