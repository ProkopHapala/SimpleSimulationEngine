#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in vec4 gl_FragCoord;
noperspective in vec3 world_pos;

// OUT --- Ouput data
out vec3 color;

// UNI --- Values that stay constant for the whole mesh.
// uniform vec4 baseColor;

uniform vec3 camPos;
//uniform mat4 camMat;

void main(){
	//color = baseColor.xyz;
	//color =  sin(world_pos*vec3(0.2,10.0,0.2))*0.5 + 0.5;
	
	vec3 ray   = world_pos - camPos;
	//color    = vec3(1.0,1.0,1.0)*sin(length(ray));
	float dist = length(ray);
	float sd   = sin(dist)*0.5+0.5;
	color      = vec3(0.05)*log(dist) + vec3(0.0,1.0,0.0)*(1/(1+sd*sd*1600.0));
	
	//gl_FragDepth = -dist;
	//gl_FragDepth = gl_FragCoord.z;
}
