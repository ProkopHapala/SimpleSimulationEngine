#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in        vec4 gl_FragCoord;
noperspective in vec3 world_pos;

// OUT --- Ouput data
out vec3 color;

uniform vec3 camPos;

void main(){

	float iso = step( 0.9, fract(world_pos.y) ); 
	float h = world_pos.y*0.05;
	//vec3 d=world_pos-camPos;
	//float l  = ( log( dot(d,d) )-6.0 )*0.4;
	color = vec3(h,0.5-0.5*h,iso);
	//color = vec3(l*c);

	//gl_FragDepth = -dist;
	//gl_FragDepth = gl_FragCoord.z;
}
