#version 330 core

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/

// IN --- Interpolated values from the vertex shaders
smooth in vec4 gl_FragCoord;
in vec2 gl_PointCoord;

// OUT --- Ouput data
out vec3 color;

// UNI --- Values that stay constant for the whole mesh.
//uniform vec4 baseColor;
uniform vec3 freq;

float lorentz(float f){ return 1.0/(1.0+f*f);};

void main(){
	//color = sin(gl_FragCoord.xyz*vec3(1.0,1.0,1000.0));
	
	vec2 uv = (gl_PointCoord-vec2(0.5,0.5) )*16.0;
	color   = vec3( 1.0-lorentz(uv.x), 1.0-2.0*lorentz(length(uv)-8.0), 1.0-lorentz(uv.y) );
	
	gl_FragDepth = gl_FragCoord.z;
}
