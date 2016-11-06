
#version 150

in  vec2 in_Position;

uniform float depth;

void main(void) {
	gl_Position = vec4(  in_Position,      depth,  1.0 );
	//gl_Position = vec4(  in_Position,      0.0,  1.0 );  // in front of everything
	//gl_Position = vec4(  in_Position,   0.9999,  1.0 );  // behind everything
}
