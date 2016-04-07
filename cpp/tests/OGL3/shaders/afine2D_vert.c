
#version 330

in  vec2 in_Position;

uniform mat2 afineMat;
uniform vec2 origin;

void main(void) {
	gl_Position = vec4(  (afineMat * in_Position) + origin,  0.0,  1.0 );
};
