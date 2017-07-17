#version 330 core

// IN --- Interpolated values from the vertex shaders
smooth in vec4 gl_FragCoord;
smooth in vec3 fragColor;
noperspective in vec3 world_pos;
in float logz;

// OUT --- Ouput data
out     vec3 color;

// UNI --- Values that stay constant for the whole mesh
uniform vec3 light_dir;
uniform vec3 camPos;

void main(){
	color = fragColor;
	
	// http://outerra.blogspot.cz/2012/11/maximizing-depth-buffer-range-and.html
	gl_FragDepth = logz;	
}
