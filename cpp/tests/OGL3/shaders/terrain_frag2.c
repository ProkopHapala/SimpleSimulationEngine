// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
// http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter08.html
// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch18.html

#version 330 core

smooth        in vec4 gl_FragCoord;

layout(location = 0) out vec3 color;

noperspective  in vec3 fragPos_world;

void main(){

	//float c = gl_FragCoord.y-10;
	//float c = fragPos_world.y+1.0;
	//color   = vec3(c,c,c);
	//color = fragPos_world*vec3(1.0,1.0,-0.01) + vec3(0.0,1.0,0.0);
	color = vec3( sin(100*fragPos_world.y), sin(10*fragPos_world.y), sin(fragPos_world.y) );

}


