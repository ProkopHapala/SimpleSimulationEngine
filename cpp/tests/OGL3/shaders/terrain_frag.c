// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
// http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter08.html
// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch18.html

#version 330 core

layout(location = 0) out vec3 color;


uniform float hmax; 
uniform float dhmax;

uniform vec2      light_pos;
uniform vec2      cam_pos;
uniform vec2      size;
uniform sampler2D texRGB;

void main(){
	vec3 hRay = normalize( gl_FragCoord - cam_pos );
	float t0 = ( cam_pos.y - hmax )/hRay.y;
	vec3  p = cam_pos + t0 * hRay;

	rxy = norm(hRay.xz);

	while( t<10000.0 ){
		//
		vec3 rgb = textureLod( texRGB, p.xz/size, 0 ).r;
		float h = textureLod( texRGB, p.xz/size, 0 ).r;
		if( h < p.y ){
			float dt = (p.y-h)/(hRay.z+dhmax);
			p += dt * hRay;
		}else{
			color = rgb; 
		}
	}

	//color   = textureLod( texRGB, uv, 0 ).rgb;
}


