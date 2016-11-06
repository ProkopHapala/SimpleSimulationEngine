// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch01.html
// http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter08.html
// http://http.developer.nvidia.com/GPUGems3/gpugems3_ch18.html

#version 330 core

layout(location = 0) out vec3 color;

noperspective  in vec3 fragPos_world;

uniform int   maxiter;
uniform float hrange; 
uniform float h0; 
uniform float dhmax;
uniform float dtmin;

uniform vec3      light_pos;
uniform vec3      cam_pos;
uniform vec2      size;
uniform vec2      tx0;
uniform sampler2D texRGB;


float height( vec2 p){
	return hrange*textureLod( texRGB, (p+tx0)/size, 0 ).r + h0;
	//return 0.5*(cos(p.x)*cos(p.y)+1.0) + h0;
}

float binSearch( vec3 hRay, vec3 p0, float tmax ){
	float dt  = tmax*0.5;
	float t   = dt;
	float err = 1e+8;      
	while(dt>dtmin*0.01){
		vec3 p = p0 + hRay * t;
		float h = height( p.xz );
		err = p.y - h;
		dt *= 0.5;
		if( err>0 ){ 
			t += dt;
		}else{  
			t -= dt;
		};		
	};
	return t;
}

void main(){

	vec3 hRay = fragPos_world - cam_pos;
	hRay      = normalize( hRay );
	float hmax = h0 + hrange;
	float t0  = ( hmax - cam_pos.y )/hRay.y;
	if (t0<0) discard;
	
	float t     = 0;
	float dt    = t0; 

	int   i     = 0;
	for(i; i<maxiter; i++){
		float t_  = t + dt;
		vec3  p   = cam_pos + t_ * hRay;
		float h   = height( p.xz );
		if( h > p.y ) break;
		dt  = (p.y-h)/(dhmax-hRay.y);
		//if( dt < dtmin ) break;
		dt  = max(dtmin, dt);
		t   = t_; 
	}
	if(i>=maxiter) discard;
	t = binSearch( hRay, cam_pos+t*hRay, dt )+t;
	vec3 p  = cam_pos + t * hRay;
	float c = (p.y - h0)/hrange;
	//float c = log(t*0.05)*0.5;
	color = vec3(c,c,c);

}


