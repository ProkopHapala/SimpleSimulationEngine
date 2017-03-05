
// Processing specific input
//uniform float time;
#version 330 core

layout(location = 0) out vec3 color;

uniform vec2  resolution;

uniform mat3 camMat;
uniform vec3 camPos;

uniform vec4  sphere;
uniform vec3  light_dir;

uniform int     nholes;
uniform vec4[8]  holes;

//const vec3 light_dir = normalize( vec3( 1, 1, 2 ) );
const vec3 clr_diffuse  = vec3( 0.9, 0.8, 0.7 );
const vec3 clr_specular = vec3( 1.0, 1.0, 1.0 );
const vec3 clr_ambient  = vec3( 0.1, 0.2, 0.3 );

// More info here: http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

vec2 rayPointDist( vec3 ray0, vec3 hRay, vec3 point ){
	vec3 pt  = point - ray0;	
	float t  = dot( hRay, pt );
    pt      -= t * hRay;
	float r2 = dot( pt, pt );  
	return vec2( t, r2 );
}

float raySphere( vec3 ray0, vec3 hRay, vec3 center, float R ){
	vec2 res = rayPointDist( ray0, hRay, center );
	float dr2 = R*R - res.y;
	if( dr2 > 0.0 ){
		return res.x - sqrt( dr2 );
	}else{
		return 1e+8;
	}
}

bool checkHoles( vec3 hit ){ 
	for( int i=0; i<nholes; i++ ){
		vec3 dhit = hit - holes[i].xyz;
		float rh2 = dot(dhit,dhit);
		float r   = holes[i].w;
		if( rh2 < r*r  ) return true;
	}
	return false;
}

float raySphereHole( vec3 ray0, vec3 hRay, vec3 center, float R, out vec3 hit, out vec3 normal ){
	vec2 res = rayPointDist( ray0, hRay, center );
	float dr2 = R*R - res.y;
	if( dr2 > 0.0 ){
		float dr = sqrt( dr2 );
		float t  = res.x - dr;
		hit    = ray0  + t*hRay;
		normal = hit - center;
		if( checkHoles( normal ) ) {
			t = res.x + dr;
			hit  = ray0  + t*hRay;
			normal = hit - center;
			if( checkHoles( normal ) ) {
				return 1e+8;
			}
			normal *=-1;
		}
		return t;
	}
	return 1e+8;
}

vec3 sphereNormal( float t, vec3 ray0, vec3 hRay, vec3 center ){
	return ray0 - center + t*hRay;
}

void main( void ){
	vec2 q    = (gl_FragCoord.xy/resolution.xy) - vec2(0.5,0.5); 
	vec3 hRay = camMat * vec3( q.xy, 1.0 );
	hRay      = normalize( hRay ); 

	vec3 hit,normal;

	float t =raySphereHole( camPos, hRay, sphere.xyz, sphere.w, hit, normal );
	//float t   = raySphere( camPos, hRay, sphere.xyz, sphere.w );
	//normal    = sphereNormal( t, camPos, hRay, sphere.xyz );

	if( (t > 1e+7) || (t<0.0) ){
		gl_FragDepth = 0.9999999;
		color = vec3( 0.0, 0.0, 0.0 );
	}else{
		normal = normalize( normal );
		float c_difuse  = -dot( normal, light_dir );
		color = vec3( clr_ambient + c_difuse*clr_diffuse );

		float depth = 1.0 - 1.0/(100.0 + t );
		gl_FragDepth = depth;
	}

}


