
// Processing specific input
//uniform float time;
uniform vec2  resolution;

uniform vec4  sphere;
uniform vec3  light_dir;

//uniform float resolution;
//uniform vec2  mouse;

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

vec3 sphereNormal( float t, vec3 ray0, vec3 hRay, vec3 center ){
	return ray0 - center + t*hRay;
}


void main( void ){

	vec2 q  = (gl_FragCoord.xy/resolution.xy) - vec2(0.5,0.5);
	//vec2 q  = (gl_FragCoord.xy/resolution.xy) - (mouse.xy/resolution.xy);
    
	vec3 ray0 = vec3( 0.0, 0.0, 50.0 );	 
	vec3 hRay = ray0 + vec3( q.x, q.y, -40.0 );
	hRay = normalize( hRay ); 

	//vec3 pos = vec3( 0.0,0.0,0.0 );

	float t = raySphere( ray0, hRay, sphere.xyz, sphere.w );

	//float depth = -t;
	//float depth = (t+51.0);
	//float depth = 1.0 - 1.0/(100.0 + (t - sphere.z) );
	//float depth = 1.0 - 1.0/(100.0 + (t - sphere.z) );
	float depth = 1.0 - 1.0/(100.0 + t );
	//float depth = sphere.z;
	//float depth = 0.5;
	//float depth = clamp( 100.0/t, 0.0, 1.0 );

	if( t > 1e+5 ){
		gl_FragDepth = 0.9999999;
		gl_FragColor = vec4( 0.0, 0.0, 0.0, 1.0 );
		//discard;
	}else{
		vec3 normal = sphereNormal( t, ray0, hRay, sphere.xyz );
		normal = normalize( normal );
		float c_difuse   = -dot( normal, light_dir );
		//float c_specular = dot( normal, light_dir-hRay  ) * 0.1;  // c_specular = c_specular*c_specular; c_specular = c_specular*c_specular; c_specular = c_specular*c_specular;
		//gl_FragColor=vec4( clr_ambient + c_difuse*clr_diffuse + c_specular*clr_specular, 1.0 );
		gl_FragColor = vec4( clr_ambient + c_difuse*clr_diffuse, 1.0 );
		gl_FragDepth = depth;
		//gl_FragDepth = 2.0;
		//gl_FragDepth = depth;
		//gl_FragColor=vec4( c_specular*clr_specular, 1.0 );
	}

	//gl_FragColor = vec4( depth, depth, depth, 1.0 );
	//gl_FragDepth = depth;


	//gl_FragDepth = 1.0/(1.0+t);

	//gl_FragColor=vec4( 1.0, 1.0, 1.0, 1.0 );
	//gl_FragColor=vec4( sin(gl_FragCoord.xy/resolution), 1.0, 1.0 );

}


