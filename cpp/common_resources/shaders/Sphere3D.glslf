#version 330 core

#define CUSTOM_DEPTH_0

// http://www.opengl-tutorial.org/beginners-tutorials/tutorial-8-basic-shading/
// https://www.gamedev.net/forums/topic/591110-geometry-shader-point-sprites-to-spheres/
// More info here: http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
// https://gamedev.stackexchange.com/questions/16588/computing-gl-fragdepth/16605#16605

//in        vec2 gl_PointCoord;
smooth in vec4 gl_FragCoord;
in        vec4 obj;
in        vec3 fpos_world;
out       vec4 gl_FragColor;

//uniform mat3 camRot;
uniform vec3 camPos;


const float inscribedRadius = 0.755761;     // https://en.wikipedia.org/wiki/Regular_icosahedron
const vec3  light_pos    = normalize( vec3( 1, 1, 2 ) );
const vec3  light_clr    = vec3( 1.0,  0.9,  0.8  );
const vec3  ambient_clr  = vec3( 0.15, 0.20, 0.25 );
const vec3  obj_clr      = vec3( 0.8,  0.8,  0.8  );
const vec2  obj_specular = vec2( 400.0, 0.5 );

// ==== functions

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

float lorenz(float x){ return 1/(1+x*x); }

/*
vec2 phong( vec3 nomral, vec3 hRay, vec3 light_dir ){ 
    // normal does not have to be normalized
    vec3  vrl = hRay+light_dir;
    vec3  iL2_vrl = 1/dot(vrl,vrl);
    vec3  iL2_nor = 1/dot(nomral, nomral);
    float D       = dot(normal, light_dir );
    float S       = dot(normal, vrl);
    return vec2(D*D,S*S*iL2_vrl)*L2_nor; 
}
*/

// ==== main

void main(){
	vec3 ray0 = camPos; 
	vec3 hRay = fpos_world - camPos; //hRay = normalize( hRay );
	float lZ = length(hRay); hRay *= (1/lZ); 
	//float lz2 = dot(hRay,hRay); hRay *= (1/sqrt(lz2));

	float t = raySphere( ray0, hRay, obj.xyz, obj.w*inscribedRadius );

	if( t > 1e+5 ){
       discard;
	}else{
		vec3 normal = sphereNormal( t, ray0, hRay, obj.xyz );   
		normal      = normalize( normal );
		float cD    = -dot( normal, light_pos );
		float cS    = -dot( normal, normalize(hRay+light_pos) ); 
		//float cS    = -dot( normal, normalize(hRay) ); 
		cS          = lorenz( (cS-1)*obj_specular.x )*obj_specular.y;
		gl_FragColor= vec4( (cD + cS)*obj_clr, 1.0 );
#ifdef CUSTOM_DEPTH_1
        //gl_FragDepth = gl_FragCoord.z;
        //gl_FragDepth = log(lZ)*0.1;
        gl_FragDepth = log(t)*0.1;
#endif 
	}
}

