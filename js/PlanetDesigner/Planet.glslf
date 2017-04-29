
// For inspiration
// https://www.shadertoy.com/view/4s23WK
// https://www.shadertoy.com/view/lsl3RH
// https://www.shadertoy.com/view/lllXz4



/*
// Contemplating

 - According to schodeinger smoke we can create vorticity field by inversion of some complex "psi" function (zeros become poles = vortex cores )
 - http://multires.caltech.edu/pubs/SchrodingersSmoke.pdf   
 - the complex function "psi" can be found as fractal turbulence function  

*/




#define PI 3.14159275358979
#define POSITIVE_INF 1e+8
#define SKY_DIST     1e+7
#define PREC_STEP    1e-5



vec2 mul_complx( vec2 a, vec2 b ){
    return vec2( a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x );
}


/*
#define HASHSCALE1 .1031
#define HASHSCALE3 vec3(.1031, .1030, .0973)
#define HASHSCALE4 vec4(1031, .1030, .0973, .1099)

// Hash without Sine   https://www.shadertoy.com/view/4djSRW

float hash12(vec2 p){
	vec3 p3  = fract(vec3(p.xyx) * HASHSCALE1);
    p3      += dot(p3, p3.yzx + 19.19);
    return fract((p3.x + p3.y) * p3.z);
}

vec2 hash22(vec2 p){
	vec3 p3 = fract(vec3(p.xyx) * HASHSCALE3);
    p3     += dot(p3, p3.yzx + 19.19);
    return fract((p3.xx+p3.yz)*p3.zy);
}

vec4 hash42(vec2 p){
	vec4 p4 = fract(vec4(p.xyxy) * HASHSCALE4);
    p4     += dot(p4, p4.wzxy + 19.19);
    return fract((p4.xxyz+p4.yzzw)*p4.zywx);
}

vec3 hashSin33( vec3 p ){
	p = vec3( dot(p,vec3(127.1,311.7, 74.7)),
			  dot(p,vec3(269.5,183.3,246.1)),
			  dot(p,vec3(113.5,271.9,124.6)));
	return fract(sin(p)*43758.5453123);
}

vec2 hashSin22( vec2 p ){
	p = vec2( dot(p,vec2(127.1,311.7)),
			  dot(p,vec2(269.5,183.3)));
	return fract(sin(p)*43758.5453123);
}

vec4 surface_color( vec2 uv ){
    //return vec4( sin(uv*100.0), 1.0,1.0 );
    return vec4( vec3(hash12(uv)) ,1.0 );
}

*/


#define NUM_OCTAVES 8

const mat2 const_rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50) );

float random (in vec2 _st) { 
    return fract(sin(dot(_st.xy, vec2(12.9898,78.233))) * 43758.54531237);
}

// Based on Morgan McGuire @morgan3d
// https://www.shadertoy.com/view/4dS3Wd
float noise (in vec2 _st) {
    vec2 i = floor(_st);
    vec2 f = fract(_st);

    // Four corners in 2D of a tile
    float v00 = random(i                 );
    float v10 = random(i + vec2(1.0, 0.0));
    float v01 = random(i + vec2(0.0, 1.0));
    float v11 = random(i + vec2(1.0, 1.0));

    vec2 u  = f * f * (3.0 - 2.0 * f);
    vec2 mu = 1.0 - u; 

    return mu.y *( mu.x*v00 + u.x*v10 ) +
            u.y *( mu.x*v01 + u.x*v11 );
}

float fbm ( in vec2 _st) {
    float v = 0.0;
    float a = 0.5;
    vec2 shift = vec2(20.0);
    for (int i = 0; i < NUM_OCTAVES; ++i) {
        v    += a * noise(_st);
        _st   = const_rot * _st * 2.2 + shift;
        a    *= 0.5;
    }
    return v;
}

vec2 fbm_dist( in vec2 _st ){

    vec2 v      = _st;
    vec2 shift  = vec2(10.0);
    float a     = 0.3;
    for (int i  = 0; i < NUM_OCTAVES; ++i) {
        float phi = noise(_st) * 6.28 * 0.5;
        v      += a * vec2( cos(phi), sin(phi) );
        _st     = const_rot * _st * 1.6 + shift;
        a      *= 0.8;
    }

//T-URBULENCE_FUNCTION
    return v;
}

// Structs
struct Ray     { vec3 o;  vec3 d;  };
struct Sphere  { vec3  p; float R; };

// === Ray 

vec3 point( Ray r, float t ){ return r.o + t*r.d; }

// === Sphere

float dist2( Sphere sph, vec3 p ){
    vec3 dp = p - sph.p;
    return dot(dp,dp)-(sph.R*sph.R);
}

float ray_t( Sphere sph, Ray ray ){
    vec3 op   = sph.p - ray.o;
    float b   = dot(op, ray.d);
    float det = b*b - dot(op,op) + sph.R*sph.R;
    if (det<0.0) return POSITIVE_INF;
    det       = sqrt(det);
    float t   = b - det; 
    if (t < 0.0) t = b + det;
    return t;
}

vec3 normal( Sphere sph, vec3 p ){ return (p-sph.p)/sph.R; }

// === Main - RayTrace

uniform vec2  resolution;
uniform float time;
uniform mat3  camMat;

#define MAX_BOUNCES 2
float gamma = 2.2;

struct DirectionalLight{  vec3 d; vec3 c; };
struct Material{ vec3 color; float gloss; };

DirectionalLight sunLight = DirectionalLight( normalize(vec3(1.0, 0.5, 0.5)), vec3(1.0) );

float Lorenz( float x ){ return 1.0/(1.0+x*x); }
    
void main( ){	

    mat3 camMat_ = camMat;
	vec3 uvz = vec3( 2.0 * gl_FragCoord.xy / resolution.xy - 1.0, 5.0 );
	
	/*
    vec3 d  = normalize(vec3(resolution.x/resolution.y * uvz.x, uvz.y, -uvz.z ) );
    Ray ray = Ray(camMat_*vec3(0.0, 0.0, 10.0 ), camMat_*d);
    
    // === Scene  ( planet, rings?, moon? )
	
	Sphere planet = Sphere( vec3(0.0,0.0,0.0), 2.0 );
	float t = ray_t( planet, ray );
	vec3 p = ray.o + ray.d*t;
	
	vec2 uv = vec2( p.y/planet.R, atan(p.x,p.z) );
	
    if( t<SKY_DIST ){
        //gl_FragColor = vec4( sin(t*4.0),sin(t),sin(t*16.0), 1.0 );
        //gl_FragColor = vec4( sin(p*100.0), 1.0 );
        //gl_FragColor = surface_color(uv);

    }else{
        discard;
    }
    
    */
    
    //float c = hash12( hash22( uvz.xy + 15454.0 ) );
    //float c = hashSin22( uvz.xy ).x;
    //float c  = fbm( uvz.xy );
    
    
    //vec2  uv  = fbm_dist( uvz.xy );
    
    vec2  uv = uvz.xy;
    //uv += sin(uv*6.0 + 1.0);
    vec2  psi = sin( uv*10.0 + mul_complx(sin(uv*20.0),uv) ); 
    //gl_FragColor =  vec4( 0.01/(psi*psi), 0.5, 1.0 );
    
    float phi = 1.0/( length(psi) + 0.01 );
    vec2  d   = psi*phi;
    vec2  ur  = vec2( cos(phi), sin(phi) );
    
    d = mul_complx( d, ur );
    d = mul_complx( d, ur );
    d = mul_complx( d, ur );
    
    gl_FragColor =  vec4( sin((uv+d*0.1)*vec2(1.0,5.0)), 0.5, 1.0 );   
    
    //float c   = ( uvz.xy );
    
    //float c = hash12( uvz.xy + 15454.0 );
    //gl_FragColor =  vec4( c,c,c,1.0 ); 
}


