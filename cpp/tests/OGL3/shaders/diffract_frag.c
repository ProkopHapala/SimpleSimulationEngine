
#version 150

uniform int  nsource;
uniform int  window_size; 
uniform vec4 sources [100];

uniform float wave_length;
uniform float distance;
uniform float aperture; 
uniform float phase;
uniform float Iscale; 
uniform float scale; 


precision highp float;

//out vec4 gl_FragColor;

// ======== axuliary

int  half_size = window_size/2; 

float aperture2  = aperture*aperture; 
float freq       = 6.28318530718/wave_length;
float spot_size  = wave_length * distance / aperture;
float rho_scale2 = 1/(2*spot_size*spot_size) ; 

const vec2 ur = vec2(   1.0,  0.0 );  
const vec2 ug = vec2(  -0.5, -0.86602540378 );
const vec2 ub = vec2(  -0.5,  0.86602540378 );

void main(void) {

	vec2 pos     = ( gl_FragCoord.xy - vec2(half_size,half_size) ) * scale; 
	float  z0sq  = distance * distance;
	float iz0sq  = 1/z0sq;

	float rho2   = dot( pos, pos );
	float env   = exp( -rho2*rho_scale2 ); 

	int hit = -1;
	vec2 sum = vec2( 0.0, 0.0 );
	for( int i=0; i<nsource; i++ ){
		vec4 si    = sources[ i ];
		vec2 d     = si.xy - pos; 		
		float r2   = dot( d, d );
		if( r2 < aperture2 ){
			hit=1;
			break;
		}
		float a     = r2*iz0sq; 
		float dr    = distance*a*( 0.5 + a*( -0.125 + a*( 0.0625 + a* -0.0390625 ) ) ); 
		float phi   = freq * dr;
		vec2  eik   = vec2(  cos( phi ),  sin( phi ) );
		vec2 wave   = env*vec2(  eik.x*si.z - eik.y*si.w,    eik.y*si.z + eik.x*si.w );
		sum += wave;
	}

	if( hit>0 ){
		gl_FragColor = vec4( 0.5,0.5,0.5,1.0 );
	}else{
		sum*=Iscale;
		gl_FragColor = vec4( abs(sum.x), abs(sum.y), abs(sum.x-sum.y), 1.0 );
	}


}


