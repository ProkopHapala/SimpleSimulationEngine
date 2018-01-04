#version 330 core

in      vec2      fUV;
uniform sampler2D texture_1; 

out vec4 gl_FragColor;

//uniform vec4 baseColor;
const vec3 baseColor = vec3(0.0,0.0,1.0);

//uniform dAAmin;
//uniform dAAmax;

const float dAAmin = 0.3;
const float dAAmax = 0.7;

void main(){
	vec4 mask = textureLod( texture_1, fUV, 0 );
	//gl_FragColor = vec4( mask.rrr, 1.0 );
	float d = mask.r;
	if( d<dAAmin ){
	    discard;
	}else{
	    gl_FragColor = vec4( baseColor, smoothstep( dAAmin, dAAmax, d ) );
	};
}


