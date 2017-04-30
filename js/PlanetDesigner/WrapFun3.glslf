
uniform vec2  resolution;
uniform float time;
uniform mat3  camMat;


//uniform float color_DEBUG;
uniform vec3 color_DEBUG;
uniform vec2 freq0;
uniform vec2 freqFac;
uniform vec2 scl; 
uniform vec2 amp;

vec3 lorenz( vec2 p ){
	float D  = 1.0/(1.0+dot(p,p));
	return vec3( p*(-2.0*D*D), D );
}

precision highp float;

//const float accod_coefs[4] = { 1.0, -1.0, 1.0, -1.0 };
//const float accod_coefs[4] = float[4]( 1.0, -1.0, 1.0, -1.0 );

vec3 accord_1( vec2 p ){
    p *= 6.0;
    vec2 sp  = sin(p); vec2 cp  = cos(p);
    vec2 spn = sp;     vec2 cpn = cp;
    vec3 fd = vec3(
        cpn.x*spn.y,
        spn.x*cpn.y,
        spn.x*spn.y
    );
    for(int i=1; i<4; i++){
        vec2 cpn_;
        cpn_ = cpn*cp - spn*sp;
        spn  = cpn*sp + spn*cp;
        cpn=cpn_;
        fd -= vec3(  
            cpn.x*spn.y,
            spn.x*cpn.y,
            spn.x*spn.y /float(i+1)
        );
    }
    fd.xy*=6.0;
    return fd;
    
}

void warp_sin( inout vec2 p, inout mat2 dT, vec2 frq ){
	p =p*scl;
	vec2 phi = p.yx*frq;
	vec2 f   =     amp*sin(phi); 
	vec2 df  = frq*amp*cos(phi);
	p       += f;
	dT       = dT*mat2( scl.x, df.x,  // jacobian transform
	                    df.y, scl.y   ); 
}

mat2 warping( inout vec2 p ){
    mat2 dT = mat2(1.0);
    vec2  freq = freq0;
    for( int i=0; i<__NWRAPS__; i++ ){ 
        //warp_sin2( p, dT, freq );
        warp_sin( p, dT, freq );
        freq *= freqFac;
    }
    return dT;
}

vec3 warped_function( vec2 p ){
    
	float fscale  = 12.0; 
	p *= fscale;	
	mat2 dT = warping( p );	
	vec3 fd = lorenz( p / fscale );
	fd.xy   = (dT * fd.xy);	
	
	//vec3 fd = lorenz( p*5.0  );
	//vec3 fd = accord_1( p );
	return fd;
}
    
void main( ){	
	vec2 p = 2.0 * gl_FragCoord.xy / resolution.xy - 1.0;
	vec3 o = warped_function( p );
    if( p.y > 0.0){
        float dd = 0.001;
        o.x = warped_function( p + vec2(dd,0.0) ).z - o.z;
        o.y = warped_function( p + vec2(0.0,dd) ).z - o.z;
        o.xy /= dd;
    };
    
    gl_FragColor =  __VIEW_FUNC__;
    //gl_FragColor =  vec4( vec3(o.z*0.5)+0.5, 1.0 );
    //gl_FragColor =  vec4( vec3(o.xy*0.02,0.0)+vec3(o.z*0.5)+0.2, 1.0 );
    //gl_FragColor =  vec4( o*vec3(0.05,0.05,0.5) +0.5, 1.0 );
    //gl_FragColor   =  vec4( o.yx*0.05+0.5,0.5, 1.0 ); 
    //gl_FragColor =  vec4( vec3(o.z)*0.5+0.5, 1.0 );   
    
    //gl_FragColor =  vec4( vec3(color_DEBUG), 1.0 );
    //gl_FragColor =  vec4( color_DEBUG, 1.0 );
}


