
uniform vec2  resolution;
uniform vec3  light_dir;

uniform int  natoms;
uniform vec4 atoms [100];
uniform vec4 coefs [100];

const vec3 clr_diffuse  = vec3( 0.9, 0.8, 0.7 );
const vec3 clr_specular = vec3( 1.0, 1.0, 1.0 );
const vec3 clr_ambient  = vec3( 0.1, 0.2, 0.3 );

//const float cutoff     = 6.0;
const float cutoff2    = 16.;
const float exp_cutoff = 0.00247875217; // exp(-cutoff2 );

const float iso   = 0.2;
const float tmin  = 50.0;
const float tmax  = 60.0;
const float dtmin = 0.01;
const float dtmax = 0.1;

const float dderiv = 0.01;

vec2 rayPointDist( vec3 ray0, vec3 hRay, vec3 point ){
	vec3 pt  = point - ray0;	
	float t  = dot( hRay, pt );
    pt      -= t * hRay;
	float r2 = dot( pt, pt );  
	return vec2( t, r2 );
}

float evalFunc( vec3 pos ){
	float sum = 0.0;
	for(int i=0; i<natoms; i++){
		vec3 dr     = pos - atoms[i].xyz;
		float r2    = dot(dr,dr); 
		float R     = atoms[i].w;     // width is last element of atom
		float arg2  = r2/(R*R);
		if( arg2 < cutoff2 ){
			float radial = exp( -arg2 ) - exp_cutoff;   // Gaussians
			sum += radial*( coefs[i].w + dot( dr, coefs[i].xyz ) );
		}
	}
	return sum;
}

float bisec( vec3 ray0, vec3 hRay, float t, float dt ){
	while( dt > dtmin ){
		dt *=0.5;
		float thalf = t+dt;
		float val = evalFunc( ray0 + thalf*hRay );
		if( abs(val) < iso ) t = thalf;
	}
	return t;
}

float rayMarch( vec3 ray0, vec3 hRay ){
	float dt = dtmax;
	float t  = tmin;
	while( t<tmax ){
		float val = evalFunc( ray0 + t*hRay );
		if( abs(val) > iso )	break; 
		t += dt;
	}
	//return t;
	return bisec( ray0, hRay, t-dt, dt );
}

vec3 getNormal( vec3 hitpos, float val0 ){
	return vec3(
		evalFunc( hitpos + vec3( dderiv, 0.0, 0.0 ) ) - val0,
		evalFunc( hitpos + vec3( 0.0, dderiv, 0.0 ) ) - val0,
		evalFunc( hitpos + vec3( 0.0, 0.0, dderiv ) ) - val0
	) / dderiv;
}

float integrate( vec3 ray0, vec3 hRay ){
	float dt = dtmax;
	float t  = tmin;
	float sum = 0.0;
	while( t<tmax ){
		sum += evalFunc( ray0 + t*hRay );
		t   += dt;
	}
	return sum;
}


void main( void ){
	vec2 q    = (gl_FragCoord.xy/resolution.xy) - vec2(0.5,0.5);
	vec3 ray0 = vec3( 0.0, 0.0, -50.0 );	 
	vec3 hRay = vec3( q.x, q.y,  10.0 );
	hRay = normalize( hRay ); 

	float t = rayMarch( ray0, hRay );
	vec3  hitpos = ray0 + t*hRay;

	float  val0  = evalFunc ( hitpos ); 
	vec3  normal = getNormal( hitpos, val0 );
	normal       = normalize( normal ); if( val0<0.0 ) normal*=-1.0;
	float c_diff = dot ( light_dir,      normal );  
	//float c_spec = dot      ( light_dir-hRay, normal );  c_spec*=c_spec; c_spec*=c_spec; c_spec*=c_spec;
	//float c = c_diff + c_spec;
	float c = 0.1 + max( 0.0, c_diff );

	if( t < (tmax-dtmax) ){
		if( val0 > 0.0 ){
			gl_FragColor = vec4( 0.0, 0.5*c,   c, 1.0 );
		}else{
			gl_FragColor = vec4( c,   0.5*c, 0.0, 1.0 );
		};
		
	}else{
		gl_FragColor = vec4( 0.0, 0.0, 0.0, 1.0 );
	}

}


