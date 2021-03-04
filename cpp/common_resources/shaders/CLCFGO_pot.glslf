

uniform vec3  camPos;
uniform vec3  lookAt;
//uniform vec3  up;
uniform float zoom;

uniform vec2  resolution;

uniform int  natoms;
uniform vec4 atoms [100];
uniform vec4 coefs [100];

uniform float iso;
uniform vec3 color;
uniform vec3 light_dir;

//const float iso = 0.2;
//const vec3 light_dir = vec3(1.0,0.0,0.0);

//const float dt_default = 0.010;
const float dt_max   = 0.200;
const float dt_min   = 0.025;
const float dt_slope = 0.5;


//const float cutoff     = 6.0;
const float cutoff2    = 36.;
const float lor_cutoff = 0.02702702702; // 1/(1 + cutoff2 );

float LCAO( in vec3 pos ){
	float sum = 0.0;
	const float const_Coulomb = 14.3996;
	for(int i=0; i<natoms; i++){
		vec3 dr     = pos - atoms[i].xyz;
		float r2    = dot(dr,dr); 
		float R     = atoms[i].w;     // width is last element of atom
		float arg2  = r2/(R*R);
		if( arg2 < cutoff2 ){
			float radial = 1./( 1. + arg2 ) - lor_cutoff;   // Lorenzian
			sum += radial*( coefs[i].w + dot( dr, coefs[i].xyz ) )*const_Coulomb;
		}
	}
	return sum;
}

float Func( in vec3 pos ){
    float f = LCAO(pos);
    return iso-f*f;
    //return abs(LCAO(pos))-iso;
}


vec2 sphIntersect( in vec3 ro, in vec3 rd, in vec4 sph ){
	vec3 oc = ro - sph.xyz;
	float b = dot( oc, rd );
	float c = dot( oc, oc ) - sph.w*sph.w;
	float h = b*b - c;
	if( h<0.0 ) return vec2(-1.0,-1.0);
    float  sh =  sqrt( h );
	return vec2( -b-sh, -b+sh);
}

vec2 rot(in vec2 v,in vec2 u){ return vec2(v.x*u.x-v.y*u.y,v.x*u.y+v.y*u.x); }

float dirDeriv( vec3 pos, vec3 hat, float d ){ return (Func(pos+hat*d)-Func(pos))/d; }

vec3 grad( vec3 pos, float d ){ 
    float f0 = Func(pos);
    vec2 dv = vec2(0.0,d);
    return vec3( Func(pos+dv.yxx)-f0, Func(pos+dv.xyx)-f0, Func(pos+dv.xxy)-f0 )/d;
}

vec3 phong( vec3 rd, vec3 normal, vec3 color ){
	//normal = normalize(normal);
    float D = dot(light_dir,normal);
    float S = dot(normal,normalize(light_dir-rd) )-1.0;
    return color*( vec3(0.2,0.2,0.25)     // ambient
         + 0.8*clamp(D,0.0,1.0) )             // difuse
         + vec3(1.0/(S*S*25600.0+1.0));   // Specular
}

float curvatureShading( vec3 pos, vec3 normal, float d ){
    //return ( dirDeriv(pos,normal,d) - dirDeriv(pos,normal,-d) )/d;
    vec3 binor   = normalize(vec3(0.0,1.0,0.0));
    vec3 u       = normalize( cross(binor,normal) );
    vec3 v       = cross(u,normal);
    float f0     = Func(pos);
    float d2f_uu = Func(pos+u*d) + Func(pos+u*-d) - 2.0*f0;
    float d2f_vv = Func(pos+u*d) + Func(pos+u*-d) - 2.0*f0;
    return (d2f_uu + d2f_vv)/(d*d);
}

float linRoot( float t1, float dt, float f1, float f2 ){ return dt*f1/(f1-f2); }

float bisec( float t1, float t2, float f1, float f2, vec3 ro, vec3 rd ){ 
    for(int i=0; i<8; i++){
        float t   = t1 + (t2-t1)*f1/(f1-f2);
        float f   = Func(ro + t*rd);
        if( f<0.0 ){t2=t;f2=f;}
        else       {t1=t;f1=f;}
    }
    return t1 + (t2-t1)*f1/(f1-f2);
}

float optStep(float of, float f){
	float dt=dt_max;
	if( (of-f)>0.0 ){ dt = clamp(0.5*dt_max*f/(of-f),dt_min,dt_max); };
    return dt;
}

//==== Main

//void mainImage( out vec4 fragColor, in vec2 fragCoord ){
void main( void ){

    vec3 ro = camPos;
    vec3 rd = lookAt - camPos;
    vec3 lf = cross(rd,vec3(0.,1.,0.));  lf = normalize( lf );
    vec3 up = cross(rd,lf);              up = normalize( up );
    vec2 q  = ((gl_FragCoord.xy/resolution.xy) - vec2(0.5,0.5))*zoom;
    //hRay += lf*q.x + up*q.y; // perspective
    ro   += lf*q.x + up*q.y; // orthographic
    rd = normalize( rd );
    
    vec4 sph = vec4( lookAt,  4.0 );
    vec2 tt = sphIntersect( ro, rd, sph ); // bounding body
    
    vec4 clr = vec4(1.0,1.0,1.0,0.5);
    if( tt.x>0.0 ){
        //float t = tt.x;
        float t = 0.0;
        //float it = 0.0;
        float f,of;
        float dt = 0.2;
        //float dt = (tt.y-tt.x)/64.0;
        //for(int i=0; i<128; i++){ // overall Ray march
        for(int i=0; i<32; i++){ // overall Ray march
            f = Func(ro + t*rd);
            if(f<0.0)break;
            dt = optStep(of,f);
            of = f;
            t += dt; 
            //it+= (1.0/64.0)*5.0;
            if( (t>tt.y) ) break;
        }
        if(f<0.0){              // precise root and normal
            t    = bisec( t-dt, t, of, f, ro, rd );
            vec3 pos    = ro+rd*t;
            vec3 normal = normalize(grad( pos, 0.001));
            //vec3 color = vec3(1.0,1.0,1.0);
            //if( LCAO(pos)>0. ){ color = vec3(0.,0.5,1.); }
            //else              { color = vec3(1.,0.5,0.); }
            clr.rgb     = phong( rd, normal, color );
            //clr.rgb *= 1.5-(t-2.0)*0.5;
            //clr.rgb = color;
            float depth = 1.0 - 1.0/(100.0 + t );
            gl_FragDepth = depth;
            gl_FragColor = clr;
            return;
        }
    }
    discard;
    //gl_FragDepth = 0.9999999;
    //gl_FragColor = vec4( 1.0, 1.0, 1.0, 1.0 );
    //discard;
    //fragColor = clr;
    //gl_FragColor = clr;
}




